#!/usr/bin/env python3
"""Predict plant species probabilities for grid crops using the pybioclip API.

This script takes one image, a directory of images, a glob of images, or the
repository's data/ image tree; splits each image into an N x N grid; and runs
each crop through BioCLIP 2 with a custom species list. Integration across
crops is intentionally left for a later step.

Examples:
    python scripts/predict_grid_species.py \
        --data-root data \
        --plot-id SCBI_008 \
        --output-csv outputs/grid_predictions/SCBI_008.csv \
        --grid-size 3

    python scripts/predict_grid_species.py \
        --images data/CPER_001 \
        --species-list assets/species_list/CPER_labels.csv \
        --output-csv outputs/grid_predictions/CPER_001.csv \
        --grid-size 4
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from pathlib import Path
from typing import Iterable

from PIL import Image


IMAGE_SUFFIXES = {".bmp", ".gif", ".jpeg", ".jpg", ".png", ".tif", ".tiff", ".webp"}
CROP_METADATA_FIELDS = [
    "image_path",
    "relative_image_path",
    "site_id",
    "year",
    "site_year",
    "plot_id",
    "subplot_id",
    "subplot_base",
    "subplot_row",
    "subplot_col",
    "image_date",
    "is_straightened",
    "crop_id",
    "grid_size",
    "grid_row",
    "grid_col",
    "x_min",
    "y_min",
    "x_max",
    "y_max",
]
PREDICTION_FIELDS = [
    *CROP_METADATA_FIELDS,
    "species_count",
    "probabilities",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split plant images into grid crops and predict species for each "
            "crop using BioCLIP 2 via the pybioclip Python API."
        )
    )
    parser.add_argument(
        "--images",
        nargs="+",
        default=[],
        help=(
            "Image path(s), directory path(s), or glob pattern(s). Quote globs "
            "so the script can expand them consistently. If omitted, images "
            "are discovered recursively from --data-root."
        ),
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data"),
        help="Root directory for the repository image tree.",
    )
    parser.add_argument(
        "--plot-id",
        action="append",
        default=[],
        help="Optional plot filter, for example SCBI_008 or CPER_001.",
    )
    parser.add_argument(
        "--species-list",
        type=Path,
        help=(
            "TaxonoPy-resolved species CSV. If omitted, the script uses "
            "assets/species_list/<site>_labels.csv when all matched images "
            "belong to one site."
        ),
    )
    parser.add_argument(
        "--species-column",
        default=None,
        help="CSV column containing species labels. Defaults to resolved_taxonomic_labels.",
    )
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument(
        "--grid-size",
        type=int,
        choices=(3, 4),
        default=3,
        help="Use 3 for a 3x3 grid or 4 for a 4x4 grid.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Number of crops to pass to BioCLIP per batch.",
    )
    parser.add_argument(
        "--device",
        default="cpu",
        help="Torch device for inference, for example 'cpu', 'cuda', or 'mps'.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover images and build crop metadata without running BioCLIP.",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def load_species(path: Path, species_column: str | None) -> list[str]:
    column = species_column or "resolved_taxonomic_labels"
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row")
        if column not in reader.fieldnames:
            raise ValueError(
                f"Column '{column}' was not found in {path}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        species = [clean(row.get(column)) for row in reader]

    unique_species = list(dict.fromkeys(label for label in species if label))
    if not unique_species:
        raise ValueError(f"No species labels found in {path}")
    return unique_species


def site_id_from_plot_id(plot_id: str) -> str:
    return plot_id.split("_", 1)[0] if "_" in plot_id else plot_id


def resolve_species_list_path(
    species_list: Path | None,
    plot_ids: list[str],
    image_paths: list[Path],
    data_root: Path,
) -> Path:
    if species_list is not None:
        return species_list

    site_ids = {site_id_from_plot_id(plot_id) for plot_id in plot_ids if plot_id}
    if not site_ids:
        site_ids = {
            metadata["site_id"]
            for metadata in (image_metadata(path, data_root) for path in image_paths)
            if metadata["site_id"]
        }

    if len(site_ids) != 1:
        raise ValueError(
            "Could not infer one site-level species list. Pass --species-list "
            "explicitly, or filter images to one site/plot."
        )

    site_id = next(iter(site_ids))
    inferred_path = Path("assets/species_list") / f"{site_id}_labels.csv"
    if not inferred_path.exists():
        raise FileNotFoundError(
            f"Inferred species list does not exist: {inferred_path}. "
            "Pass --species-list explicitly if you want a different label set."
        )
    return inferred_path


def image_metadata(path: Path, data_root: Path | None = None) -> dict[str, str]:
    normalized_stem = re.sub(r"\s+", "", path.stem)
    match = re.match(
        r"^(?P<plot_id>[A-Z]{4}_\d{3})_PlantDiversity_"
        r"(?P<subplot_base>\d+)_(?P<subplot_row>\d+)_(?P<subplot_col>\d+)_"
        r"(?P<image_date>\d{8})(?P<suffix>-straightened)?$",
        normalized_stem,
    )

    plot_id = match.group("plot_id") if match else ""
    site_id = plot_id.split("_", 1)[0] if plot_id else ""
    subplot_base = match.group("subplot_base") if match else ""
    subplot_row = match.group("subplot_row") if match else ""
    subplot_col = match.group("subplot_col") if match else ""
    subplot_id = (
        f"{subplot_base}_{subplot_row}_{subplot_col}"
        if subplot_base and subplot_row and subplot_col
        else ""
    )
    image_date = match.group("image_date") if match else ""
    is_straightened = "true" if match and match.group("suffix") else "false"

    site_year = ""
    year = ""
    parts = path.parts
    for part in parts:
        site_year_match = re.match(r"^(?P<site>[A-Z]{4})_(?P<year>\d{4})$", part)
        if site_year_match:
            site_year = part
            year = site_year_match.group("year")
            if not site_id:
                site_id = site_year_match.group("site")
            break

    if not year and image_date:
        year = image_date[:4]

    relative_image_path = str(path)
    if data_root:
        try:
            relative_image_path = str(path.resolve().relative_to(data_root.resolve()))
        except ValueError:
            relative_image_path = str(path)

    return {
        "relative_image_path": relative_image_path,
        "site_id": site_id,
        "year": year,
        "site_year": site_year,
        "plot_id": plot_id,
        "subplot_id": subplot_id,
        "subplot_base": subplot_base,
        "subplot_row": subplot_row,
        "subplot_col": subplot_col,
        "image_date": image_date,
        "is_straightened": is_straightened,
    }


def path_matches_filters(
    path: Path,
    data_root: Path,
    plot_ids: set[str],
) -> bool:
    metadata = image_metadata(path, data_root)
    if plot_ids and metadata["plot_id"] not in plot_ids:
        return False
    return True


def expand_images(
    inputs: Iterable[str],
    data_root: Path,
    plot_ids: set[str],
) -> list[Path]:
    image_paths: list[Path] = []
    input_items = list(inputs)

    if input_items:
        for item in input_items:
            matches = [Path(path) for path in glob.glob(item)]
            candidates = matches if matches else [Path(item)]

            for path in candidates:
                if path.is_dir():
                    image_paths.extend(
                        sorted(
                            child
                            for child in path.rglob("*")
                            if child.suffix.lower() in IMAGE_SUFFIXES
                        )
                    )
                elif path.is_file() and path.suffix.lower() in IMAGE_SUFFIXES:
                    image_paths.append(path)
    else:
        image_paths.extend(
            sorted(
                child
                for child in data_root.rglob("*")
                if child.suffix.lower() in IMAGE_SUFFIXES
            )
        )

    image_paths = [
        path
        for path in image_paths
        if path_matches_filters(path, data_root, plot_ids)
    ]

    unique_paths = list(dict.fromkeys(path.resolve() for path in image_paths))
    if not unique_paths:
        raise ValueError("No readable image files matched the image inputs/filters")
    return unique_paths


def grid_boxes(width: int, height: int, grid_size: int) -> list[tuple[int, int, int, int, int, int]]:
    boxes: list[tuple[int, int, int, int, int, int]] = []
    for row in range(grid_size):
        y_min = round(row * height / grid_size)
        y_max = round((row + 1) * height / grid_size)
        for col in range(grid_size):
            x_min = round(col * width / grid_size)
            x_max = round((col + 1) * width / grid_size)
            boxes.append((row, col, x_min, y_min, x_max, y_max))
    return boxes


def make_crops(
    image_paths: list[Path],
    grid_size: int,
    data_root: Path,
) -> tuple[list[Image.Image], list[dict[str, object]]]:
    crops: list[Image.Image] = []
    crop_meta: list[dict[str, object]] = []

    for image_path in image_paths:
        image_info = image_metadata(image_path, data_root)
        with Image.open(image_path) as image:
            image = image.convert("RGB")
            width, height = image.size
            for row, col, x_min, y_min, x_max, y_max in grid_boxes(
                width, height, grid_size
            ):
                crop = image.crop((x_min, y_min, x_max, y_max))
                crop_id = f"{image_path.stem}_r{row}_c{col}"

                crops.append(crop.copy())
                crop_meta.append(
                    {
                        "image_path": str(image_path),
                        **image_info,
                        "crop_id": crop_id,
                        "grid_size": grid_size,
                        "grid_row": row,
                        "grid_col": col,
                        "x_min": x_min,
                        "y_min": y_min,
                        "x_max": x_max,
                        "y_max": y_max,
                    }
                )

    return crops, crop_meta


def predict_crops(
    crops: list[Image.Image],
    crop_meta: list[dict[str, object]],
    species: list[str],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    import torch

    if sys.version_info >= (3, 12):
        torch.compile = lambda model, *_, **__: model

    from bioclip.predict import CustomLabelsClassifier

    classifier = CustomLabelsClassifier(cls_ary=species, device=args.device)
    predictions = classifier.predict(crops, k=len(species), batch_size=args.batch_size)

    probabilities_by_crop: dict[int, list[dict[str, object]]] = {
        index: [] for index in range(len(crops))
    }
    for prediction in predictions:
        crop_index = int(prediction["file_name"])
        probabilities_by_crop[crop_index].append(
            {
                "classification": prediction["classification"],
                "probability": prediction["score"],
            }
        )

    rows: list[dict[str, object]] = []
    for crop_index, meta in enumerate(crop_meta):
        probabilities = probabilities_by_crop[crop_index]
        rows.append(
            {
                **{field: meta[field] for field in CROP_METADATA_FIELDS},
                "species_count": len(species),
                "probabilities": json.dumps(probabilities, ensure_ascii=False),
            }
        )
    return rows


def write_predictions(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=PREDICTION_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    image_paths = expand_images(
        args.images,
        args.data_root,
        set(args.plot_id),
    )
    species_list = resolve_species_list_path(
        args.species_list,
        args.plot_id,
        image_paths,
        args.data_root,
    )
    species = load_species(species_list, args.species_column)
    crops, crop_meta = make_crops(image_paths, args.grid_size, args.data_root)
    if args.dry_run:
        print(f"Using species list: {species_list}")
        print(f"Loaded {len(species)} species labels")
        print(f"Matched {len(image_paths)} image(s)")
        print(f"Prepared {len(crops)} crop(s)")
        print("Dry run complete; BioCLIP was not loaded")
        return 0

    prediction_rows = predict_crops(crops, crop_meta, species, args)
    write_predictions(args.output_csv, prediction_rows)

    print(f"Loaded {len(species)} species labels")
    print(f"Using species list: {species_list}")
    print(f"Processed {len(image_paths)} image(s) into {len(crops)} crop(s)")
    print(f"Wrote predictions to {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
