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
import hashlib
import glob
import json
import re
import sys
from pathlib import Path
from typing import Iterable

from PIL import Image


IMAGE_SUFFIXES = {".bmp", ".gif", ".jpeg", ".jpg", ".png", ".tif", ".tiff", ".webp"}
DEFAULT_SPECIES_COLUMN = "resolved_taxonomic_labels"
DEFAULT_PROBABILITY_COLUMN_NAME = "resolved_labels"
DEFAULT_TEXT_EMBEDDING_CACHE_DIR = Path("outputs/text_embeddings")
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
        help=f"CSV column containing species labels. Defaults to {DEFAULT_SPECIES_COLUMN}.",
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
        "--text-embedding-cache-dir",
        type=Path,
        default=DEFAULT_TEXT_EMBEDDING_CACHE_DIR,
        help=(
            "Directory for cached BioCLIP text embeddings. Default: "
            f"{DEFAULT_TEXT_EMBEDDING_CACHE_DIR}"
        ),
    )
    parser.add_argument(
        "--no-text-embedding-cache",
        action="store_true",
        help="Recompute text embeddings instead of reading/writing the cache.",
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


def effective_species_column(species_column: str | None) -> str:
    return species_column or DEFAULT_SPECIES_COLUMN


def make_unique_column_names(
    names: list[str],
    reserved_names: set[str] | None = None,
) -> list[str]:
    reserved_names = set(reserved_names or set())
    used = set(reserved_names)
    unique_names: list[str] = []
    counts: dict[str, int] = {}

    for name in names:
        base_name = name or "species_probability"
        candidate = base_name
        while candidate in used:
            counts[base_name] = counts.get(base_name, 1) + 1
            candidate = f"{base_name}_{counts[base_name]}"
        used.add(candidate)
        unique_names.append(candidate)

    return unique_names


def load_species(path: Path, species_column: str | None) -> tuple[list[str], list[str]]:
    column = effective_species_column(species_column)
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row")
        if column not in reader.fieldnames:
            raise ValueError(
                f"Column '{column}' was not found in {path}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        probability_name_column = (
            DEFAULT_PROBABILITY_COLUMN_NAME
            if DEFAULT_PROBABILITY_COLUMN_NAME in reader.fieldnames
            else column
        )
        species_rows = [
            (
                clean(row.get(column)),
                clean(row.get(probability_name_column)) or clean(row.get(column)),
            )
            for row in reader
        ]

    unique_species_by_label: dict[str, str] = {}
    for label, probability_column in species_rows:
        if label and label not in unique_species_by_label:
            unique_species_by_label[label] = probability_column

    unique_species = list(unique_species_by_label.keys())
    if not unique_species:
        raise ValueError(f"No species labels found in {path}")
    probability_columns = [
        unique_species_by_label[label]
        for label in unique_species
    ]
    probability_columns = make_unique_column_names(
        probability_columns,
        reserved_names=set(PREDICTION_FIELDS),
    )
    return unique_species, probability_columns


def text_embedding_cache_metadata(
    species: list[str],
    species_list: Path,
    species_column: str,
    model_str: str,
) -> dict[str, object]:
    return {
        "version": 1,
        "model_str": model_str,
        "species_list": str(species_list.resolve()),
        "species_column": species_column,
        "species_count": len(species),
        "species": species,
    }


def text_embedding_cache_path(
    cache_dir: Path,
    metadata: dict[str, object],
) -> Path:
    digest = hashlib.sha256(
        json.dumps(metadata, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    species_list_stem = Path(str(metadata["species_list"])).stem
    species_column = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(metadata["species_column"]))
    return cache_dir / f"{species_list_stem}_{species_column}_{digest}.pt"


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
    species_probability_columns: list[str],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    import torch

    if sys.version_info >= (3, 12):
        torch.compile = lambda model, *_, **__: model

    from bioclip.predict import CustomLabelsClassifier

    class CachedCustomLabelsClassifier(CustomLabelsClassifier):
        def __init__(
            self,
            cls_ary: list[str],
            cache_dir: Path | None,
            species_list: Path,
            species_column: str,
            **kwargs: object,
        ) -> None:
            self.text_embedding_cache_dir = cache_dir
            self.species_list = species_list
            self.species_column = species_column
            self.text_embedding_cache_status = "disabled"
            self.text_embedding_cache_path = None
            super().__init__(cls_ary=cls_ary, **kwargs)

        def _get_txt_embeddings(self, classnames: list[str]):  # type: ignore[override]
            if self.text_embedding_cache_dir is None:
                return super()._get_txt_embeddings(classnames)

            metadata = text_embedding_cache_metadata(
                species=list(classnames),
                species_list=self.species_list,
                species_column=self.species_column,
                model_str=getattr(self, "model_str", "unknown"),
            )
            cache_path = text_embedding_cache_path(
                self.text_embedding_cache_dir,
                metadata,
            )
            self.text_embedding_cache_path = cache_path

            if cache_path.exists():
                cached = torch.load(cache_path, map_location=self.device)
                if (
                    isinstance(cached, dict)
                    and cached.get("metadata") == metadata
                    and "txt_embeddings" in cached
                ):
                    txt_embeddings = cached["txt_embeddings"].to(self.device)
                    if txt_embeddings.shape[1] == len(classnames):
                        self.text_embedding_cache_status = "hit"
                        return txt_embeddings

            txt_embeddings = super()._get_txt_embeddings(classnames)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            torch.save(
                {
                    "metadata": metadata,
                    "txt_embeddings": txt_embeddings.detach().cpu(),
                },
                cache_path,
            )
            self.text_embedding_cache_status = "miss"
            return txt_embeddings

    cache_dir = None if args.no_text_embedding_cache else args.text_embedding_cache_dir
    classifier = CachedCustomLabelsClassifier(
        cls_ary=species,
        cache_dir=cache_dir,
        species_list=args.resolved_species_list_path,
        species_column=args.resolved_species_column,
        device=args.device,
    )
    if classifier.text_embedding_cache_path is not None:
        print(
            "Text embedding cache: "
            f"{classifier.text_embedding_cache_status} "
            f"({classifier.text_embedding_cache_path})"
        )
    predictions = classifier.predict(crops, k=len(species), batch_size=args.batch_size)

    species_index = {label: index for index, label in enumerate(species)}
    probabilities_by_crop: dict[int, list[float]] = {
        index: [0.0] * len(species) for index in range(len(crops))
    }
    for prediction in predictions:
        crop_index = int(prediction["file_name"])
        species_position = species_index[prediction["classification"]]
        probabilities_by_crop[crop_index][species_position] = prediction["score"]

    rows: list[dict[str, object]] = []
    for crop_index, meta in enumerate(crop_meta):
        probability_vector = probabilities_by_crop[crop_index]
        rows.append(
            {
                **{field: meta[field] for field in CROP_METADATA_FIELDS},
                "species_count": len(species),
                **dict(zip(species_probability_columns, probability_vector)),
            }
        )
    return rows


def write_predictions(
    path: Path,
    rows: list[dict[str, object]],
    species_probability_columns: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[*PREDICTION_FIELDS, *species_probability_columns],
        )
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
    species, species_probability_columns = load_species(
        species_list,
        args.species_column,
    )
    args.resolved_species_list_path = species_list
    args.resolved_species_column = effective_species_column(args.species_column)
    crops, crop_meta = make_crops(image_paths, args.grid_size, args.data_root)
    if args.dry_run:
        print(f"Using species list: {species_list}")
        print(f"Loaded {len(species)} species labels")
        print(f"Matched {len(image_paths)} image(s)")
        print(f"Prepared {len(crops)} crop(s)")
        print("Dry run complete; BioCLIP was not loaded")
        return 0

    prediction_rows = predict_crops(
        crops,
        crop_meta,
        species,
        species_probability_columns,
        args,
    )
    write_predictions(args.output_csv, prediction_rows, species_probability_columns)

    print(f"Loaded {len(species)} species labels")
    print(f"Using species list: {species_list}")
    print(f"Processed {len(image_paths)} image(s) into {len(crops)} crop(s)")
    print(f"Wrote predictions to {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
