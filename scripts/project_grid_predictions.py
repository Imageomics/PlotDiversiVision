#!/usr/bin/env python3
"""Project grid-level species predictions back onto their original images.

Examples:
    python scripts/project_grid_predictions.py

    python scripts/project_grid_predictions.py \
        --predictions temp_results/toy_run \
        --output-dir outputs/projected_predictions \
        --threshold 0.1
"""

from __future__ import annotations

import argparse
import csv
import textwrap
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


METADATA_COLUMNS = {
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
    "species_count",
}
REQUIRED_COLUMNS = {"image_path", "x_min", "y_min", "x_max", "y_max"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate original images with above-threshold grid predictions."
    )
    parser.add_argument(
        "--predictions",
        type=Path,
        default=Path("temp_results"),
        help="Prediction CSV file or directory to search recursively (default: temp_results).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("temp_results/projected_figures"),
        help="Directory for projected PNG figures (default: temp_results/projected_figures).",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data"),
        help="Fallback root for relative_image_path values (default: data).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Annotate species with a probability strictly above this value (default: 0.1).",
    )
    return parser.parse_args()


def prediction_csvs(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(path.rglob("*.csv"))
    raise FileNotFoundError(f"Prediction input does not exist: {path}")


def resolve_image_path(row: dict[str, str], data_root: Path) -> Path:
    image_path = Path(row["image_path"])
    if image_path.is_file():
        return image_path

    relative_path = row.get("relative_image_path", "")
    if relative_path:
        fallback_path = data_root / relative_path
        if fallback_path.is_file():
            return fallback_path

    raise FileNotFoundError(
        f"Could not find source image '{image_path}' or '{relative_path}' under {data_root}"
    )


def probability_columns(fieldnames: list[str] | None) -> list[str]:
    if fieldnames is None:
        raise ValueError("Prediction CSV has no header row")
    missing = REQUIRED_COLUMNS - set(fieldnames)
    if missing:
        raise ValueError(f"Prediction CSV is missing required columns: {', '.join(sorted(missing))}")
    return [column for column in fieldnames if column not in METADATA_COLUMNS]


def above_threshold_species(
    row: dict[str, str], columns: list[str], threshold: float
) -> list[str]:
    matches: list[tuple[str, float]] = []
    for species in columns:
        try:
            probability = float(row.get(species, ""))
        except ValueError:
            continue
        if probability > threshold:
            matches.append((species, probability))
    return [
        f"{species} ({probability:.3f})"
        for species, probability in sorted(matches, key=lambda item: item[1], reverse=True)
    ]


def load_font(size: int) -> ImageFont.ImageFont:
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size=size)
    except OSError:
        return ImageFont.load_default()


def annotation_lines(species: list[str], crop_width: int, crop_height: int) -> tuple[list[str], ImageFont.ImageFont]:
    wrapped = [line for name in species for line in textwrap.wrap(name, width=26)]
    font_size = min(24, max(10, crop_width // 26), max(10, crop_height // (len(wrapped) + 2)))
    return wrapped, load_font(font_size)


def annotate_image(image_path: Path, rows: list[dict[str, str]], columns: list[str], threshold: float) -> Image.Image:
    with Image.open(image_path) as source:
        image = source.convert("RGB")

    draw = ImageDraw.Draw(image, "RGBA")
    outline_width = max(2, min(image.size) // 500)
    for row in rows:
        x_min, y_min, x_max, y_max = (int(float(row[field])) for field in ("x_min", "y_min", "x_max", "y_max"))
        draw.rectangle((x_min, y_min, x_max, y_max), outline=(255, 215, 0, 255), width=outline_width)

        species = above_threshold_species(row, columns, threshold)
        if not species:
            continue
        lines, font = annotation_lines(species, x_max - x_min, y_max - y_min)
        text = "\n".join(lines)
        left, top = x_min + outline_width + 3, y_min + outline_width + 3
        bbox = draw.multiline_textbbox((left, top), text, font=font, spacing=2)
        draw.rounded_rectangle((bbox[0] - 3, bbox[1] - 3, bbox[2] + 3, bbox[3] + 3), radius=2, fill=(0, 0, 0, 185))
        draw.multiline_text((left, top), text, font=font, fill=(255, 255, 255, 255), spacing=2)
    return image


def output_path(csv_path: Path, predictions: Path, image_path: Path, output_dir: Path) -> Path:
    csv_relative = csv_path.relative_to(predictions) if predictions.is_dir() else Path(csv_path.stem)
    return output_dir / csv_relative.with_suffix("") / f"{image_path.stem}.png"


def project_csv(csv_path: Path, predictions: Path, output_dir: Path, data_root: Path, threshold: float) -> int:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = probability_columns(reader.fieldnames)
        grouped_rows: dict[Path, list[dict[str, str]]] = defaultdict(list)
        for row in reader:
            grouped_rows[resolve_image_path(row, data_root)].append(row)

    for image_path, rows in grouped_rows.items():
        destination = output_path(csv_path, predictions, image_path, output_dir)
        destination.parent.mkdir(parents=True, exist_ok=True)
        annotate_image(image_path, rows, columns, threshold).save(destination)
    return len(grouped_rows)


def main() -> None:
    args = parse_args()
    if not 0 <= args.threshold <= 1:
        raise ValueError("--threshold must be between 0 and 1")

    csv_paths = prediction_csvs(args.predictions)
    if not csv_paths:
        raise FileNotFoundError(f"No CSV files found in {args.predictions}")

    figure_count = 0
    for csv_path in csv_paths:
        figure_count += project_csv(
            csv_path, args.predictions, args.output_dir, args.data_root, args.threshold
        )
    print(f"Wrote {figure_count} projected figures to {args.output_dir}")


if __name__ == "__main__":
    main()
