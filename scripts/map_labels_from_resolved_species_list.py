#!/usr/bin/env python3
"""Create subplot-level test label files from resolved plot species lists."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


DEFAULT_SOURCE_CSV = Path("assets/NEON_plotData.csv")
DEFAULT_SPECIES_LIST_DIR = Path("assets/species_list")
DEFAULT_OUTPUT_DIR = Path("assets/test_labels")

RESOLVED_COLUMNS = [
    "resolved_labels",
    "resolved_scientific_names",
    "resolved_taxonomic_labels",
    "taxonopy_resolution_statuses",
]

OUTPUT_COLUMNS = [
    "plotID",
    "subplotID",
    "original_labels",
    *RESOLVED_COLUMNS,
    "label_count",
    "unmapped_original_labels",
]

SKIP_LABELS = {"", "NA", "N/A", "UNKNOWN", "UNKNOWN SP.", "UNKNOWN SPECIES"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a subplot-level test label file for one NEON plot by "
            "looking up its labels in an existing TaxonoPy-passed species list."
        )
    )
    parser.add_argument(
        "--source-csv",
        type=Path,
        default=DEFAULT_SOURCE_CSV,
        help=f"NEON plot data CSV. Default: {DEFAULT_SOURCE_CSV}",
    )
    parser.add_argument(
        "--plot-id",
        required=True,
        help="NEON plot ID to export, for example CPER_001 or SCBI_008.",
    )
    parser.add_argument(
        "--resolved-species-list",
        type=Path,
        help=(
            "TaxonoPy-passed species list. Default: "
            "assets/species_list/<plot-id>_labels.csv"
        ),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        help=(
            "Output subplot label file. Default: "
            "assets/test_labels/<plot-id>_subplot_labels.csv"
        ),
    )
    parser.add_argument(
        "--label-delimiter",
        default=";",
        help="Delimiter used for multi-label fields. Default: ';'",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header row")
        return [dict(row) for row in reader]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def default_species_list_path(plot_id: str) -> Path:
    return DEFAULT_SPECIES_LIST_DIR / f"{plot_id}_labels.csv"


def default_output_path(plot_id: str) -> Path:
    return DEFAULT_OUTPUT_DIR / f"{plot_id}_subplot_labels.csv"


def binomial_name(label: str) -> str:
    parts = label.split()
    if len(parts) < 2:
        return label
    return " ".join(parts[:2])


def natural_key(value: str) -> list[object]:
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def is_species_row(row: dict[str, str], plot_id: str) -> bool:
    if clean(row.get("plotID")) != plot_id:
        return False
    if clean(row.get("taxonRank")).lower() != "species":
        return False
    label = clean(row.get("scientificName"))
    return label.upper() not in SKIP_LABELS


def load_lookup(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    rows = read_csv(path)
    lookup: dict[str, dict[str, str]] = {}
    binomial_lookup: dict[str, dict[str, str]] = {}

    for row in rows:
        source_label = clean(row.get("label"))
        if not source_label:
            continue
        mapped = {column: clean(row.get(column)) for column in RESOLVED_COLUMNS}
        lookup[source_label] = mapped
        binomial_lookup.setdefault(binomial_name(source_label), mapped)

    if not lookup:
        raise ValueError(f"No usable labels found in {path}")
    return lookup, binomial_lookup


def map_label(
    label: str,
    lookup: dict[str, dict[str, str]],
    binomial_lookup: dict[str, dict[str, str]],
) -> dict[str, str] | None:
    if label in lookup:
        return lookup[label]
    return binomial_lookup.get(binomial_name(label))


def collect_subplot_labels(
    rows: list[dict[str, str]],
    plot_id: str,
) -> dict[str, set[str]]:
    subplot_labels: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if not is_species_row(row, plot_id):
            continue
        subplot_id = clean(row.get("subplotID"))
        label = clean(row.get("scientificName"))
        if subplot_id:
            subplot_labels[subplot_id].add(label)
    return subplot_labels


def make_output_rows(
    subplot_labels: dict[str, set[str]],
    plot_id: str,
    lookup: dict[str, dict[str, str]],
    binomial_lookup: dict[str, dict[str, str]],
    delimiter: str,
) -> tuple[list[dict[str, str]], int]:
    output_rows: list[dict[str, str]] = []
    unmapped_count = 0

    for subplot_id in sorted(subplot_labels, key=natural_key):
        labels = sorted(subplot_labels[subplot_id], key=natural_key)
        mapped: list[dict[str, str]] = []
        unmapped: list[str] = []

        for label in labels:
            mapped_label = map_label(label, lookup, binomial_lookup)
            if mapped_label is None:
                unmapped.append(label)
            else:
                mapped.append(mapped_label)

        unmapped_count += len(unmapped)
        output_row = {
            "plotID": plot_id,
            "subplotID": subplot_id,
            "original_labels": delimiter.join(labels),
            "label_count": str(len(labels)),
            "unmapped_original_labels": delimiter.join(unmapped),
        }
        for column in RESOLVED_COLUMNS:
            values = [item[column] for item in mapped if item[column]]
            output_row[column] = delimiter.join(values)
        output_rows.append(output_row)

    return output_rows, unmapped_count


def main() -> int:
    args = parse_args()
    species_list = args.resolved_species_list or default_species_list_path(args.plot_id)
    output_csv = args.output_csv or default_output_path(args.plot_id)

    source_rows = read_csv(args.source_csv)
    lookup, binomial_lookup = load_lookup(species_list)
    subplot_labels = collect_subplot_labels(source_rows, args.plot_id)
    if not subplot_labels:
        raise ValueError(f"No species-level subplot labels found for {args.plot_id}")

    rows, unmapped_count = make_output_rows(
        subplot_labels=subplot_labels,
        plot_id=args.plot_id,
        lookup=lookup,
        binomial_lookup=binomial_lookup,
        delimiter=args.label_delimiter,
    )
    write_csv(output_csv, rows)

    original_label_count = sum(len(labels) for labels in subplot_labels.values())
    print(f"Wrote test labels to {output_csv}")
    print(f"Plot: {args.plot_id}")
    print(f"Subplots: {len(rows)}")
    print(f"Original subplot labels: {original_label_count}")
    print(f"Unmapped labels: {unmapped_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
