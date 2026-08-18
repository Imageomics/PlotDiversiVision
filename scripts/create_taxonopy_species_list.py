#!/usr/bin/env python3
"""Create a TaxonoPy-passed species list from NEON plot data.

The workflow keeps TaxonoPy inputs and raw resolver output in ``outputs/`` and
writes only the final compact label file to ``assets/species_list/``.

Example:
    python scripts/create_taxonopy_species_list.py \
        --source-csv assets/NEON_plotData.csv \
        --plot-id SCBI_008 \
        --name SCBI_008 \
        --taxonopy-bin taxonopy
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path

from taxonopy_resolver import TaxonopyConfig, resolve_species_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract unique species labels for one or more NEON plots, run "
            "TaxonoPy, and place only the final compact label file under "
            "assets/species_list."
        )
    )
    parser.add_argument("--source-csv", required=True, type=Path)
    parser.add_argument("--name", required=True, help="Output dataset stem.")
    parser.add_argument(
        "--plot-id",
        action="append",
        default=[],
        help=(
            "NEON plotID value, for example SCBI_008. Use one plotID by "
            "default; repeated values create a combined list when needed."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/species_list"),
        help="Directory for extracted labels and TaxonoPy intermediate files.",
    )
    parser.add_argument(
        "--assets-dir",
        type=Path,
        default=Path("assets/species_list"),
        help="Directory for final compact TaxonoPy-passed label CSVs.",
    )
    parser.add_argument(
        "--resolved-species-list",
        type=Path,
        default=None,
        help=(
            "Optional existing TaxonoPy-passed species list. When provided, "
            "the extracted plot labels are mapped by lookup and TaxonoPy is "
            "not run again."
        ),
    )
    parser.add_argument("--taxonopy-bin", default="taxonopy")
    parser.add_argument(
        "--gnverifier-bin-dir",
        type=Path,
        default=Path("outputs/tools/gnverifier/bin"),
        help="Optional directory containing a local gnverifier executable.",
    )
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--full-rerun", action="store_true")
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def binomial_key(label: str) -> str:
    parts = label.split()
    if len(parts) < 2:
        return ""
    return f"{parts[0]} {parts[1]}"


def extract_species(args: argparse.Namespace, output_dir: Path) -> tuple[Path, int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    labels: set[str] = set()
    rows_seen = 0
    species_rows_used = 0

    plot_ids = set(args.plot_id)
    if not plot_ids:
        raise ValueError("Provide at least one --plot-id")

    with args.source_csv.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if plot_ids and clean(row.get("plotID")) not in plot_ids:
                continue

            rows_seen += 1
            if clean(row.get("taxonRank")) != "species":
                continue

            label = clean(row.get("scientificName"))
            if not label or label == "NA" or label.lower().startswith("unknown"):
                continue
            labels.add(label)
            species_rows_used += 1

    sorted_labels = sorted(labels)
    csv_path = output_dir / f"{args.name}_species_labels.csv"

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["label"])
        writer.writeheader()
        for label in sorted_labels:
            writer.writerow({"label": label})

    return csv_path, rows_seen, species_rows_used


def run_taxonopy(args: argparse.Namespace, input_csv: Path, output_dir: Path) -> Path:
    final_csv = args.assets_dir / f"{args.name}_labels.csv"
    work_dir = output_dir / "taxonopy"
    cache_dir = work_dir / "cache"
    env = os.environ.copy()
    if args.gnverifier_bin_dir.exists():
        env["PATH"] = f"{args.gnverifier_bin_dir.resolve()}{os.pathsep}{env['PATH']}"
    env["HOME"] = str((work_dir / "home").resolve())

    args.assets_dir.mkdir(parents=True, exist_ok=True)
    resolve_species_csv(
        TaxonopyConfig(
            input_csv=input_csv,
            output_csv=final_csv,
            work_dir=work_dir,
            label_column="label",
            source_dataset=args.name,
            taxonopy_bin=args.taxonopy_bin,
            cache_dir=cache_dir,
            full_rerun=args.full_rerun,
            batch_size=args.batch_size,
            env=env,
        )
    )
    return final_csv


def map_from_resolved_species_list(args: argparse.Namespace, input_csv: Path) -> Path:
    if args.resolved_species_list is None:
        raise ValueError("--resolved-species-list is required for lookup mapping")

    final_csv = args.assets_dir / f"{args.name}_labels.csv"
    with args.resolved_species_list.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{args.resolved_species_list} has no header row")
        lookup = {clean(row.get("label")): row for row in reader if clean(row.get("label"))}
    binomial_lookup = {
        binomial_key(label): row
        for label, row in lookup.items()
        if binomial_key(label)
    }

    output_rows: list[dict[str, str]] = []
    unmapped = 0
    with input_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = clean(row.get("label"))
            resolved = lookup.get(label) or binomial_lookup.get(binomial_key(label))
            if not resolved:
                unmapped += 1
                continue
            output_rows.append(
                {
                    "label": label,
                    "resolved_labels": clean(resolved.get("resolved_labels")),
                    "resolved_scientific_names": clean(
                        resolved.get("resolved_scientific_names")
                    ),
                    "resolved_taxonomic_labels": clean(
                        resolved.get("resolved_taxonomic_labels")
                    ),
                    "taxonopy_resolution_statuses": clean(
                        resolved.get("taxonopy_resolution_statuses")
                    ),
                }
            )

    args.assets_dir.mkdir(parents=True, exist_ok=True)
    with final_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "label",
                "resolved_labels",
                "resolved_scientific_names",
                "resolved_taxonomic_labels",
                "taxonopy_resolution_statuses",
            ],
        )
        writer.writeheader()
        writer.writerows(output_rows)

    if unmapped:
        print(f"Unmapped labels from existing resolved list: {unmapped}")
    return final_csv


def main() -> int:
    args = parse_args()
    if shutil.which(args.taxonopy_bin) is None:
        raise RuntimeError(
            f"Could not find '{args.taxonopy_bin}'. Install it with: "
            "pip install taxonopy"
        )

    output_dir = args.output_dir / args.name
    csv_path, rows_seen, species_rows_used = extract_species(args, output_dir)
    if args.resolved_species_list:
        final_csv = map_from_resolved_species_list(args, csv_path)
    else:
        final_csv = run_taxonopy(args, csv_path, output_dir)

    with csv_path.open(newline="", encoding="utf-8") as f:
        unique_labels = sum(1 for _ in csv.DictReader(f))

    print(f"Rows seen: {rows_seen}")
    print(f"Species rows used: {species_rows_used}")
    print(f"Unique labels: {unique_labels}")
    print(f"Intermediate CSV labels: {csv_path}")
    print(f"Final labels: {final_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
