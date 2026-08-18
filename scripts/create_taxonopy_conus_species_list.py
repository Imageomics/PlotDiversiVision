#!/usr/bin/env python3
"""Create a TaxonoPy-passed species list from a single CONUS state."""

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
            "Extract unique accepted plant names for one state from the CONUS "
            "plant list, run TaxonoPy, and write the final labels to assets."
        )
    )
    parser.add_argument(
        "--source-csv",
        type=Path,
        default=Path("assets/conus_plant_lists_accepted.csv"),
    )
    parser.add_argument("--state", required=True, help="State name, for example Colorado.")
    parser.add_argument(
        "--name",
        default=None,
        help="Output dataset stem. Defaults to <state>_conus in lowercase.",
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


def default_name(state: str) -> str:
    return f"{state.strip().lower().replace(' ', '_')}_conus"


def extract_species(args: argparse.Namespace, output_dir: Path) -> tuple[Path, int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    labels: set[str] = set()
    rows_seen = 0
    rows_used = 0
    target_state = args.state.casefold()

    with args.source_csv.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if clean(row.get("state")).casefold() != target_state:
                continue

            rows_seen += 1
            label = clean(row.get("Scientific.Name.with.Author"))
            if not label:
                continue
            labels.add(label)
            rows_used += 1

    if not labels:
        raise ValueError(f"No species labels found for state '{args.state}'")

    csv_path = output_dir / f"{args.name}_species_labels.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["label"])
        writer.writeheader()
        for label in sorted(labels):
            writer.writerow({"label": label})

    return csv_path, rows_seen, rows_used


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


def main() -> int:
    args = parse_args()
    args.name = args.name or default_name(args.state)
    if shutil.which(args.taxonopy_bin) is None:
        raise RuntimeError(
            f"Could not find '{args.taxonopy_bin}'. Install it with: "
            "pip install taxonopy"
        )

    output_dir = args.output_dir / args.name
    input_csv, rows_seen, rows_used = extract_species(args, output_dir)
    final_csv = run_taxonopy(args, input_csv, output_dir)

    with input_csv.open(newline="", encoding="utf-8") as f:
        unique_labels = sum(1 for _ in csv.DictReader(f))

    print(f"Rows seen for state: {rows_seen}")
    print(f"Rows with labels: {rows_used}")
    print(f"Unique labels: {unique_labels}")
    print(f"Intermediate CSV labels: {input_csv}")
    print(f"Final labels: {final_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
