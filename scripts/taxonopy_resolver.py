"""Shared TaxonoPy resolution helpers for species-list creation scripts."""

from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
RESOLVED_COLUMNS = (
    "resolution_status",
    "resolution_path",
    "resolution_strategy",
    "final_query_term",
    "final_query_rank",
    "final_data_source_id",
)


@dataclass
class TaxonopyConfig:
    input_csv: Path
    output_csv: Path
    work_dir: Path
    label_column: str
    source_dataset: str
    taxonopy_bin: str
    cache_dir: Path | None = None
    full_rerun: bool = False
    batch_size: int | None = None
    default_kingdom: str = "Plantae"
    env: dict[str, str] | None = None


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


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def stable_uuid(source_id: str, label: str) -> str:
    digest = hashlib.sha1(f"{source_id}\0{label}".encode("utf-8")).hexdigest()
    return (
        f"{digest[:8]}-{digest[8:12]}-{digest[12:16]}-"
        f"{digest[16:20]}-{digest[20:32]}"
    )


def rank_values(row: dict[str, str], label: str, default_kingdom: str) -> dict[str, str]:
    values = {rank: clean(row.get(rank)) for rank in RANKS}
    if not values["kingdom"]:
        values["kingdom"] = default_kingdom
    if not values["species"]:
        values["species"] = label
    if not values["genus"] and " " in label:
        values["genus"] = label.split()[0]
    return values


def prepare_taxonopy_input(config: TaxonopyConfig) -> tuple[Path, list[dict[str, str]]]:
    rows = read_csv(config.input_csv)
    prepared_by_key: dict[tuple[str, tuple[tuple[str, str], ...]], dict[str, str]] = {}

    for row_index, row in enumerate(rows, start=1):
        label = clean(row.get(config.label_column))
        if not label:
            continue

        ranks = rank_values(row, label, config.default_kingdom)
        key = (label, tuple(sorted(ranks.items())))
        if key in prepared_by_key:
            continue

        source_id = str(row_index)
        prepared_by_key[key] = {
            "uuid": stable_uuid(source_id, label),
            "scientific_name": label,
            "common_name": clean(row.get("common_name") or row.get("common")),
            "source_dataset": config.source_dataset,
            "source_id": source_id,
            "original_label": label,
            **ranks,
        }

    if not prepared_by_key:
        raise ValueError("No usable labels were found in the input CSV")

    taxonopy_input = config.work_dir / "taxonopy_input.csv"
    fieldnames = [
        "uuid",
        *RANKS,
        "scientific_name",
        "common_name",
        "source_dataset",
        "source_id",
    ]
    prepared_rows = list(prepared_by_key.values())
    write_csv(taxonopy_input, prepared_rows, fieldnames)
    return taxonopy_input, prepared_rows


def run_taxonopy(config: TaxonopyConfig, taxonopy_input: Path) -> Path:
    if shutil.which(config.taxonopy_bin) is None:
        raise RuntimeError(
            f"Could not find '{config.taxonopy_bin}'. Install it with: "
            "pip install taxonopy"
        )

    resolved_dir = config.work_dir / "taxonopy_resolved"
    cmd = [
        config.taxonopy_bin,
        "resolve",
        "--input",
        str(taxonopy_input),
        "--output-dir",
        str(resolved_dir),
        "--output-format",
        "csv",
    ]
    if config.cache_dir:
        cmd[1:1] = ["--cache-dir", str(config.cache_dir)]
    if config.full_rerun:
        cmd.append("--full-rerun")
    if config.batch_size:
        cmd.extend(["--batch-size", str(config.batch_size)])

    subprocess.run(cmd, check=True, env=config.env)
    return resolved_dir


def load_resolved_rows(resolved_dir: Path) -> list[dict[str, str]]:
    files = sorted(resolved_dir.glob("*.resolved.csv"))
    if not files:
        raise FileNotFoundError(f"No TaxonoPy CSV outputs found in {resolved_dir}")

    rows: list[dict[str, str]] = []
    for path in files:
        rows.extend(read_csv(path))
    return rows


def build_mapping(
    prepared_rows: list[dict[str, str]], resolved_rows: list[dict[str, str]]
) -> dict[str, dict[str, str]]:
    input_by_uuid = {row["uuid"]: row for row in prepared_rows}
    mapping: dict[str, dict[str, str]] = {}

    for resolved in resolved_rows:
        uuid = clean(resolved.get("uuid"))
        original = input_by_uuid.get(uuid, {})
        original_label = clean(original.get("original_label"))
        if not original_label:
            original_label = clean(resolved.get("scientific_name"))

        resolved_scientific = clean(resolved.get("scientific_name"))
        resolved_ranks = {rank: clean(resolved.get(rank)) for rank in RANKS}
        if resolved_ranks["species"]:
            resolved_ranks["species"] = resolved_ranks["species"].split()[-1]
        taxonomic_label = " ".join(
            value for value in (resolved_ranks[rank] for rank in RANKS) if value
        )

        mapping[original_label] = {
            "label": original_label,
            "resolved_labels": (
                " ".join([resolved_ranks["genus"], resolved_ranks["species"]]) or resolved_scientific
            ),
            "resolved_scientific_names": resolved_scientific,
            "resolved_taxonomic_labels": taxonomic_label,
            "taxonopy_resolution_statuses": clean(resolved.get("resolution_status")),
            **{f"resolved_{rank}": value for rank, value in resolved_ranks.items()},
            **{col: clean(resolved.get(col)) for col in RESOLVED_COLUMNS},
        }

    return mapping


def write_resolved_species_list(
    input_csv: Path,
    output_csv: Path,
    label_column: str,
    mapping: dict[str, dict[str, str]],
) -> None:
    output_rows: list[dict[str, str]] = []
    resolved_labels_set = set()
    for row in read_csv(input_csv):
        label = clean(row.get(label_column))
        resolved = mapping.get(label)

        # Uncomment the following lines if you want to force the inclusion of all resolved labels
        # if not resolved:
        #     continue

        if (
                not resolved or
                resolved["resolved_labels"] in resolved_labels_set or
                resolved["taxonopy_resolution_statuses"] == "FAILED_FORCED_INPUT" or
                len(resolved["resolved_taxonomic_labels"].split()) < 5
            ):
            continue
        resolved_labels_set.add(resolved["resolved_labels"])
        output_rows.append(
            {
                "label": label,
                "resolved_labels": resolved["resolved_labels"],
                "resolved_scientific_names": resolved["resolved_scientific_names"],
                "resolved_taxonomic_labels": resolved["resolved_taxonomic_labels"],
                "taxonopy_resolution_statuses": resolved[
                    "taxonopy_resolution_statuses"
                ],
            }
        )

    write_csv(
        output_csv,
        output_rows,
        [
            "label",
            "resolved_labels",
            "resolved_scientific_names",
            "resolved_taxonomic_labels",
            "taxonopy_resolution_statuses",
        ],
    )


def resolve_species_csv(config: TaxonopyConfig) -> Path:
    config.work_dir.mkdir(parents=True, exist_ok=True)
    resolved_dir = config.work_dir / "taxonopy_resolved"
    if resolved_dir.exists() and config.full_rerun:
        shutil.rmtree(resolved_dir)

    taxonopy_input, prepared_rows = prepare_taxonopy_input(config)
    resolved_dir = run_taxonopy(config, taxonopy_input)
    mapping = build_mapping(prepared_rows, load_resolved_rows(resolved_dir))
    write_resolved_species_list(
        config.input_csv,
        config.output_csv,
        config.label_column,
        mapping,
    )
    return resolved_dir
