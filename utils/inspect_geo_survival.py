#!/usr/bin/env python3
"""Inspect GEO series matrix files for sample metadata and survival fields."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


SURVIVAL_KEYS = (
    "overall survival",
    "overall_survival",
    "overall.survival",
    "survival status",
    "survival_status",
    "vital status",
    "death",
    "code_os",
    "disease-free survival",
    "progression-free survival",
    "pfs",
    "dfs",
)


def clean(value: str) -> str:
    return value.strip().strip('"')


def parse_matrix(path: Path) -> dict[str, object]:
    sample_count = 0
    table_rows = 0
    table_cols = 0
    survival_fields: set[str] = set()
    titles: list[str] = []
    geo: list[str] = []
    sources: list[str] = []

    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [clean(v) for v in next(csv.reader([line], delimiter="\t"))[1:]]
                sample_count = len(titles)
            elif line.startswith("!Sample_geo_accession"):
                geo = [clean(v) for v in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                sources = [clean(v) for v in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_characteristics"):
                values = [clean(v) for v in next(csv.reader([line], delimiter="\t"))[1:]]
                for value in values:
                    key = value.split(":", 1)[0].strip().lower() if ":" in value else value.lower()
                    if any(s in key for s in SURVIVAL_KEYS):
                        survival_fields.add(key)
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_rows += 1
                if table_rows == 1:
                    table_cols = len(line.rstrip("\n").split("\t"))

    tumor_like = sum(1 for value in sources + titles if re.search(r"tumou?r|cancer|ALL|NET|GBM|glioblastoma", value, re.I))
    return {
        "file": str(path),
        "samples": sample_count,
        "table_rows": table_rows,
        "table_cols": table_cols,
        "tumor_like_labels": tumor_like,
        "survival_fields": "; ".join(sorted(survival_fields)),
        "has_survival": bool(survival_fields),
        "has_data_table": table_rows > 1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("k_15/external/geo_matrix_inspection.csv"))
    args = parser.parse_args()

    rows = [parse_matrix(path) for path in args.paths]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            Path(row["file"]).name,
            f"samples={row['samples']}",
            f"data_rows={row['table_rows']}",
            f"survival={row['has_survival']}",
        )


if __name__ == "__main__":
    main()
