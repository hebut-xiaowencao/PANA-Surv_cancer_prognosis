#!/usr/bin/env python3
"""Prepare UCSC Xena HiSeqV2 matrices for downstream DEG analysis.

The TCGA Xena HiSeqV2 matrices are log2-scale expression matrices. This script
back-transforms values with 2^x - 1, aligns columns to local omics sample order,
and writes a gene x sample CSV for each cancer.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def col_to_index(cell_ref: str) -> int:
    letters = re.sub(r"\d", "", cell_ref)
    idx = 0
    for char in letters:
        idx = idx * 26 + ord(char.upper()) - ord("A") + 1
    return idx - 1


def read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    try:
        with zf.open("xl/sharedStrings.xml") as handle:
            root = ET.parse(handle).getroot()
    except KeyError:
        return []
    strings: list[str] = []
    for si in root.iter(f"{NS_MAIN}si"):
        strings.append("".join(node.text or "" for node in si.iter(f"{NS_MAIN}t")))
    return strings


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{NS_MAIN}t"))
    value_node = cell.find(f"{NS_MAIN}v")
    if value_node is None or value_node.text is None:
        return ""
    value = value_node.text
    if cell_type == "s":
        return shared_strings[int(value)]
    return value


def row_values(row: ET.Element, shared_strings: list[str], width: int | None = None) -> list[str]:
    values = [""] * (width or 0)
    for cell in row.findall(f"{NS_MAIN}c"):
        idx = col_to_index(cell.attrib["r"])
        if idx >= len(values):
            values.extend([""] * (idx + 1 - len(values)))
        values[idx] = cell_value(cell, shared_strings)
    return values


def read_omics_sample_order(path: Path) -> list[str]:
    if path.suffix.lower() == ".xlsx":
        with zipfile.ZipFile(path) as zf:
            shared_strings = read_shared_strings(zf)
            with zf.open("xl/worksheets/sheet1.xml") as handle:
                for _, elem in ET.iterparse(handle, events=("end",)):
                    if elem.tag == f"{NS_MAIN}row":
                        header = row_values(elem, shared_strings)
                        return header[2:]
        raise ValueError(f"No rows found in {path}")

    with path.open(newline="") as handle:
        return next(csv.reader(handle))[2:]


def tcga_patient(sample: str) -> str:
    parts = sample.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else sample


def sample_key(sample: str) -> str:
    parts = sample.split("-")
    if len(parts) >= 4:
        return "-".join(parts[:4])
    return sample


def build_xena_column_map(header: list[str]) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for idx, sample in enumerate(header[1:]):
        mapping.setdefault(sample, idx + 1)
        mapping.setdefault(sample_key(sample), idx + 1)
        mapping.setdefault(tcga_patient(sample), idx + 1)
    return mapping


def summarize_values(values: list[float]) -> dict[str, float]:
    if not values:
        return {}
    values = sorted(values)
    def q(frac: float) -> float:
        return values[int((len(values) - 1) * frac)]
    return {
        "min": values[0],
        "q01": q(0.01),
        "q25": q(0.25),
        "median": q(0.5),
        "q75": q(0.75),
        "q99": q(0.99),
        "max": values[-1],
    }


def prepare_one(cancer: str, xena_gz: Path, omics_path: Path, out_csv: Path, summary_csv: Path) -> None:
    omics_samples = read_omics_sample_order(omics_path)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    raw_values: list[float] = []
    counts_values: list[float] = []
    matched_samples: list[str] = []
    selected_indices: list[int] = []

    with gzip.open(xena_gz, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        xena_map = build_xena_column_map(header)
        for sample in omics_samples:
            idx = xena_map.get(sample) or xena_map.get(sample_key(sample)) or xena_map.get(tcga_patient(sample))
            if idx is not None:
                matched_samples.append(sample)
                selected_indices.append(idx)

        with out_csv.open("w", newline="") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["Gene"] + matched_samples)
            genes = 0
            for row in reader:
                gene = row[0].split("|")[0]
                out_row: list[str] = [gene]
                genes += 1
                for idx in selected_indices:
                    if idx >= len(row) or row[idx] == "":
                        out_row.append("")
                        continue
                    raw = float(row[idx])
                    count = max(0.0, math.pow(2.0, raw) - 1.0)
                    raw_values.append(raw)
                    counts_values.append(count)
                    out_row.append(str(int(round(count))))
                writer.writerow(out_row)

    raw_summary = summarize_values(raw_values)
    count_summary = summarize_values(counts_values)
    with summary_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["field", "value"])
        writer.writerow(["cancer", cancer])
        writer.writerow(["xena_file", str(xena_gz)])
        writer.writerow(["omics_file", str(omics_path)])
        writer.writerow(["omics_samples", len(omics_samples)])
        writer.writerow(["matched_samples", len(matched_samples)])
        writer.writerow(["genes", genes])
        writer.writerow(["transform", "round(2^x - 1)"])
        for key, value in raw_summary.items():
            writer.writerow([f"log2_{key}", value])
        for key, value in count_summary.items():
            writer.writerow([f"reconstructed_{key}", value])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", required=True)
    parser.add_argument("--xena-gz", required=True, type=Path)
    parser.add_argument("--omics", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--summary-csv", required=True, type=Path)
    args = parser.parse_args()
    prepare_one(args.cancer, args.xena_gz, args.omics, args.out_csv, args.summary_csv)


if __name__ == "__main__":
    main()
