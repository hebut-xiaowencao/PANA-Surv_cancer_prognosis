#!/usr/bin/env python3
"""Prepare GSE96058 BRCA RNA data for relationship-gene validation."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DEFAULT_EXPR_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE96058&file=GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz&format=file"
)


def col_to_index(cell_ref: str) -> int:
    letters = re.sub(r"\d", "", cell_ref)
    idx = 0
    for char in letters:
        idx = idx * 26 + ord(char.upper()) - ord("A") + 1
    return idx - 1


def load_relationship_genes(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as zf:
        shared: list[str] = []
        try:
            root = ET.parse(zf.open("xl/sharedStrings.xml")).getroot()
            for si in root.iter(f"{XLSX_NS}si"):
                shared.append("".join((t.text or "") for t in si.iter(f"{XLSX_NS}t")))
        except KeyError:
            pass

        genes: set[str] = set()
        for _, row in ET.iterparse(zf.open("xl/worksheets/sheet1.xml"), events=("end",)):
            if row.tag != f"{XLSX_NS}row":
                continue
            values: list[str] = []
            for cell in row.findall(f"{XLSX_NS}c"):
                idx = col_to_index(cell.attrib["r"])
                while len(values) <= idx:
                    values.append("")
                value_node = cell.find(f"{XLSX_NS}v")
                value = "" if value_node is None or value_node.text is None else value_node.text
                if cell.attrib.get("t") == "s" and value:
                    value = shared[int(value)]
                values[idx] = value
            if values and values[0] != "gene_x":
                genes.update(v.strip() for v in values[:2] if v.strip())
            row.clear()
    return genes


def parse_series_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    samples: list[str] = []
    clinical: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                samples = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
                clinical = {sample: {"sample_id": sample} for sample in samples}
            elif line.startswith("!Sample_geo_accession"):
                values = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
                for sample, value in zip(samples, values):
                    clinical[sample]["geo_accession"] = value
            elif line.startswith("!Sample_characteristics_ch1"):
                values = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
                for sample, value in zip(samples, values):
                    if ":" not in value:
                        continue
                    key, val = value.split(":", 1)
                    key = key.strip().lower().replace(" ", "_").replace("-", "_")
                    clinical[sample][key] = val.strip()
            elif line.startswith("!series_matrix_table_begin"):
                break
    return samples, clinical


def write_clinical(clinical: dict[str, dict[str, str]], output: Path) -> None:
    fields = sorted({key for row in clinical.values() for key in row})
    preferred = [
        "sample_id",
        "geo_accession",
        "scan_b_external_id",
        "overall_survival_days",
        "overall_survival_event",
        "pam50_subtype",
        "age_at_diagnosis",
        "tumor_size",
        "lymph_node_status",
        "er_status",
        "pgr_status",
        "her2_status",
    ]
    fields = [field for field in preferred if field in fields] + [
        field for field in fields if field not in preferred
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(clinical[sample] for sample in clinical)


def filter_expression(expr_url: str, genes: set[str], samples: list[str], output: Path) -> int:
    sample_set = set(samples)
    req = urllib.request.Request(expr_url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as response:
        with gzip.GzipFile(fileobj=response) as gz:
            text = (line.decode("utf-8", errors="replace") for line in gz)
            reader = csv.reader(text)
            header = next(reader)
            keep = [idx for idx, col in enumerate(header) if idx > 0 and col in sample_set]
            output.parent.mkdir(parents=True, exist_ok=True)
            count = 0
            with output.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["Gene"] + [header[idx] for idx in keep])
                for row in reader:
                    if not row:
                        continue
                    gene = row[0].strip()
                    if gene in genes:
                        writer.writerow([gene] + [row[idx] for idx in keep])
                        count += 1
            return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relationships", type=Path, default=Path("relationships.xlsx"))
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--expr-url", default=DEFAULT_EXPR_URL)
    parser.add_argument("--out-dir", type=Path, default=Path("k_15/external/BRCA_GSE96058"))
    args = parser.parse_args()

    genes = load_relationship_genes(args.relationships)
    samples, clinical = parse_series_matrix(args.series_matrix)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_clinical(clinical, args.out_dir / "GSE96058_clinical.csv")
    matched = filter_expression(args.expr_url, genes, samples, args.out_dir / "GSE96058_relationship_gene_expression.csv")

    with (args.out_dir / "GSE96058_prepare_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"relationship_genes={len(genes)}\n")
        handle.write(f"clinical_samples={len(samples)}\n")
        handle.write(f"matched_expression_genes={matched}\n")

    print(f"relationship_genes={len(genes)}")
    print(f"clinical_samples={len(samples)}")
    print(f"matched_expression_genes={matched}")


if __name__ == "__main__":
    main()
