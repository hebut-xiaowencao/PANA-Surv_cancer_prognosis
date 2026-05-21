#!/usr/bin/env python3
"""Prepare paired BRCA expression and methylation matrices from GSE20713.

The two omics are normalized separately with row-wise z-scores before they are
concatenated, so expression and methylation keep equal footing.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET


XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


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


def split_symbols(value: str) -> list[str]:
    symbols: list[str] = []
    for chunk in re.split(r"///|;|,", value or ""):
        symbol = chunk.strip()
        if symbol and symbol != "---":
            symbols.append(symbol)
    return symbols


def load_gpl570_probe_map(path: Path, genes: set[str]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = None
        for line in handle:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                reader = csv.DictReader(handle, delimiter="\t")
                break
        if not in_table or reader is None:
            return mapping
        for row in reader:
            probe_id = row.get("ID", "")
            if probe_id == "!platform_table_end":
                break
            matched = sorted(set(split_symbols(row.get("Gene symbol", ""))) & genes)
            if matched:
                mapping[probe_id] = matched
    return mapping


def load_gpl8490_probe_map(path: Path, genes: set[str]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("IlmnID,"):
                header = next(csv.reader([line]))
                break
        else:
            return mapping
        reader = csv.DictReader(handle, fieldnames=header)
        for row in reader:
            probe_id = row.get("IlmnID", "")
            matched = sorted(set(split_symbols(row.get("Symbol", ""))) & genes)
            if probe_id and matched:
                mapping[probe_id] = matched
    return mapping


def patient_from_title(title: str) -> str | None:
    match = re.search(r"patient\s+(P2?_\w+|P_\w+)", title)
    if not match:
        return None
    patient = match.group(1)
    if re.search(r"_N\d+", patient):
        return None
    return patient.replace("P2_", "P_")


def parse_samples(path: Path) -> tuple[list[str], list[str], list[str], dict[str, dict[str, str]]]:
    titles: list[str] = []
    geo: list[str] = []
    sources: list[str] = []
    characteristics: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                sources = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                values = [v.strip().strip('"') for v in next(csv.reader([line], delimiter="\t"))[1:]]
                for value in values:
                    if ":" in value:
                        key = value.split(":", 1)[0].strip().lower().replace(".", "_").replace(" ", "_")
                        characteristics.setdefault(key, [])
                        break
                if values and ":" in values[0]:
                    key = values[0].split(":", 1)[0].strip().lower().replace(".", "_").replace(" ", "_")
                    characteristics[key] = [v.split(":", 1)[1].strip() if ":" in v else "" for v in values]
            elif line.startswith("!series_matrix_table_begin"):
                break
    patients = [patient_from_title(title) for title in titles]
    keep = [
        i for i, patient in enumerate(patients)
        if patient and i < len(sources) and sources[i] == "Breast tumor"
    ]
    clinical: dict[str, dict[str, str]] = {}
    for i in keep:
        patient = patients[i]
        if not patient:
            continue
        clinical[patient] = {
            "patient_id": patient,
            "geo_accession": geo[i] if i < len(geo) else "",
            "title": titles[i] if i < len(titles) else "",
        }
        for key, values in characteristics.items():
            if i < len(values):
                clinical[patient][key] = values[i]
    return [patients[i] for i in keep if patients[i]], [titles[i] for i in keep], [geo[i] for i in keep], clinical


def read_gene_matrix(path: Path, probe_map: dict[str, list[str]], wanted_patients: list[str]) -> dict[str, list[float]]:
    patients, _, _, _ = parse_samples(path)
    patient_to_col = {patient: i for i, patient in enumerate(patients)}
    keep_cols = [patient_to_col[p] for p in wanted_patients]

    sums: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(wanted_patients))
    counts: dict[str, list[int]] = defaultdict(lambda: [0] * len(wanted_patients))
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                break
        if not in_table:
            return {}
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        for row in reader:
            if not row or row[0] == "!series_matrix_table_end":
                break
            probe_id = row[0].strip().strip('"')
            genes = probe_map.get(probe_id)
            if not genes:
                continue
            values = []
            for col in keep_cols:
                raw = row[col + 1].strip().strip('"')
                try:
                    values.append(float(raw))
                except ValueError:
                    values.append(float("nan"))
            for gene in genes:
                for i, value in enumerate(values):
                    if not math.isnan(value):
                        sums[gene][i] += value
                        counts[gene][i] += 1

    matrix: dict[str, list[float]] = {}
    for gene, total in sums.items():
        row = []
        for value, count in zip(total, counts[gene]):
            row.append(value / count if count else float("nan"))
        if sum(not math.isnan(v) for v in row) >= max(3, len(wanted_patients) // 2):
            matrix[gene] = row
    return matrix


def zscore_rows(matrix: dict[str, list[float]]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for gene, values in matrix.items():
        clean = [v for v in values if not math.isnan(v)]
        if not clean:
            continue
        mean = sum(clean) / len(clean)
        var = sum((v - mean) ** 2 for v in clean) / len(clean)
        sd = math.sqrt(var) or 1.0
        result[gene] = [0.0 if math.isnan(v) else (v - mean) / sd for v in values]
    return result


def write_matrix(path: Path, rows: list[tuple[str, list[float]]], samples: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Feature"] + samples)
        for feature, values in rows:
            writer.writerow([feature] + [f"{v:.8g}" for v in values])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--relationships", type=Path, default=Path("relationships.xlsx"))
    parser.add_argument("--expr-matrix", type=Path, required=True)
    parser.add_argument("--meth-matrix", type=Path, required=True)
    parser.add_argument("--expr-annot", type=Path, required=True)
    parser.add_argument("--meth-annot", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("k_15/external/BRCA_GSE20713"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    relationship_genes = load_relationship_genes(args.relationships)
    expr_patients, _, expr_geo, expr_clinical = parse_samples(args.expr_matrix)
    meth_patients, _, meth_geo, meth_clinical = parse_samples(args.meth_matrix)
    paired = sorted(set(expr_patients) & set(meth_patients), key=lambda x: int(re.sub(r"\D", "", x)))

    expr_map = load_gpl570_probe_map(args.expr_annot, relationship_genes)
    meth_map = load_gpl8490_probe_map(args.meth_annot, relationship_genes)
    expr_matrix = zscore_rows(read_gene_matrix(args.expr_matrix, expr_map, paired))
    meth_matrix = zscore_rows(read_gene_matrix(args.meth_matrix, meth_map, paired))
    common_genes = sorted(set(expr_matrix) & set(meth_matrix))

    expr_rows = [(f"{gene}__geneExp", expr_matrix[gene]) for gene in common_genes]
    meth_rows = [(f"{gene}__methylation", meth_matrix[gene]) for gene in common_genes]
    write_matrix(args.out_dir / "GSE20713_relationship_multiomics_zscore.csv", expr_rows + meth_rows, paired)
    write_matrix(args.out_dir / "GSE20713_relationship_geneExp_zscore.csv", expr_rows, paired)
    write_matrix(args.out_dir / "GSE20713_relationship_methylation_zscore.csv", meth_rows, paired)

    geo_lookup_expr = dict(zip(expr_patients, expr_geo))
    geo_lookup_meth = dict(zip(meth_patients, meth_geo))
    with (args.out_dir / "GSE20713_paired_samples.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "patient_id",
            "expression_geo_accession",
            "methylation_geo_accession",
            "t_os",
            "e_os",
            "t_rfs",
            "e_rfs",
            "er_status",
            "her2_status",
            "grade",
            "node",
            "subtypeihc",
            "subtypege",
        ])
        for patient in paired:
            clinical = expr_clinical.get(patient, meth_clinical.get(patient, {}))
            writer.writerow([
                patient,
                geo_lookup_expr.get(patient, ""),
                geo_lookup_meth.get(patient, ""),
                clinical.get("t_os", ""),
                clinical.get("e_os", ""),
                clinical.get("t_rfs", ""),
                clinical.get("e_rfs", ""),
                clinical.get("er_status", clinical.get("er", "")),
                clinical.get("her2_status", clinical.get("her2", "")),
                clinical.get("grade", ""),
                clinical.get("node", ""),
                clinical.get("subtypeihc", clinical.get("subtype_ihc", "")),
                clinical.get("subtypege", ""),
            ])

    with (args.out_dir / "GSE20713_multiomics_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"relationship_genes={len(relationship_genes)}\n")
        handle.write(f"paired_tumor_samples={len(paired)}\n")
        handle.write(f"expression_relationship_genes={len(expr_matrix)}\n")
        handle.write(f"methylation_relationship_genes={len(meth_matrix)}\n")
        handle.write(f"common_multiomics_genes={len(common_genes)}\n")
        handle.write(f"concatenated_features={len(expr_rows) + len(meth_rows)}\n")
        handle.write("normalization=geneExp and methylation were row-zscored separately before concatenation\n")

    print(f"paired_tumor_samples={len(paired)}")
    print(f"expression_relationship_genes={len(expr_matrix)}")
    print(f"methylation_relationship_genes={len(meth_matrix)}")
    print(f"common_multiomics_genes={len(common_genes)}")
    print(f"concatenated_features={len(expr_rows) + len(meth_rows)}")


if __name__ == "__main__":
    main()
