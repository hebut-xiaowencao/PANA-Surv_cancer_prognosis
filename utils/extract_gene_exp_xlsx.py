#!/usr/bin/env python3
"""Extract gene expression rows from the first sheet of an omics xlsx file."""

from __future__ import annotations

import argparse
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
        parts = []
        for text_node in si.iter(f"{NS_MAIN}t"):
            parts.append(text_node.text or "")
        strings.append("".join(parts))
    return strings


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        parts = [node.text or "" for node in cell.iter(f"{NS_MAIN}t")]
        return "".join(parts)

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


def extract_gene_exp(xlsx_path: Path, output_path: Path) -> None:
    with zipfile.ZipFile(xlsx_path) as zf:
        shared_strings = read_shared_strings(zf)
        with zf.open("xl/worksheets/sheet1.xml") as handle:
            context = ET.iterparse(handle, events=("end",))
            header: list[str] | None = None
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w", encoding="utf-8") as out:
                for _, elem in context:
                    if elem.tag != f"{NS_MAIN}row":
                        continue
                    values = row_values(elem, shared_strings, len(header) if header else None)
                    if header is None:
                        header = values
                        out.write("\t".join(header) + "\n")
                    elif len(values) >= 2 and values[1] == "geneExp":
                        if len(values) < len(header):
                            values.extend([""] * (len(header) - len(values)))
                        out.write("\t".join(values[: len(header)]) + "\n")
                    elem.clear()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    extract_gene_exp(args.xlsx, args.output)


if __name__ == "__main__":
    main()
