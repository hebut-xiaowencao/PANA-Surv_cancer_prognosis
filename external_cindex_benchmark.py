#!/usr/bin/env python3
"""Run C-index-only survival benchmarks on external BRCA cohorts."""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import torch


XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def c_index(duration: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    pairs = conc = tied = 0
    n = len(duration)
    for i in range(n):
        for j in range(i + 1, n):
            if event[i] == 1 and duration[i] < duration[j]:
                pairs += 1
                conc += risk[i] > risk[j]
                tied += risk[i] == risk[j]
            elif event[j] == 1 and duration[j] < duration[i]:
                pairs += 1
                conc += risk[j] > risk[i]
                tied += risk[j] == risk[i]
    return float((conc + 0.5 * tied) / pairs) if pairs else float("nan")


def col_to_index(cell_ref: str) -> int:
    letters = re.sub(r"\d", "", cell_ref)
    idx = 0
    for char in letters:
        idx = idx * 26 + ord(char.upper()) - ord("A") + 1
    return idx - 1


def load_edges(path: Path) -> list[tuple[str, str]]:
    with zipfile.ZipFile(path) as zf:
        shared: list[str] = []
        try:
            root = ET.parse(zf.open("xl/sharedStrings.xml")).getroot()
            for si in root.iter(f"{XLSX_NS}si"):
                shared.append("".join((t.text or "") for t in si.iter(f"{XLSX_NS}t")))
        except KeyError:
            pass
        edges: list[tuple[str, str]] = []
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
            if len(values) >= 2 and values[0] != "gene_x":
                edges.append((values[0].strip(), values[1].strip()))
            row.clear()
    return edges


def parse_time(value: str) -> float:
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return float("nan")
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return float("nan")
    number = float(match.group(0))
    if "y" in text.lower():
        return number * 365.25
    return number


def read_matrix(path: Path, id_col: str) -> tuple[list[str], list[str], np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = header[1:]
        features = []
        rows = []
        for row in reader:
            if not row:
                continue
            features.append(row[0])
            rows.append([float(x) if x not in ("", "NA") else 0.0 for x in row[1:]])
    matrix = np.asarray(rows, dtype=np.float32).T
    return samples, features, matrix


def row_zscore_matrix(x: np.ndarray) -> np.ndarray:
    mean = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, keepdims=True)
    sd[sd == 0] = 1
    return (x - mean) / sd


def load_gse96058(root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    samples, features, x = read_matrix(root / "GSE96058_relationship_gene_expression.csv", "Gene")
    x = row_zscore_matrix(x)
    clinical = {}
    with (root / "GSE96058_clinical.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            clinical[row["sample_id"]] = row
    keep = []
    duration = []
    event = []
    for i, sample in enumerate(samples):
        row = clinical.get(sample)
        if not row:
            continue
        t = parse_time(row.get("overall_survival_days", ""))
        e = row.get("overall_survival_event", "")
        if math.isnan(t) or t <= 0 or e in ("", "NA"):
            continue
        keep.append(i)
        duration.append(t)
        event.append(float(e))
    return x[keep], np.asarray(duration, dtype=np.float32), np.asarray(event, dtype=np.float32), features


def load_gse20713(root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    samples, features, x = read_matrix(root / "GSE20713_relationship_multiomics_zscore.csv", "Feature")
    clinical = {}
    with (root / "GSE20713_paired_samples.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            clinical[row["patient_id"]] = row
    keep = []
    duration = []
    event = []
    for i, sample in enumerate(samples):
        row = clinical.get(sample)
        if not row:
            continue
        t = parse_time(row.get("t_os", ""))
        e = row.get("e_os", "")
        if math.isnan(t) or t <= 0 or e in ("", "NA"):
            continue
        keep.append(i)
        duration.append(t)
        event.append(float(e))
    return x[keep], np.asarray(duration, dtype=np.float32), np.asarray(event, dtype=np.float32), features


def feature_gene(feature: str) -> str:
    return feature.split("__", 1)[0]


def graph_smooth(x: np.ndarray, features: list[str], edges: list[tuple[str, str]]) -> np.ndarray:
    feature_index = {f: i for i, f in enumerate(features)}
    gene_to_features: dict[str, list[int]] = {}
    for idx, feature in enumerate(features):
        gene_to_features.setdefault(feature_gene(feature), []).append(idx)
    adj = np.eye(len(features), dtype=np.float32)
    for a, b in edges:
        for ia in gene_to_features.get(a, []):
            for ib in gene_to_features.get(b, []):
                same_omic = ("__" not in features[ia] or "__" not in features[ib] or
                             features[ia].split("__", 1)[1] == features[ib].split("__", 1)[1])
                if same_omic:
                    adj[ia, ib] = 1.0
                    adj[ib, ia] = 1.0
    degree = adj.sum(axis=1, keepdims=True)
    degree[degree == 0] = 1
    return x @ (adj / degree).T


def graph_smooth_topk(x: np.ndarray, features: list[str], edges: list[tuple[str, str]], k: int) -> np.ndarray:
    feature_index = {f: i for i, f in enumerate(features)}
    gene_to_features: dict[str, list[int]] = {}
    for idx, feature in enumerate(features):
        gene_to_features.setdefault(feature_gene(feature), []).append(idx)

    centered = x - x.mean(axis=0, keepdims=True)
    norm = np.linalg.norm(centered, axis=0)
    norm[norm == 0] = 1
    normalized = centered / norm

    neighbors: dict[int, set[int]] = {i: {i} for i in range(len(features))}
    for a, b in edges:
        for ia in gene_to_features.get(a, []):
            for ib in gene_to_features.get(b, []):
                same_omic = ("__" not in features[ia] or "__" not in features[ib] or
                             features[ia].split("__", 1)[1] == features[ib].split("__", 1)[1])
                if same_omic:
                    neighbors[ia].add(ib)
                    neighbors[ib].add(ia)

    adj = np.zeros((len(features), len(features)), dtype=np.float32)
    for i, idxs in neighbors.items():
        candidates = [j for j in idxs if j != i]
        if candidates:
            scores = np.abs(normalized[:, candidates].T @ normalized[:, i])
            order = np.argsort(scores)[-k:]
            selected = [candidates[j] for j in order]
        else:
            selected = []
        adj[i, i] = 1.0
        for j in selected:
            adj[i, j] = 1.0

    degree = adj.sum(axis=1, keepdims=True)
    degree[degree == 0] = 1
    return x @ (adj / degree).T


def make_folds(n: int, k: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return np.array_split(idx, k)


class MLP(torch.nn.Module):
    def __init__(self, input_dim: int, hidden: list[int], dropout: float = 0.1, out_dim: int = 1):
        super().__init__()
        layers = []
        last = input_dim
        for h in hidden:
            layers += [torch.nn.Linear(last, h), torch.nn.ReLU(), torch.nn.Dropout(dropout)]
            last = h
        layers.append(torch.nn.Linear(last, out_dim))
        self.net = torch.nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def cox_loss(risk: torch.Tensor, duration: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    order = torch.argsort(duration, descending=True)
    risk = risk.view(-1)[order]
    event = event.view(-1)[order]
    log_cumsum = torch.logcumsumexp(risk, dim=0)
    observed = event > 0
    if observed.sum() == 0:
        return risk.sum() * 0
    return -((risk - log_cumsum)[observed]).mean()


def discrete_targets(duration: np.ndarray, event: np.ndarray, bins: np.ndarray) -> np.ndarray:
    return np.searchsorted(bins, duration, side="right").clip(0, len(bins))


def deephit_loss(logits: torch.Tensor, target: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    prob = torch.softmax(logits, dim=1)
    event_mask = event > 0
    eps = 1e-7
    loss_event = -torch.log(prob[torch.arange(len(target)), target].clamp_min(eps))[event_mask]
    surv = torch.flip(torch.cumsum(torch.flip(prob, dims=[1]), dim=1), dims=[1])
    loss_cens = -torch.log(surv[torch.arange(len(target)), target].clamp_min(eps))[~event_mask]
    parts = []
    if loss_event.numel():
        parts.append(loss_event.mean())
    if loss_cens.numel():
        parts.append(loss_cens.mean())
    return sum(parts) / len(parts)


def mtlr_loss(logits: torch.Tensor, target: torch.Tensor, event: torch.Tensor) -> torch.Tensor:
    hazard = torch.sigmoid(logits)
    idx = torch.arange(logits.shape[1], device=logits.device).view(1, -1)
    t = target.view(-1, 1)
    e = event.view(-1, 1) > 0
    y = (idx == t).float()
    cens_mask = idx >= t
    event_loss = torch.nn.functional.binary_cross_entropy(hazard, y, reduction="none").sum(1)
    cens_loss = -torch.log((1 - hazard).clamp_min(1e-7))
    cens_loss = (cens_loss * cens_mask.float()).sum(1)
    return torch.where(e.view(-1), event_loss, cens_loss).mean()


def fit_model(
    method: str,
    x_train: np.ndarray,
    duration_train: np.ndarray,
    event_train: np.ndarray,
    x_test: np.ndarray,
    duration_test: np.ndarray,
    event_test: np.ndarray,
    seed: int,
) -> float:
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    input_dim = x_train.shape[1]
    discrete = method in {"DeepHit", "MTLR"}
    if method == "CoxNN":
        model = MLP(input_dim, [128, 32], 0.05, 1)
    elif method == "DeepSurv":
        model = MLP(input_dim, [256, 128, 32], 0.2, 1)
    elif method == "LAGProg":
        model = MLP(input_dim, [], 0.0, 1)
    elif method == "GraphSurv":
        model = MLP(input_dim, [128, 32], 0.1, 1)
    elif method.startswith("ADLAF_k"):
        model = MLP(input_dim, [256, 64, 16], 0.15, 1)
    elif method in {"DeepHit", "MTLR"}:
        model = MLP(input_dim, [128, 64], 0.1, 12)
        bins = np.quantile(duration_train[event_train > 0], np.linspace(0.05, 0.95, 11))
        bins = np.unique(bins)
        model.net[-1] = torch.nn.Linear(model.net[-1].in_features, len(bins) + 1)
    else:
        raise ValueError(method)

    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    xt = torch.tensor(x_train, dtype=torch.float32)
    xv = torch.tensor(x_test, dtype=torch.float32)
    dt = torch.tensor(duration_train, dtype=torch.float32)
    et = torch.tensor(event_train, dtype=torch.float32)
    if discrete:
        tt = torch.tensor(discrete_targets(duration_train, event_train, bins), dtype=torch.long)
    best = None
    best_loss = float("inf")
    bad = 0
    for _ in range(120):
        model.train()
        opt.zero_grad()
        out = model(xt)
        if method == "DeepHit":
            loss = deephit_loss(out, tt, et)
        elif method == "MTLR":
            loss = mtlr_loss(out, tt, et)
        else:
            loss = cox_loss(out, dt, et)
        loss.backward()
        opt.step()
        value = float(loss.detach())
        if value < best_loss - 1e-4:
            best_loss = value
            best = {k: v.detach().clone() for k, v in model.state_dict().items()}
            bad = 0
        else:
            bad += 1
        if bad >= 15:
            break
    if best is not None:
        model.load_state_dict(best)
    model.eval()
    with torch.no_grad():
        out = model(xv)
        if method in {"DeepHit", "MTLR"}:
            prob = torch.softmax(out, dim=1) if method == "DeepHit" else torch.sigmoid(out)
            risk = prob.sum(dim=1).cpu().numpy()
        else:
            risk = out.view(-1).cpu().numpy()
    return c_index(duration_test, event_test, risk)


def benchmark_dataset(
    name: str,
    x: np.ndarray,
    duration: np.ndarray,
    event: np.ndarray,
    features: list[str],
    edges,
    folds: int,
    seed: int,
    adlaf_ks: list[int],
):
    graph_x = graph_smooth(x, features, edges)
    augmented = {
        "CoxNN": x,
        "DeepSurv": x,
        "DeepHit": x,
        "MTLR": x,
        "LAGProg": graph_x,
        "GraphSurv": np.concatenate([x, graph_x], axis=1),
    }
    for k in adlaf_ks:
        graph_x_k = graph_smooth_topk(x, features, edges, k)
        augmented[f"ADLAF_k{k}"] = np.concatenate([x, graph_x_k, x - graph_x_k], axis=1)
    fold_ids = make_folds(len(duration), folds, seed)
    rows = []
    for method, mat in augmented.items():
        for fold, test_idx in enumerate(fold_ids, start=1):
            train_idx = np.setdiff1d(np.arange(len(duration)), test_idx)
            mu = mat[train_idx].mean(axis=0, keepdims=True)
            sd = mat[train_idx].std(axis=0, keepdims=True)
            sd[sd == 0] = 1
            x_train = (mat[train_idx] - mu) / sd
            x_test = (mat[test_idx] - mu) / sd
            cidx = fit_model(method, x_train, duration[train_idx], event[train_idx], x_test, duration[test_idx], event[test_idx], seed + fold)
            rows.append({"dataset": name, "method": method, "fold": fold, "c_index": cidx, "n": len(duration), "events": int(event.sum())})
            print(name, method, fold, f"{cidx:.4f}", flush=True)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/home/zhouyj/project2"))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--adlaf-ks", nargs="+", type=int, default=[15])
    parser.add_argument("--output", type=Path, default=Path("k_15/external/external_cindex_benchmark.csv"))
    args = parser.parse_args()

    edges = load_edges(args.root / "relationships.xlsx")
    datasets = {
        "GSE96058_SCANB_RNA": load_gse96058(args.root / "k_15/external/BRCA_GSE96058"),
        "GSE20713_BRCA_multiomics": load_gse20713(args.root / "k_15/external/BRCA_GSE20713"),
    }
    rows = []
    for name, (x, duration, event, features) in datasets.items():
        rows.extend(benchmark_dataset(name, x, duration, event, features, edges, args.folds, args.seed, args.adlaf_ks))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "method", "fold", "c_index", "n", "events"])
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output.with_name(args.output.stem + "_summary.csv")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "method", "mean_c_index", "std_c_index", "folds"])
        writer.writeheader()
        for dataset in sorted({r["dataset"] for r in rows}):
            for method in sorted({r["method"] for r in rows}):
                vals = [r["c_index"] for r in rows if r["dataset"] == dataset and r["method"] == method]
                if vals:
                    writer.writerow({
                        "dataset": dataset,
                        "method": method,
                        "mean_c_index": float(np.mean(vals)),
                        "std_c_index": float(np.std(vals)),
                        "folds": len(vals),
                    })


if __name__ == "__main__":
    main()
