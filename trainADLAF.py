import csv
import warnings
import argparse

from config import BASE_PATH, DATA_DIR, DEFAULT_CANCERS, K
from Models.ADLAFtrain import train
from utils.DataProcessing import CancerDataset
import torch_geometric.transforms as T

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

def append_csv(path, row, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        if write_header:
            w.writeheader()
        w.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancers", nargs="+", default=DEFAULT_CANCERS)
    parser.add_argument("--layers", nargs="+", type=int, default=[2, 3, 4])
    parser.add_argument("--lr", type=float, default=1e-3)
    args = parser.parse_args()

    cancers = args.cancers
    LAYER_LIST = args.layers
    LR_INIT = args.lr

    LOG_CSV = BASE_PATH / f"k_{K}" / "gcnlayer_sweep_plateau.csv"
    HEADER = ["cancer", "gcn_layers", "lr_init", "status", "best_fold", "test_cindex", "error"]

    value = 0.7

    for cancer_name in cancers:
        data = CancerDataset(
            root=str(DATA_DIR / cancer_name),
            transform=T.ToSparseTensor()
        )

        for gcn_layers in LAYER_LIST:
            print(f"Dataset: {cancer_name} | gcn_layers: {gcn_layers} | lr_init: {LR_INIT}")

            try:
                res = train(
                    cancer_name=cancer_name,
                    data=data,
                    lr=LR_INIT,
                    value=value,
                    gcn_layers=gcn_layers
                )

                row = {
                    "cancer": res.get("cancer", cancer_name),
                    "gcn_layers": int(res.get("gcn_layers", gcn_layers)),
                    "lr_init": float(res.get("lr", LR_INIT)),
                    "status": res.get("status", "OK"),
                    "best_fold": res.get("best_fold", ""),
                    "test_cindex": res.get("test_cindex", ""),
                    "error": res.get("error", ""),
                }
                append_csv(LOG_CSV, row, HEADER)

            except Exception as e:
                row = {
                    "cancer": cancer_name,
                    "gcn_layers": int(gcn_layers),
                    "lr_init": float(LR_INIT),
                    "status": "ERROR",
                    "best_fold": "",
                    "test_cindex": "",
                    "error": f"{type(e).__name__}: {str(e)}",
                }
                append_csv(LOG_CSV, row, HEADER)
                continue


