# encoding: utf-8
import argparse
import os.path
import csv
from pathlib import Path

import torch
from config import BASE_PATH, K, LOG_DIR, TRAIN_DIR, PRETRAIN_DIR, cancer_data_dir, clinic_info_path
from utils.Indicators import concordance_index
from torch_geometric.loader import DataLoader
import numpy as np
from utils.DataProcessing import CancerDataset
from utils.Tensorboard import SummaryWriter
import torch_geometric.transforms as T

import random

import warnings

# 忽略FutureWarning和UserWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


torch.set_printoptions(profile='full')

writer = SummaryWriter(log_dir=str(LOG_DIR))

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

set_seed(40)

def get_augmented_features(data, cvae_model, concat, device):
    augmented_features_list = []
    for _ in range(concat):  
        z = torch.randn([data.x.size(0), cvae_model.latent_size]).to(device)  
        augmented_features = cvae_model.inference(z, data.x).detach()
        augmented_features_list.append(augmented_features)
    return torch.cat(augmented_features_list, dim=1)  

def default_checkpoint(cancer_name, layer, lr_tag, fold):
    return TRAIN_DIR / cancer_name / f"L{layer}" / lr_tag / f"{cancer_name}_fold{fold}_best.pth"


def run(cancer_name, data, checkpoint_path, output_path):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    clinic_info = torch.load(str(clinic_info_path(cancer_name)), map_location='cpu', weights_only=False)
    event, duration = clinic_info['event'], clinic_info['duration']

    loader = DataLoader(dataset=data, batch_size=1, shuffle=False)

    model = torch.load(str(checkpoint_path), map_location=device, weights_only=False).to(device)
    cvae_model = torch.load(str(PRETRAIN_DIR / cancer_name / f"{cancer_name}_KEGG_pretrain.pth"), weights_only=False).to(
        device=device)
    cvae_model.eval()

    Cindex, risks = evaluate(model=model, cvae_model=cvae_model, dataloader=loader,
                             test_time=duration, device=device)
    print(f"{cancer_name} C-index: {Cindex:.4f}")

    risks_cpu = risks.cpu()
    median = risks_cpu.median()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer_csv = csv.writer(f)
        writer_csv.writerow(["sample_index", "risk_score", "risk_group", "duration", "event"])
        for index, item in enumerate(risks_cpu.numpy()):
            group = "high_risk" if item[0] > median else "low_risk"
            writer_csv.writerow([index, float(item[0]), group, float(np.squeeze(duration[index])), int(np.squeeze(event[index]))])
    print(f"Saved risk groups: {output_path}")


def evaluate(model, cvae_model, dataloader, test_time, device):
    model.eval()
    with torch.no_grad():
        risks = torch.zeros([test_time.shape[0], 1], dtype=torch.float, device=device)
        for id, graphs in enumerate(dataloader):
            graphs = graphs.to(device)
            augmented_features = get_augmented_features(data=graphs, cvae_model=cvae_model, concat=1,
                                                        device=device)
            graphs.x = torch.cat((graphs.x, augmented_features), dim=1)  # 拼接特征

            risk, _ = model(graphs.x, graphs.adj_t)
            risk = risk.mean()
            risks[id] = risk

        risks_save = risks.detach()
        cindex = concordance_index(test_time, -risks_save.cpu().numpy())
    return cindex, risks_save


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", default="LUAD")
    parser.add_argument("--layer", type=int, default=2)
    parser.add_argument("--fold", type=int, default=1)
    parser.add_argument("--lr-tag", default="lr_1e-03")
    parser.add_argument("--checkpoint", default=None)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    cancer = args.cancer
    checkpoint = args.checkpoint
    if checkpoint is None:
        checkpoint = default_checkpoint(cancer, args.layer, args.lr_tag, args.fold)
    else:
        checkpoint = Path(os.path.abspath(checkpoint))

    output = args.output
    if output is None:
        output = BASE_PATH / f"k_{K}" / "biomarker" / cancer / f"{cancer}_risk_groups.csv"
    else:
        output = Path(os.path.abspath(output))

    data = CancerDataset(root=str(cancer_data_dir(cancer)), transform=T.ToSparseTensor())
    run(cancer_name=cancer, data=data, checkpoint_path=checkpoint, output_path=output)
