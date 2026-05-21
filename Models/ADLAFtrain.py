# encoding: utf-8
import os
import os.path
import random
import numpy as np
import torch

from config import LOG_DIR, TRAIN_DIR, PRETRAIN_DIR, clinic_info_path
from utils.sampling import KFold, resample
from torch_geometric.loader import DataLoader
from tqdm.auto import tqdm

from Models.GraphSurv import GraphSurv
from Models.CoxNN import DeepCox_LossFunc
from utils.Indicators import concordance_index
from utils.support import split_censor, sort_data
from utils.EarlyStopping import EarlyStopping
from utils.Tensorboard import SummaryWriter

import contextlib, io

torch.set_printoptions(profile="full")

writer = SummaryWriter(log_dir=str(LOG_DIR))

# =========================
# Silent helper
# =========================
@contextlib.contextmanager
def suppress_output():
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield

# =========================
# Utils
# =========================
def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

@torch.no_grad()
def get_augmented_features(data, wcvae_model, concat, device):
    wcvae_model.eval()
    out = []
    for _ in range(concat):
        z = torch.randn([data.x.size(0), wcvae_model.latent_size], device=device)
        aug = wcvae_model.inference(z, data.x).detach()
        out.append(aug)
    return torch.cat(out, dim=1)

def _to_time_tensor(time_np, device):
    t = torch.as_tensor(time_np, dtype=torch.float32, device=device)
    if t.ndim == 1:
        t = t.view(-1, 1)
    return t

def _is_finite(x: torch.Tensor) -> bool:
    return bool(torch.isfinite(x).all().item())

def _lr_tag(lr: float) -> str:
    return f"{lr:.0e}".replace("+", "")

# =========================
# Train / Evaluate
# =========================
def train(cancer_name, data, lr, value, gcn_layers: int = 2):
    """
    改动点：
    1) 使用 ReduceLROnPlateau（监控 1 - validate_Cindex）
    2) GraphSurv 支持 gcn_layers（2/3/4）
    其余流程不动
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    set_seed(42)

    # ---- load clinic info ----
    clinic_info = torch.load(str(clinic_info_path(cancer_name)), map_location="cpu", weights_only=False)
    event = np.asarray(clinic_info["event"])
    duration = np.asarray(clinic_info["duration"])

    # ---- split indices by event ----
    event_indices = np.where(event == 1)[0]
    non_event_indices = np.where(event == 0)[0]

    # ---- bootstrap sampling to form train+vali ----
    event_sample_size = round(len(event_indices) * 0.8)
    non_event_sample_size = round(len(non_event_indices) * 0.8)

    sampled_index = set()
    for _ in range(10):
        sampled_event = resample(event_indices, n_samples=event_sample_size, replace=False, random_state=42)
        sampled_non_event = resample(non_event_indices, n_samples=non_event_sample_size, replace=False, random_state=42)
        sampled = np.concatenate([sampled_event, sampled_non_event])
        sampled_index |= set(sampled)

    sampled_list = sorted(list(sampled_index))
    train_validate_data = [data[i] for i in sampled_list]
    train_validate_event = event[sampled_list]
    train_validate_duration = duration[sampled_list]

    test_index = sorted(list(set(np.arange(len(data))) - sampled_index))
    test_data = [data[i] for i in test_index]
    test_time = duration[test_index]

    # ---- censor split ----
    censored_time, censored_data, uncensored_time, uncensored_data = split_censor(
        data=train_validate_data,
        status=np.array(train_validate_event),
        time=np.array(train_validate_duration)
    )
    censored_data.extend(uncensored_data)
    censored_time = np.vstack((censored_time, uncensored_time))

    # ---- KFold ----
    n_splits = int(os.environ.get("ADLAF_N_SPLITS", "10"))
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    n_folds = kf.get_n_splits()
    max_epochs = int(os.environ.get("ADLAF_MAX_EPOCHS", "60"))

    # ---- run dir（按 layers + lr 分开存）----
    run_dir = TRAIN_DIR / cancer_name / f"L{int(gcn_layers)}" / f"lr_{_lr_tag(float(lr))}"
    os.makedirs(run_dir, exist_ok=True)

    best_vali = []

    # ===== 单一总进度条：每个 epoch 更新一次 =====
    total_steps = n_folds * max_epochs
    pbar = tqdm(total=total_steps, desc=f"{cancer_name}|L={int(gcn_layers)}|lr={float(lr):.0e}", dynamic_ncols=True)

    # ✅ KEGG only
    type1 = "KEGG"
    wcvae_model1 = torch.load(
        str(PRETRAIN_DIR / cancer_name / f"{cancer_name}_{type1}_pretrain.pth"),
        weights_only=False
    ).to(device)

    try:
        for fold_id, (train_index, validate_index) in enumerate(kf.split(X=censored_data), start=1):
            train_data = [censored_data[i] for i in train_index]
            train_time_np = censored_time[train_index]
            validate_data = [censored_data[i] for i in validate_index]
            validate_time_np = censored_time[validate_index]

            # sort train
            _, train_data, train_time_np = sort_data(train_data, train_time_np)

            # time tensor for loss
            train_time_t = _to_time_tensor(train_time_np, device)
            validate_time_t = _to_time_tensor(validate_time_np, device)

            train_loader = DataLoader(dataset=train_data, batch_size=1, shuffle=False)
            validate_loader = DataLoader(dataset=validate_data, batch_size=1, shuffle=False)
            test_loader = DataLoader(dataset=test_data, batch_size=1, shuffle=False)

            num_train = len(train_data)
            num_validate = len(validate_data)

            # ✅ 只改：传入 gcn_layers
            model = GraphSurv(input_dim=4, n_layer=int(gcn_layers)).to(device)
            loss_func = DeepCox_LossFunc()
            early_stopping = EarlyStopping(patience=30, verbose=False)

            optimizer = torch.optim.Adam(model.parameters(), lr=float(lr), weight_decay=1.2e-4)

            # ✅ NEW: ReduceLROnPlateau（监控 1 - val_cindex）
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", patience=5, factor=0.5)


            vali_list = []

            for epoch in range(max_epochs):
                if getattr(early_stopping, "early_stop", False):
                    skipped = max_epochs - epoch
                    pbar.total -= skipped
                    pbar.refresh()
                    break

                # -------- train --------
                model.train()
                risks = torch.zeros([num_train, 1], dtype=torch.float32, device=device)

                for bid, graphs in enumerate(train_loader):
                    graphs = graphs.to(device)

                    aug1 = get_augmented_features(graphs, wcvae_model1, concat=1, device=device)
                    graphs.x1 = torch.cat((graphs.x, aug1), dim=1)

                    _, h1 = model(graphs.x1, graphs.adj_t)

                    graph_embed = torch.mean(h1, dim=0, keepdim=True)
                    risk = model.CoxNN(graph_embed)
                    risks[bid] = risk

                if not _is_finite(risks):
                    raise FloatingPointError(f"[{cancer_name}] risks contains NaN/Inf")

                loss = loss_func(risks, train_time_t)
                if not _is_finite(loss):
                    raise FloatingPointError(f"[{cancer_name}] loss is NaN/Inf")

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

                train_risks_np = risks.detach().cpu().numpy()
                train_Cindex = concordance_index(train_time_np, -train_risks_np)

                # -------- validate --------
                model.eval()
                with torch.no_grad():
                    validate_risks = torch.zeros([num_validate, 1], dtype=torch.float32, device=device)
                    for bid, graphs in enumerate(validate_loader):
                        graphs = graphs.to(device)

                        aug1 = get_augmented_features(graphs, wcvae_model1, concat=1, device=device)
                        graphs.x1 = torch.cat((graphs.x, aug1), dim=1)

                        _, h1 = model(graphs.x1, graphs.adj_t)

                        graph_embed = torch.mean(h1, dim=0, keepdim=True)
                        risk = model.CoxNN(graph_embed)
                        validate_risks[bid] = risk

                    if not _is_finite(validate_risks):
                        raise FloatingPointError(f"[{cancer_name}] validate_risks contains NaN/Inf")

                    validate_loss = loss_func(validate_risks, validate_time_t)
                    if not _is_finite(validate_loss):
                        raise FloatingPointError(f"[{cancer_name}] validate_loss is NaN/Inf")

                    validate_risks_np = validate_risks.detach().cpu().numpy()
                    validate_Cindex = concordance_index(validate_time_np, -validate_risks_np)

                    # ✅ NEW: plateau step（跟 early_stopping 同指标）
                    plateau_metric = float(1.0 - validate_Cindex)
                    scheduler.step(plateau_metric)

                    writer.add_scalar(f"{cancer_name}/L{int(gcn_layers)}/lr_{_lr_tag(float(lr))}/Fold{fold_id}/val_loss",
                                      validate_loss.item(), epoch)
                    writer.add_scalar(f"{cancer_name}/L{int(gcn_layers)}/lr_{_lr_tag(float(lr))}/Fold{fold_id}/val_cindex",
                                      validate_Cindex, epoch)
                    writer.add_scalar(f"{cancer_name}/L{int(gcn_layers)}/lr_{_lr_tag(float(lr))}/Fold{fold_id}/lr",
                                      optimizer.param_groups[0]["lr"], epoch)

                    with suppress_output():
                        early_stopping(
                            val_loss=1 - validate_Cindex,
                            model=model,
                            path=str(run_dir / f"{cancer_name}_fold{fold_id}_best.pth")
                        )

                vali_list.append(validate_Cindex)

                cur_lr = optimizer.param_groups[0]["lr"]
                pbar.set_postfix({
                    "fold": f"{fold_id}/{n_folds}",
                    "ep": f"{epoch+1}/{max_epochs}",
                    "trL": f"{loss.item():.4f}",
                    "trC": f"{train_Cindex:.4f}",
                    "vaL": f"{validate_loss.item():.4f}",
                    "vaC": f"{validate_Cindex:.4f}",
                    "lr":  f"{cur_lr:.6g}",
                })
                pbar.update(1)

            fold_best = float(np.max(vali_list)) if len(vali_list) else float("-inf")
            best_vali.append(fold_best)

        best_fold = int(np.argmax(best_vali)) + 1
        best_ckpt = run_dir / f"{cancer_name}_fold{best_fold}_best.pth"
        model_best = torch.load(str(best_ckpt), weights_only=False).to(device)

        test_Cindex, _ = evaluate(
            cancer_name=cancer_name,
            model=model_best,
            wcvae_model1=wcvae_model1,
            dataloader=test_loader,
            test_time=np.array(test_time),
            device=device
        )

        pbar.set_postfix({"best_fold": best_fold, "testC": f"{test_Cindex:.4f}"})
        return {
            "cancer": cancer_name,
            "lr": float(lr),
            "gcn_layers": int(gcn_layers),
            "best_fold": best_fold,
            "test_cindex": float(test_Cindex),
            "status": "OK",
            "error": ""
        }

    finally:
        pbar.close()


def evaluate(cancer_name, model, wcvae_model1, dataloader, test_time, device):
    model.eval()
    with torch.no_grad():
        risks = torch.zeros([len(test_time), 1], dtype=torch.float32, device=device)

        for bid, graphs in enumerate(dataloader):
            graphs = graphs.to(device)

            aug1 = get_augmented_features(graphs, wcvae_model1, concat=1, device=device)
            graphs.x1 = torch.cat((graphs.x, aug1), dim=1)

            _, h1 = model(graphs.x1, graphs.adj_t)

            graph_embed = torch.mean(h1, dim=0, keepdim=True)
            risk = model.CoxNN(graph_embed)
            risks[bid] = risk

        cindex = concordance_index(test_time, -risks.detach().cpu().numpy())
    return float(cindex), risks.detach()
