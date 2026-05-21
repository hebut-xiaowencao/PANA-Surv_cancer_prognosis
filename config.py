import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = Path(os.environ.get("ADLAF_BASE_PATH", PROJECT_ROOT)).resolve()
K = int(os.environ.get("ADLAF_K", "15"))

RUN_ROOT = BASE_PATH / f"k_{K}"
DATA_DIR = RUN_ROOT / "data"
TRAIN_DIR = RUN_ROOT / "Train"
PRETRAIN_DIR = RUN_ROOT / "Pretrain"
INDICATOR_DIR = RUN_ROOT / "indicator"
LOG_DIR = RUN_ROOT / "log"
RELATIONSHIPS_PATH = Path(
    os.environ.get("ADLAF_RELATIONSHIPS_PATH", BASE_PATH / "data" / "relationships.xlsx")
).resolve()


DEFAULT_CANCERS = ["BLCA", "BRCA", "CESC", "COAD", "HNSC", "LGG", "LUAD", "MESO", "SARC", "SKCM"]


def cancer_data_dir(cancer_name: str) -> Path:
    return DATA_DIR / cancer_name


def cancer_train_dir(cancer_name: str) -> Path:
    return TRAIN_DIR / cancer_name


def cancer_pretrain_dir(cancer_name: str) -> Path:
    return PRETRAIN_DIR / cancer_name


def clinic_info_path(cancer_name: str) -> Path:
    return cancer_data_dir(cancer_name) / "processed" / f"{cancer_name}_info.pkl"
