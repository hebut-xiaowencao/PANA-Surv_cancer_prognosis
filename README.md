# PANA

PANA is a Python/R workflow for pan-cancer multi-omics survival analysis. It builds graph-structured cancer samples from gene expression and DNA methylation data, pretrains WCVAE feature generators, trains ADLAF survival models, and exports patient risk groups for downstream biomarker analysis.

## Repository Layout

```text
PANA/
  config.py                         # shared paths and default cancer list
  preADLAF.py                       # WCVAE pretraining
  trainADLAF.py                     # ADLAF model training
  CaseStudyADLAF.py                 # inference and risk-group export
  DEG_XENA_counts.R                 # DESeq2 differential-expression analysis
  Models/                           # model definitions
  utils/                            # data processing and helper scripts
  requirements.txt                  # Python dependencies
```

By default, outputs are written under `k_15/` at the project root:

```text
k_15/
  data/<CANCER>/raw/                # input Excel omics files
  data/<CANCER>/processed/          # PyTorch Geometric processed graphs
  Pretrain/<CANCER>/                # WCVAE checkpoints
  Train/<CANCER>/                   # ADLAF checkpoints
  biomarker/<CANCER>/               # risk groups and DEG outputs
```

## Requirements

Install Python dependencies:

```bash
cd /home/zhouyj/project2
pip install -r PANA/requirements.txt
```

The main Python dependencies are:

- `numpy`
- `pandas`
- `scipy`
- `torch`
- `torch-geometric`
- `tensorboard`
- `matplotlib`
- `networkx`
- `tqdm`
- `openpyxl`

For downstream DEG analysis, install R packages:

```r
install.packages(c("BiocManager", "ggplot2", "pheatmap"))
BiocManager::install(c("DESeq2", "biomaRt", "clusterProfiler"))
```

The current tested environment in this workspace is:

- Python `3.13.2`
- R `4.4.3`
- Bioconductor `3.20`
- DESeq2 `1.46.0`
- biomaRt `2.62.1`

## Input Data

Each cancer type should have a raw omics Excel file under:

```text
k_15/data/<CANCER>/raw/
```

The loader reads the first file in that `raw/` directory. The file should contain:

- `Platform`: marks rows as `geneExp` or `methylation`
- `GeneSymbol`: gene symbol for each feature row
- event and survival duration rows at the top, as expected by `utils/DataProcessing.py`
- sample columns for gene expression and methylation values

PANA also needs a gene relationship file. By default it is expected at:

```text
data/relationships.xlsx
```

You can override this path with `ADLAF_RELATIONSHIPS_PATH`.

## Configuration

`PANA/config.py` controls the default paths. The most useful environment variables are:

```bash
export ADLAF_BASE_PATH=/home/zhouyj/project2
export ADLAF_K=15
export ADLAF_RELATIONSHIPS_PATH=/home/zhouyj/project2/data/relationships.xlsx
```

If these are not set, PANA uses the project root as `ADLAF_BASE_PATH` and `k_15/` as the run directory.

Default cancer types are:

```text
BLCA BRCA CESC COAD HNSC LGG LUAD MESO SARC SKCM
```

## Quick Start

Run commands from the project root.

### 1. Pretrain WCVAE Models

```bash
PYTHONPATH=PANA python PANA/preADLAF.py
```

This processes available cancer datasets and saves KEGG/KNN WCVAE checkpoints under:

```text
k_15/Pretrain/<CANCER>/
```

Useful options:

```bash
PYTHONPATH=PANA python PANA/preADLAF.py \
  --latent_size 10 \
  --pretrain_lr 0.001 \
  --total_iterations 10000 \
  --batch_size 64 \
  --seed 42
```

### 2. Train ADLAF Survival Models

Train all default cancers with 2, 3, and 4 GCN layers:

```bash
PYTHONPATH=PANA python PANA/trainADLAF.py
```

Train selected cancers:

```bash
PYTHONPATH=PANA python PANA/trainADLAF.py \
  --cancers LUAD CESC LGG MESO \
  --layers 2 3 4 \
  --lr 0.001
```

Training summaries are appended to:

```text
k_15/gcnlayer_sweep_plateau.csv
```

Model checkpoints are saved under:

```text
k_15/Train/<CANCER>/L<LAYER>/lr_1e-03/
```

### 3. Export Risk Groups

After training, run inference for a cancer type and checkpoint:

```bash
PYTHONPATH=PANA python PANA/CaseStudyADLAF.py \
  --cancer LUAD \
  --layer 2 \
  --fold 1 \
  --lr-tag lr_1e-03
```

By default, this reads:

```text
k_15/Train/LUAD/L2/lr_1e-03/LUAD_fold1_best.pth
k_15/Pretrain/LUAD/LUAD_KEGG_pretrain.pth
```

and writes:

```text
k_15/biomarker/LUAD/LUAD_risk_groups.csv
```

You can also pass an explicit checkpoint and output path:

```bash
PYTHONPATH=PANA python PANA/CaseStudyADLAF.py \
  --cancer LUAD \
  --checkpoint /path/to/checkpoint.pth \
  --output /path/to/LUAD_risk_groups.csv
```

The risk-group CSV contains:

- `sample_index`
- `risk_score`
- `risk_group`
- `duration`
- `event`

## DEG Analysis

`DEG_XENA_counts.R` compares high-risk and low-risk groups using reconstructed Xena HiSeqV2 count matrices.

Example:

```bash
Rscript PANA/DEG_XENA_counts.R \
  --cancer LUAD \
  --project-root /home/zhouyj/project2 \
  --k 15
```

Default inputs:

```text
k_15/xena_hiseqv2/LUAD_xena_hiseqv2_reconstructed_counts.csv
k_15/biomarker/LUAD/LUAD_risk_groups.csv
```

Default output directory:

```text
k_15/biomarker/LUAD/DEG_XENA_counts/
```

Main outputs include:

- `<CANCER>_DEG_all.csv`
- `<CANCER>_DEG_up.csv`
- `<CANCER>_DEG_down.csv`
- `<CANCER>_DEG_summary.csv`
- `<CANCER>_volcano.png`
- `<CANCER>_top_DEG_heatmap.png`

## Utility Scripts

Some helper scripts in `PANA/utils/` are useful for preparing external data:

- `prepare_xena_hiseqv2_counts.py`: prepare UCSC Xena HiSeqV2 matrices for DEG analysis
- `prepare_gse20713_multiomics.py`: prepare paired GEO multi-omics data
- `prepare_gse96058_brca.py`: prepare BRCA GSE96058 data
- `inspect_geo_survival.py`: inspect GEO survival fields

Read the script headers and arguments before running them, because expected input files differ by dataset.

## Notes

- Run Python commands with `PYTHONPATH=PANA` so imports such as `from config import ...` resolve correctly.
- If raw data change, remove the corresponding `k_15/data/<CANCER>/processed/` files to force PyTorch Geometric to rebuild processed graphs.
- GPU is used automatically when CUDA is available; otherwise PANA runs on CPU.
- The package uses public database data and computational models. Platform names such as Illumina Infinium HumanMethylation450 BeadChip and UNC Illumina HiSeq are data-source platforms, not software packages with independent version numbers.
