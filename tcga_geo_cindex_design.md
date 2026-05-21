# C-index Benchmark Design

## Cohorts

- TCGA cohort: use the existing 10 local cancers under `k_15/data/*/processed`.
- GEO cohort: use the five multi-omics datasets listed in `PANA/geo_multiomics_manifest.csv`.
- Do not count GSE96058/SCAN-B as one of the GEO multi-omics datasets because it is RNA-only.

## Models

- CoxNN
- DeepSurv
- GraphSurv
- LAGProg
- DeepHit
- MTLR
- ADLAF, fixed at `k = 15`

## Metrics

- Only C-index is reported.
- For multi-omics GEO data, normalize each omics block independently by feature-wise z-score before concatenating expression and methylation.
- For TCGA, use the local processed matrices that were generated with `k = 15`.

## Current GEO Processing State

- GSE20713_BRCA is already prepared.
- GSE75537_75538_OTSCC matrix files are downloaded; CpG-to-gene annotation is pending.
- GSE34861_34937_BALL matrix files are downloaded; custom platform feature-to-gene annotation is pending.
- GSE211483_211485_GEP_NEN metadata and mRNA matrix are downloaded; methylation processed matrix download is pending.
- GSE37815_37816_BLCA matrix files are downloaded; tumor sample filtering and probe-to-gene annotation are pending.

## Backup GEO Candidate

- GSE121720_121721_GBM has RNA/WGBS metadata with OS/PFS, but WGBS methylation is chromosome-level and needs heavier gene-level aggregation, so it is kept as a backup instead of the first five runnable GEO datasets.
