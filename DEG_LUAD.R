#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

user_lib <- "/home/zhouyj/R/x86_64-conda-linux-gnu-library/4.4"
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- which(args == name)
  if (length(hit) == 0 || hit == length(args)) default else args[[hit + 1]]
}

cancer <- get_arg("--cancer", "LUAD")
project_root <- normalizePath(get_arg("--project-root", "/home/zhouyj/project2"), mustWork = TRUE)
k <- get_arg("--k", "15")
base_dir <- file.path(project_root, paste0("k_", k))
raw_xlsx <- get_arg("--raw-xlsx", file.path(base_dir, "data", cancer, "raw", paste0("omicsdata_", cancer, ".xlsx")))
risk_csv <- get_arg("--risk-csv", file.path(base_dir, "biomarker", cancer, paste0(cancer, "_risk_groups.csv")))
out_dir <- get_arg("--out-dir", file.path(base_dir, "biomarker", cancer, "DEG"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
expr_tsv <- get_arg("--expr-tsv", file.path(out_dir, paste0(cancer, "_geneExp.tsv")))
logfc_cutoff <- as.numeric(get_arg("--logfc-cutoff", "1"))
p_cutoff <- as.numeric(get_arg("--p-cutoff", "0.05"))
p_column <- get_arg("--p-column", "P.Value")

message("Cancer: ", cancer)
message("Raw omics: ", raw_xlsx)
message("Risk groups: ", risk_csv)
message("Output: ", out_dir)
message("DEG threshold: |effect| > ", logfc_cutoff, ", ", p_column, " < ", p_cutoff)

extract_script <- file.path(project_root, "PANA", "utils", "extract_gene_exp_xlsx.py")
if (!file.exists(expr_tsv) || file.info(expr_tsv)$mtime < file.info(raw_xlsx)$mtime) {
  message("Extracting geneExp matrix: ", expr_tsv)
  status <- system2("python", c(extract_script, "--xlsx", raw_xlsx, "--output", expr_tsv))
  if (status != 0) {
    stop("Failed to extract geneExp rows from xlsx.")
  }
}

omics <- data.table::fread(expr_tsv, data.table = FALSE, check.names = FALSE)
sample_ids <- colnames(omics)[-(1:2)]

expr <- as.matrix(omics[, sample_ids, drop = FALSE])
mode(expr) <- "numeric"
rownames(expr) <- make.unique(omics$GeneSymbol)

risk <- read.csv(risk_csv, stringsAsFactors = FALSE)
risk$sample_id <- sample_ids[risk$sample_index + 1]
risk <- risk[match(colnames(expr), risk$sample_id), ]
if (any(is.na(risk$risk_group))) {
  stop("Risk groups could not be aligned to expression columns.")
}
risk$risk_group <- factor(risk$risk_group, levels = c("low_risk", "high_risk"))
write.csv(risk, file.path(out_dir, paste0(cancer, "_risk_groups_with_sample_id.csv")), row.names = FALSE)

design <- model.matrix(~ 0 + risk$risk_group)
colnames(design) <- c("low_risk", "high_risk")
fit <- lmFit(expr, design)
contrast <- makeContrasts(low_vs_high = low_risk - high_risk, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg$Gene <- rownames(deg)
if (!p_column %in% colnames(deg)) {
  stop("Unknown p-column: ", p_column, ". Use P.Value or adj.P.Val.")
}
deg$plot_p <- deg[[p_column]]
deg$status <- ifelse(deg$logFC > logfc_cutoff & deg[[p_column]] < p_cutoff, "upregulated",
                     ifelse(deg$logFC < -logfc_cutoff & deg[[p_column]] < p_cutoff, "downregulated", "not_significant"))
deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]

up <- deg[deg$status == "upregulated", ]
down <- deg[deg$status == "downregulated", ]
summary_df <- data.frame(
  cancer = cancer,
  samples = ncol(expr),
  genes = nrow(expr),
  high_risk = sum(risk$risk_group == "high_risk"),
  low_risk = sum(risk$risk_group == "low_risk"),
  logfc_cutoff = logfc_cutoff,
  p_column = p_column,
  p_cutoff = p_cutoff,
  upregulated = nrow(up),
  downregulated = nrow(down)
)

write.csv(deg, file.path(out_dir, paste0(cancer, "_DEG_all.csv")), row.names = FALSE)
write.csv(up, file.path(out_dir, paste0(cancer, "_DEG_up.csv")), row.names = FALSE)
write.csv(down, file.path(out_dir, paste0(cancer, "_DEG_down.csv")), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, paste0(cancer, "_DEG_summary.csv")), row.names = FALSE)
writeLines(up$Gene, file.path(out_dir, paste0(cancer, "_up_genes.txt")))
writeLines(down$Gene, file.path(out_dir, paste0(cancer, "_down_genes.txt")))
print(summary_df)

volcano <- ggplot(deg, aes(x = logFC, y = -log10(plot_p))) +
  geom_point(aes(color = status), alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c(
    downregulated = "#f0a73a",
    upregulated = "#504099",
    not_significant = "gray70"
  )) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "steelblue") +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black") +
  labs(title = paste0(cancer, " Differential Gene Analysis"),
       x = "effect size (low risk vs high risk)",
       y = paste0("-log10(", p_column, ")")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(out_dir, paste0(cancer, "_volcano.png")), volcano, width = 8, height = 6, dpi = 300)

heat_genes <- head(deg$Gene[deg$status != "not_significant"], 50)
if (length(heat_genes) >= 2) {
  expr_heat <- expr[heat_genes, , drop = FALSE]
  sorted_samples <- risk$sample_id[order(risk$risk_group)]
  expr_heat <- expr_heat[, sorted_samples, drop = FALSE]
  expr_heat <- t(scale(t(expr_heat)))
  expr_heat[is.na(expr_heat)] <- 0
  expr_heat[expr_heat > 2.5] <- 2.5
  expr_heat[expr_heat < -2.5] <- -2.5
  annotation <- data.frame(Risk = risk$risk_group[match(sorted_samples, risk$sample_id)])
  rownames(annotation) <- sorted_samples
  png(file.path(out_dir, paste0(cancer, "_top_DEG_heatmap.png")), width = 1800, height = 1400, res = 220)
  pheatmap(expr_heat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           annotation_col = annotation,
           color = colorRampPalette(c("#2e59a7", "#f7f7f7", "#e60012"))(100),
           breaks = seq(-2.5, 2.5, length.out = 101),
           border_color = NA)
  dev.off()
}

enrich_terms <- function(gene_ids, universe_ids, term_to_genes, term_names, entrez_to_symbol, min_overlap = 2) {
  gene_ids <- intersect(unique(gene_ids), universe_ids)
  universe_ids <- unique(universe_ids)
  rows <- lapply(names(term_to_genes), function(term_id) {
    term_genes <- intersect(unique(term_to_genes[[term_id]]), universe_ids)
    overlap <- intersect(gene_ids, term_genes)
    if (length(overlap) < min_overlap) {
      return(NULL)
    }
    pvalue <- phyper(length(overlap) - 1,
                     length(term_genes),
                     length(universe_ids) - length(term_genes),
                     length(gene_ids),
                     lower.tail = FALSE)
    data.frame(
      ID = term_id,
      Description = unname(term_names[term_id]),
      GeneRatio = paste0(length(overlap), "/", length(gene_ids)),
      BgRatio = paste0(length(term_genes), "/", length(universe_ids)),
      pvalue = pvalue,
      Count = length(overlap),
      geneID = paste(unique(entrez_to_symbol[overlap]), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame())
  }
  result$p.adjust <- p.adjust(result$pvalue, method = "BH")
  result <- result[order(result$p.adjust, result$pvalue), ]
  rownames(result) <- NULL
  result
}

save_enrich_dotplot <- function(result, path, title) {
  if (nrow(result) == 0) {
    return()
  }
  plot_df <- head(result, 20)
  plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
  dotplot <- ggplot(plot_df, aes(x = Count, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_continuous(low = "#d7191c", high = "#2c7bb6", trans = "reverse") +
    labs(title = title, x = "Gene count", y = NULL, color = "BH adjusted P") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(path, dotplot, width = 9, height = 7, dpi = 300)
}

if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("GO.db", quietly = TRUE) &&
    requireNamespace("KEGGREST", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(GO.db)
    library(KEGGREST)
  })

  universe_map <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = unique(omics$GeneSymbol),
                                        keytype = "SYMBOL",
                                        columns = c("SYMBOL", "ENTREZID"))
  universe_map <- universe_map[!is.na(universe_map$ENTREZID), ]
  universe_ids <- unique(universe_map$ENTREZID)
  entrez_to_symbol <- universe_map$SYMBOL
  names(entrez_to_symbol) <- universe_map$ENTREZID

  go_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = universe_ids,
                                  keytype = "ENTREZID",
                                  columns = c("GO", "ONTOLOGY"))
  go_map <- go_map[!is.na(go_map$GO) & go_map$ONTOLOGY == "BP", ]
  go_terms <- split(go_map$ENTREZID, go_map$GO)
  go_names <- vapply(names(go_terms), function(go_id) {
    term <- GO.db::GOTERM[[go_id]]
    if (is.null(term)) go_id else AnnotationDbi::Term(term)
  }, character(1))

  kegg_links <- KEGGREST::keggLink("pathway", "hsa")
  kegg_gene_ids <- sub("^hsa:", "", names(kegg_links))
  kegg_path_ids <- sub("^path:", "", as.character(kegg_links))
  keep <- kegg_gene_ids %in% universe_ids
  kegg_terms <- split(kegg_gene_ids[keep], kegg_path_ids[keep])
  kegg_names <- KEGGREST::keggList("pathway", "hsa")
  names(kegg_names) <- sub("^path:", "", names(kegg_names))
  kegg_names <- sub(" - Homo sapiens \\(human\\)$", "", kegg_names)

  gene_sets <- list(up = unique(up$Gene), down = unique(down$Gene))
  for (direction in names(gene_sets)) {
    genes_for_enrich <- gene_sets[[direction]]
    if (length(genes_for_enrich) < 5) {
      next
    }
    gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = genes_for_enrich,
                                      keytype = "SYMBOL",
                                      columns = c("SYMBOL", "ENTREZID"))
    gene_ids <- unique(gene_map$ENTREZID[!is.na(gene_map$ENTREZID)])
    go_result <- enrich_terms(gene_ids, universe_ids, go_terms, go_names, entrez_to_symbol)
    kegg_result <- enrich_terms(gene_ids, universe_ids, kegg_terms, kegg_names, entrez_to_symbol)
    write.csv(go_result, file.path(out_dir, paste0(cancer, "_GO_BP_", direction, ".csv")), row.names = FALSE)
    write.csv(kegg_result, file.path(out_dir, paste0(cancer, "_KEGG_", direction, ".csv")), row.names = FALSE)
    save_enrich_dotplot(go_result, file.path(out_dir, paste0(cancer, "_GO_BP_", direction, "_dotplot.png")),
                        paste0(cancer, " GO BP ", direction))
    save_enrich_dotplot(kegg_result, file.path(out_dir, paste0(cancer, "_KEGG_", direction, "_dotplot.png")),
                        paste0(cancer, " KEGG ", direction))
  }
} else {
  message("Skipping GO/KEGG enrichment: org.Hs.eg.db, AnnotationDbi, GO.db, or KEGGREST is not installed.")
  message("Use generated *_up_genes.txt and *_down_genes.txt after installing those packages.")
}
