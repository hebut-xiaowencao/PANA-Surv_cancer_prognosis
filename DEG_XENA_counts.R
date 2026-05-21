#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- which(args == name)
  if (length(hit) == 0 || hit == length(args)) default else args[[hit + 1]]
}

cancer <- get_arg("--cancer", "LUAD")
project_root <- normalizePath(get_arg("--project-root", "/home/zhouyj/project2"), mustWork = TRUE)
k <- get_arg("--k", "15")
base_dir <- file.path(project_root, paste0("k_", k))
counts_csv <- get_arg("--counts-csv", file.path(base_dir, "xena_hiseqv2", paste0(cancer, "_xena_hiseqv2_reconstructed_counts.csv")))
risk_csv <- get_arg("--risk-csv", file.path(base_dir, "biomarker", cancer, paste0(cancer, "_risk_groups.csv")))
out_dir <- get_arg("--out-dir", file.path(base_dir, "biomarker", cancer, "DEG_XENA_counts"))
logfc_cutoff <- as.numeric(get_arg("--logfc-cutoff", "1"))
padj_cutoff <- as.numeric(get_arg("--padj-cutoff", "0.05"))
min_count <- as.numeric(get_arg("--min-count", "10"))
min_samples <- as.numeric(get_arg("--min-samples", "3"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Cancer: ", cancer)
message("Counts: ", counts_csv)
message("Risk groups: ", risk_csv)
message("Output: ", out_dir)

counts_df <- read.csv(counts_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Gene" %in% colnames(counts_df)) {
  stop("Counts CSV must contain a Gene column.")
}

count_mat <- as.matrix(counts_df[, setdiff(colnames(counts_df), "Gene"), drop = FALSE])
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
rownames(count_mat) <- counts_df$Gene
count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)

risk <- read.csv(risk_csv, stringsAsFactors = FALSE)
if (!all(c("sample_index", "risk_group") %in% colnames(risk))) {
  stop("Risk CSV must contain sample_index and risk_group columns.")
}

sample_ids <- colnames(count_mat)
if ("sample_id" %in% colnames(risk)) {
  risk$sample_id <- risk$sample_id
} else {
  max_index <- max(risk$sample_index, na.rm = TRUE) + 1
  if (max_index > length(sample_ids)) {
    warning("Risk sample_index exceeds reconstructed count columns; matching available sample positions only.")
  }
  risk$sample_id <- sample_ids[risk$sample_index + 1]
}
risk <- risk[risk$sample_id %in% colnames(count_mat), ]
risk <- risk[match(colnames(count_mat), risk$sample_id), ]
risk <- risk[!is.na(risk$risk_group), ]
count_mat <- count_mat[, risk$sample_id, drop = FALSE]

if (any(is.na(risk$risk_group)) || ncol(count_mat) != nrow(risk)) {
  stop("Risk groups could not be aligned to count columns.")
}

risk$risk_group <- factor(risk$risk_group, levels = c("high_risk", "low_risk"))
keep <- rowSums(count_mat >= min_count) >= min_samples
filtered_mat <- count_mat[keep, , drop = FALSE]
message("Genes before filter: ", nrow(count_mat), "; after filter: ", nrow(filtered_mat))

coldata <- data.frame(row.names = colnames(filtered_mat), risk_group = risk$risk_group)
dds <- DESeqDataSetFromMatrix(countData = filtered_mat, colData = coldata, design = ~ risk_group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("risk_group", "low_risk", "high_risk"))
deg <- as.data.frame(res)
deg$Gene <- rownames(deg)
deg <- deg[order(deg$padj, deg$pvalue), ]
deg$status <- ifelse(!is.na(deg$padj) & deg$padj < padj_cutoff & deg$log2FoldChange > logfc_cutoff, "upregulated",
                     ifelse(!is.na(deg$padj) & deg$padj < padj_cutoff & deg$log2FoldChange < -logfc_cutoff, "downregulated", "not_significant"))
deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]

up <- deg[deg$status == "upregulated", ]
down <- deg[deg$status == "downregulated", ]
summary_df <- data.frame(
  cancer = cancer,
  samples = ncol(filtered_mat),
  genes_before_filter = nrow(count_mat),
  genes_after_filter = nrow(filtered_mat),
  high_risk = sum(risk$risk_group == "high_risk"),
  low_risk = sum(risk$risk_group == "low_risk"),
  logfc_cutoff = logfc_cutoff,
  padj_cutoff = padj_cutoff,
  min_count = min_count,
  min_samples = min_samples,
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

plot_df <- deg
plot_df$plot_p <- ifelse(is.na(plot_df$padj) | plot_df$padj <= 0, min(plot_df$padj[plot_df$padj > 0], na.rm = TRUE), plot_df$padj)
volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(plot_p))) +
  geom_point(aes(color = status), alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c(
    downregulated = "#f0a73a",
    upregulated = "#504099",
    not_significant = "gray70"
  )) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "steelblue") +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black") +
  labs(title = paste0(cancer, " DESeq2 DEG from Xena reconstructed counts"),
       x = "log2 fold change (low risk vs high risk)",
       y = "-log10(BH adjusted P)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(out_dir, paste0(cancer, "_volcano.png")), volcano, width = 8, height = 6, dpi = 300)

up_heat <- head(up$Gene[up$Gene %in% rownames(count_mat)], 25)
down_heat <- head(down$Gene[down$Gene %in% rownames(count_mat)], 25)
heat_genes <- unique(c(up_heat, down_heat))
if (length(heat_genes) < 2) {
  heat_genes <- head(deg$Gene[deg$status != "not_significant" & deg$Gene %in% rownames(count_mat)], 50)
}
if (length(heat_genes) >= 2) {
  sorted_samples <- risk$sample_id[order(risk$risk_group)]
  expr_heat <- log2(count_mat[heat_genes, sorted_samples, drop = FALSE] + 1)
  expr_heat <- t(scale(t(expr_heat)))
  expr_heat[is.na(expr_heat)] <- 0
  expr_heat[expr_heat > 2.5] <- 2.5
  expr_heat[expr_heat < -2.5] <- -2.5
  annotation <- data.frame(
    Risk = factor(ifelse(risk$risk_group[match(sorted_samples, risk$sample_id)] == "high_risk", "High", "Low"),
                  levels = c("High", "Low"))
  )
  rownames(annotation) <- sorted_samples
  annotation_colors <- list(Risk = c(High = "#f46d00", Low = "#00b050"))
  png(file.path(out_dir, paste0(cancer, "_top_DEG_heatmap.png")), width = 1800, height = 1400, res = 220)
  pheatmap(expr_heat,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           annotation_col = annotation,
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("#2f61b5", "#efc0cf", "#d7191c"))(120),
           breaks = seq(-2.5, 2.5, length.out = 121),
           border_color = NA,
           main = "row z-score of log2 expression")
  dev.off()
}

enrich_terms <- function(gene_ids, universe_ids, term_to_genes, term_names, entrez_to_symbol, min_overlap = 2) {
  gene_ids <- intersect(unique(gene_ids), universe_ids)
  universe_ids <- unique(universe_ids)
  if (length(gene_ids) == 0) {
    return(data.frame())
  }
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
      GeneRatioNum = length(overlap) / length(gene_ids),
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
  plot_df <- plot_df[order(plot_df$GeneRatioNum, plot_df$p.adjust), ]
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)
  dotplot <- ggplot(plot_df, aes(x = GeneRatioNum, y = Description, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    scale_color_continuous(low = "#d7191c", high = "#2c7bb6") +
    scale_size_continuous(range = c(2.2, 8)) +
    labs(title = title, x = "GeneRatio", y = NULL, color = "p.adjust") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.35),
      panel.grid.minor = element_blank()
    )
  ggsave(path, dotplot, width = 9, height = 7, dpi = 300)
}

if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("GO.db", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(GO.db)
  })

  universe_symbols <- unique(rownames(filtered_mat))
  universe_map <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = universe_symbols,
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

  kegg_terms <- list()
  kegg_names <- character(0)
  if (requireNamespace("KEGGREST", quietly = TRUE)) {
    suppressPackageStartupMessages(library(KEGGREST))
    kegg_ok <- tryCatch({
      kegg_links <- KEGGREST::keggLink("pathway", "hsa")
      kegg_gene_ids <- sub("^hsa:", "", names(kegg_links))
      kegg_path_ids <- sub("^path:", "", as.character(kegg_links))
      keep_kegg <- kegg_gene_ids %in% universe_ids
      kegg_terms <<- split(kegg_gene_ids[keep_kegg], kegg_path_ids[keep_kegg])
      kegg_names_raw <- KEGGREST::keggList("pathway", "hsa")
      names(kegg_names_raw) <- sub("^path:", "", names(kegg_names_raw))
      kegg_names <<- sub(" - Homo sapiens \\(human\\)$", "", kegg_names_raw)
      TRUE
    }, error = function(e) {
      message("Skipping KEGG enrichment: ", conditionMessage(e))
      FALSE
    })
  } else {
    kegg_ok <- FALSE
    message("Skipping KEGG enrichment: KEGGREST is not installed.")
  }

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
    write.csv(go_result, file.path(out_dir, paste0(cancer, "_GO_BP_", direction, ".csv")), row.names = FALSE)
    save_enrich_dotplot(go_result, file.path(out_dir, paste0(cancer, "_GO_BP_", direction, "_dotplot.png")),
                        paste0("Top 20 GO Terms (Biological Process)"))
    if (kegg_ok) {
      kegg_result <- enrich_terms(gene_ids, universe_ids, kegg_terms, kegg_names, entrez_to_symbol)
      write.csv(kegg_result, file.path(out_dir, paste0(cancer, "_KEGG_", direction, ".csv")), row.names = FALSE)
      save_enrich_dotplot(kegg_result, file.path(out_dir, paste0(cancer, "_KEGG_", direction, "_dotplot.png")),
                          paste0("Top 20 KEGG Pathways"))
    }
  }
} else {
  message("Skipping GO enrichment: org.Hs.eg.db, AnnotationDbi, or GO.db is not installed.")
}
