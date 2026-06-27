# =============================================================================
# 16_synovium_DE_analysis.R
# -----------------------------------------------------------------------------
# Stage 2: limma DE analysis on GSE55235 synovium
#   - Filter to OA and Normal samples only (discard RA)
#   - Map probes to gene symbols via hgu133a.db
#   - Apply thresholds matching Block B: padj < 0.05, |log2FC| >= 0.585
#
# INPUT  : data/raw/GSE55235_eset.rds
# OUTPUT : data/processed/syn_de_results.rds
#          data/processed/syn_degs_filtered.rds
#          tables/table_F1_synovium_DEG_summary.csv
# =============================================================================

# Install needed packages
required <- c("limma", "hgu133a.db")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    cat("Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(limma)
  library(hgu133a.db)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("STAGE 2: SYNOVIUM DE ANALYSIS\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load eset and filter to OA and Normal only
# -----------------------------------------------------------------------------
eset <- readRDS(file.path(PROJECT_ROOT, "data/raw/GSE55235_eset.rds"))
pheno <- pData(eset)

# Re-infer condition (in case Condition_inferred wasn't saved back to eset)
pheno$Condition <- dplyr::case_when(
  grepl("healthy", tolower(pheno$title))   ~ "Normal",
  grepl("osteoarth", tolower(pheno$title)) ~ "OA",
  grepl("rheumatoid", tolower(pheno$title)) ~ "RA",
  TRUE ~ "UNCLEAR"
)

keep <- pheno$Condition %in% c("Normal", "OA")
eset_sub <- eset[, keep]
pheno_sub <- pheno[keep, ]

cat("=== FILTERED TO OA AND NORMAL ===\n")
cat("Samples retained:", ncol(eset_sub), "\n")
print(table(pheno_sub$Condition))
cat("\n")

# -----------------------------------------------------------------------------
# 2. Check normalization status
# -----------------------------------------------------------------------------
expr <- exprs(eset_sub)
cat("=== EXPRESSION MATRIX QUALITY ===\n")
cat("Dims:", nrow(expr), "x", ncol(expr), "\n")
cat("Range:", round(range(expr, na.rm = TRUE), 2), "\n")
cat("Median:", round(median(expr, na.rm = TRUE), 2), "\n")

# If values are positive (e.g. 2-15 range), data is already log2-transformed
# If values are huge (e.g. 0-60000), need to log-transform
if (max(expr, na.rm = TRUE) > 100) {
  cat("Expression appears linear-scale; applying log2 transformation...\n")
  expr <- log2(expr + 1)
} else {
  cat("Expression appears already log2-transformed (GEOquery default).\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 3. limma DE analysis: OA vs Normal
# -----------------------------------------------------------------------------
condition <- factor(pheno_sub$Condition, levels = c("Normal", "OA"))
design <- model.matrix(~ condition)
colnames(design) <- c("Intercept", "OA_vs_Normal")

cat("=== DESIGN MATRIX ===\n")
print(table(condition))
cat("\n")

cat("Running limma...\n")
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Extract DE results for OA_vs_Normal coefficient
de_probes <- topTable(fit, coef = "OA_vs_Normal",
                      number = Inf, sort.by = "none")

cat("DE probes:", nrow(de_probes), "\n\n")

# -----------------------------------------------------------------------------
# 4. Map probes to gene symbols
# -----------------------------------------------------------------------------
cat("=== MAPPING PROBES TO GENE SYMBOLS ===\n")

probe_ids <- rownames(de_probes)
symbol_map <- AnnotationDbi::select(hgu133a.db,
                                    keys = probe_ids,
                                    columns = c("SYMBOL", "ENTREZID"),
                                    keytype = "PROBEID")

# Take first symbol per probe (handles multi-mappers simply)
symbol_map <- symbol_map %>%
  dplyr::group_by(PROBEID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

de_probes$Probe <- rownames(de_probes)
de_annot <- de_probes %>%
  dplyr::left_join(symbol_map, by = c("Probe" = "PROBEID")) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::rename(Gene = SYMBOL)

cat("Probes with gene symbol:", nrow(de_annot),
    "/", nrow(de_probes),
    sprintf("(%.1f%%)\n", 100 * nrow(de_annot) / nrow(de_probes)))

# Collapse to gene-level: for multiple probes per gene, keep the one with
# smallest adjusted p-value (standard practice)
de_gene <- de_annot %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Gene, logFC, AveExpr, P.Value, adj.P.Val, t, ENTREZID, Probe)

cat("Genes after collapsing probes:", nrow(de_gene), "\n\n")

# Rename columns to match Block B convention for consistency
de_gene <- de_gene %>%
  dplyr::rename(log2FoldChange = logFC, padj = adj.P.Val, pvalue = P.Value) %>%
  dplyr::mutate(Direction = ifelse(log2FoldChange > 0, "Up_in_OA", "Down_in_OA"))

saveRDS(de_gene, file.path(PROJECT_ROOT, "data/processed/syn_de_results.rds"))
cat("Saved: data/processed/syn_de_results.rds\n\n")

# -----------------------------------------------------------------------------
# 5. Apply Block-B matching thresholds
# -----------------------------------------------------------------------------
LOG2FC_THRESHOLD <- 0.585
PADJ_THRESHOLD   <- 0.05

syn_degs <- de_gene %>%
  dplyr::filter(padj < PADJ_THRESHOLD,
                abs(log2FoldChange) >= LOG2FC_THRESHOLD)

syn_up   <- syn_degs %>% dplyr::filter(log2FoldChange >=  LOG2FC_THRESHOLD)
syn_down <- syn_degs %>% dplyr::filter(log2FoldChange <= -LOG2FC_THRESHOLD)

cat("=== SYNOVIUM OA DEG SUMMARY ===\n")
cat(sprintf("Threshold: padj < %.2f AND |log2FC| >= %.3f\n",
            PADJ_THRESHOLD, LOG2FC_THRESHOLD))
cat(sprintf("Synovium OA DEGs total: %d\n",   nrow(syn_degs)))
cat(sprintf("  Up in OA:    %d\n",   nrow(syn_up)))
cat(sprintf("  Down in OA:  %d\n",   nrow(syn_down)))
cat(sprintf("Up/Down ratio: %.2f\n\n", nrow(syn_up) / max(1, nrow(syn_down))))

# -----------------------------------------------------------------------------
# 6. Anchor gene check
# -----------------------------------------------------------------------------
anchor_check <- c("ALDH2", "SOD2", "GPX3", "GPX4", "NQO1", "HMOX1", "NFE2L2",
                  "PRDX1", "ACSL4", "SLC3A2", "TFRC", "AGO2", "SRSF3",
                  "DDX17", "DGCR8", "DROSHA", "DICER1")

anchor_in_syn <- syn_degs %>%
  dplyr::filter(Gene %in% anchor_check) %>%
  dplyr::select(Gene, log2FoldChange, padj, Direction) %>%
  dplyr::arrange(padj)

cat("=== ANCHOR GENES IN SYNOVIUM DEGs ===\n")
if (nrow(anchor_in_syn) > 0) {
  print(anchor_in_syn, row.names = FALSE, digits = 3)
} else {
  cat("No anchor genes pass synovium DEG threshold.\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# 7. Save refined synovium DEG object (matching Block B structure)
# -----------------------------------------------------------------------------
syn_refined <- list(
  oa_degs          = syn_degs,
  oa_up            = syn_up,
  oa_down          = syn_down,
  anchor_in_oa     = anchor_in_syn,
  threshold_padj   = PADJ_THRESHOLD,
  threshold_log2FC = LOG2FC_THRESHOLD,
  source           = "GSE55235 synovium (OA vs Normal, limma)"
)

saveRDS(syn_refined,
        file.path(PROJECT_ROOT, "data/processed/syn_degs_filtered.rds"))

# DEG summary table
summary_tbl <- data.frame(
  Tissue   = "Synovium (GSE55235)",
  N_genes_tested  = nrow(de_gene),
  N_DEGs          = nrow(syn_degs),
  N_up_in_OA      = nrow(syn_up),
  N_down_in_OA    = nrow(syn_down),
  N_anchor_DE     = nrow(anchor_in_syn)
)
print(summary_tbl, row.names = FALSE)
write.csv(summary_tbl,
          file.path(PROJECT_ROOT, "tables/table_F1_synovium_DEG_summary.csv"),
          row.names = FALSE)

cat("\nSaved: data/processed/syn_degs_filtered.rds\n")
cat("Saved: tables/table_F1_synovium_DEG_summary.csv\n\n")

cat("============================================\n")
cat("STAGE 2 COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")