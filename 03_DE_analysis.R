# =============================================================================
# 03_DE_analysis.R
# -----------------------------------------------------------------------------
# Paper 3 — Differential expression with batch adjustment
#
# INPUT  : data/raw/GSE114007_raw_counts.csv
#          data/processed/metadata_cleaned.csv
#
# OUTPUT : data/processed/dds_oa.rds
#          data/processed/vst_matrix.rds
#          data/processed/de_results_full.csv
#          data/processed/de_results_annotated.csv
#          tables/table_S1_DE_key_genes.csv
#          tables/table_DE_comparison_legacy_vs_batch.csv
#          logs/03_DE_analysis_sessionInfo.txt
#
# DESIGN : ~ Batch + Condition, Control as reference
# FILTER : rowSums(counts >= 10) >= 5
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("Paper 3 — STEP 3: DE Analysis (batch-adjusted)\n")
cat("============================================\n")
cat("Run started:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load raw counts
# -----------------------------------------------------------------------------
counts_raw <- read.csv(
  file.path(PROJECT_ROOT, "data/raw/GSE114007_raw_counts.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cat("Raw counts loaded:", nrow(counts_raw), "genes x",
    ncol(counts_raw) - 1, "samples\n")

# Separate gene column from count matrix
count_matrix <- as.matrix(counts_raw[, -1])
rownames(count_matrix) <- counts_raw$symbol
mode(count_matrix) <- "integer"   # DESeq2 requires integers

# -----------------------------------------------------------------------------
# 2. Load cleaned metadata and align with counts
# -----------------------------------------------------------------------------
meta <- read.csv(
  file.path(PROJECT_ROOT, "data/processed/metadata_cleaned.csv"),
  stringsAsFactors = FALSE
)

# Enforce factor levels
meta$Condition <- factor(meta$Condition, levels = c("Control", "OA"))
meta$Batch     <- factor(meta$Batch,     levels = c("Cohort_A", "Cohort_B"))

# Reorder count matrix to match metadata row order (critical for DESeq2)
count_matrix <- count_matrix[, meta$SampleID]

# Set rownames of metadata to sample IDs
rownames(meta) <- meta$SampleID

# Verify alignment
stopifnot(all(colnames(count_matrix) == rownames(meta)))
cat("Metadata aligned. Sample order verified.\n\n")

cat("=== BATCH x CONDITION (design matrix sanity) ===\n")
print(table(meta$Batch, meta$Condition))
cat("\n")

# -----------------------------------------------------------------------------
# 3. Build DESeq2 object
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = meta,
  design    = ~ Batch + Condition
)
cat("DESeq2 object built.\n")
cat("  Design:", deparse(design(dds)), "\n")
cat("  Reference for Condition:", levels(dds$Condition)[1], "\n")
cat("  Reference for Batch    :", levels(dds$Batch)[1], "\n\n")

# -----------------------------------------------------------------------------
# 4. Filter low-expressed genes
# -----------------------------------------------------------------------------
# Keep genes with at least 10 counts in at least 5 samples
keep <- rowSums(counts(dds) >= 10) >= 5
cat("=== GENE FILTER ===\n")
cat("Before filter:", nrow(dds), "genes\n")
dds <- dds[keep, ]
cat("After filter :", nrow(dds), "genes\n")
cat("Removed      :", sum(!keep), "genes\n\n")

# Check specifically: does AIFM2 survive?
cat("AIFM2 in filtered set:", "AIFM2" %in% rownames(dds), "\n")
if ("AIFM2" %in% rownames(dds)) {
  cat("  AIFM2 raw counts:\n")
  print(counts(dds)["AIFM2", ])
}
cat("\n")

# -----------------------------------------------------------------------------
# 5. Run DESeq2
# -----------------------------------------------------------------------------
cat("Running DESeq2 (this takes ~30–60 seconds)...\n")
dds <- DESeq(dds)
cat("DESeq2 complete.\n\n")

# Extract Condition effect (OA vs Control, adjusted for Batch)
res <- results(dds, name = "Condition_OA_vs_Control")

cat("=== DE RESULTS SUMMARY ===\n")
summary(res)
cat("\n")

# -----------------------------------------------------------------------------
# 6. Build annotated results table
# -----------------------------------------------------------------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("Gene") %>%
  mutate(
    Significant = padj < 0.05 & !is.na(padj),
    Direction = case_when(
      Significant & log2FoldChange > 0 ~ "Up in OA",
      Significant & log2FoldChange < 0 ~ "Down in OA",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

n_sig      <- sum(res_df$Significant, na.rm = TRUE)
n_up       <- sum(res_df$Direction == "Up in OA")
n_down     <- sum(res_df$Direction == "Down in OA")
cat(sprintf("Significant (padj<0.05): %d  (%d up, %d down)\n\n",
            n_sig, n_up, n_down))

# -----------------------------------------------------------------------------
# 7. Save DE object and results
# -----------------------------------------------------------------------------
saveRDS(dds,
        file.path(PROJECT_ROOT, "data/processed/dds_oa.rds"))

write.csv(res_df,
          file.path(PROJECT_ROOT, "data/processed/de_results_full.csv"),
          row.names = FALSE)

write.csv(res_df,
          file.path(PROJECT_ROOT, "data/processed/de_results_annotated.csv"),
          row.names = FALSE)

cat("Saved: dds_oa.rds\n")
cat("Saved: de_results_full.csv\n")
cat("Saved: de_results_annotated.csv\n\n")

# -----------------------------------------------------------------------------
# 8. VST transformation (for downstream WGCNA/correlation/ROC)
# -----------------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

saveRDS(vst_mat,
        file.path(PROJECT_ROOT, "data/processed/vst_matrix.rds"))

cat("VST transformation complete.\n")
cat("Saved: vst_matrix.rds (", nrow(vst_mat), "genes x",
    ncol(vst_mat), "samples)\n\n")

# -----------------------------------------------------------------------------
# 9. Key genes table for manuscript
# -----------------------------------------------------------------------------
key_genes <- c(
  "ALDH2", "SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX3", "GPX4",
  "NQO1", "HMOX1", "NFE2L2", "KEAP1", "TXNRD1", "TXN", "PRDX1",
  "GSTA1", "GSTP1", "AKR1B1", "PTGES3",
  "DROSHA", "DGCR8", "DICER1", "AGO1", "AGO2", "TARBP2", "XRN1",
  "DDX5", "DDX17", "TNRC6A", "HNRNPA1", "LIN28A", "KHSRP", "ADAR",
  "QKI", "SRSF1", "SRSF3",
  "SLC7A11", "SLC3A2", "ACSL4", "LPCAT3", "TFRC", "FTH1", "FTL",
  "NCOA4", "STEAP3", "VDAC2", "CISD1", "NFS1", "AIFM2"
)

key_de <- res_df %>% filter(Gene %in% key_genes)

write.csv(key_de,
          file.path(PROJECT_ROOT, "tables/table_S1_DE_key_genes.csv"),
          row.names = FALSE)
cat("Saved: tables/table_S1_DE_key_genes.csv (", nrow(key_de), "genes)\n\n")

# -----------------------------------------------------------------------------
# 10. Compare new (batch-adjusted) vs legacy DE for headline genes
# -----------------------------------------------------------------------------
legacy <- read.csv(
  file.path(PROJECT_ROOT, "data/raw/GSE114007_Deseq_file.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

headline <- c("ALDH2", "GPX3", "SOD2", "GPX4", "NQO1", "DDX17",
              "AGO2", "DICER1", "DROSHA", "DGCR8", "AIFM2")

comparison <- data.frame(Gene = headline) %>%
  left_join(
    legacy %>% select(Gene, legacy_log2FC = log2FoldChange,
                      legacy_padj   = padj),
    by = "Gene"
  ) %>%
  left_join(
    res_df %>% select(Gene, new_log2FC = log2FoldChange,
                      new_padj   = padj),
    by = "Gene"
  ) %>%
  mutate(
    delta_log2FC = new_log2FC - legacy_log2FC,
    legacy_sig   = legacy_padj < 0.05,
    new_sig      = new_padj < 0.05,
    sig_change   = case_when(
      is.na(legacy_sig) | is.na(new_sig) ~ "NA",
      legacy_sig & new_sig               ~ "Both sig",
      !legacy_sig & new_sig              ~ "Newly sig",
      legacy_sig & !new_sig              ~ "Lost sig",
      TRUE                               ~ "Neither"
    )
  )

write.csv(comparison,
          file.path(PROJECT_ROOT, "tables/table_DE_comparison_legacy_vs_batch.csv"),
          row.names = FALSE)

cat("=== LEGACY vs BATCH-ADJUSTED DE (key genes) ===\n")
print(comparison, row.names = FALSE, digits = 3)
cat("\nSaved: tables/table_DE_comparison_legacy_vs_batch.csv\n\n")

# -----------------------------------------------------------------------------
# 11. Save session info for reproducibility
# -----------------------------------------------------------------------------
sink(file.path(PROJECT_ROOT, "logs/03_DE_analysis_sessionInfo.txt"))
cat("Script: 03_DE_analysis.R\n")
cat("Run:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("============================================\n")
cat("STEP 3 COMPLETE\n")
cat("============================================\n")