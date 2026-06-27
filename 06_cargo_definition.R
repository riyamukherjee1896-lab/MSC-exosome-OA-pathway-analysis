PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
})

cargo_data <- readRDS(file.path(PROJECT_ROOT, "data/processed/cargo_miRNAs.rds"))
raw_mat      <- cargo_data$raw_mat
norm_mat     <- cargo_data$norm_mat
cargo_mirnas <- cargo_data$cargo_mirnas
sample_meta  <- cargo_data$sample_meta

cat("============================================\n")
cat("BLOCK A: ABUNDANCE CORRECTION\n")
cat("============================================\n\n")

# -----------------------------------------------------------------------------
# 1. Proper abundance ranking
# -----------------------------------------------------------------------------
# For each cargo miRNA, compute:
#   - mean raw count (treating NA as 0, since NA = not detected)
#   - median raw count (same treatment)
#   - n samples detected

raw_zeros <- raw_mat[cargo_mirnas, ]
raw_zeros[is.na(raw_zeros)] <- 0   # undetected -> 0 reads

cargo_abundance <- data.frame(
  miRNA            = cargo_mirnas,
  mean_raw_counts  = rowMeans(raw_zeros),
  median_raw_counts= apply(raw_zeros, 1, median),
  max_raw_counts   = apply(raw_zeros, 1, max),
  n_detected       = rowSums(raw_mat[cargo_mirnas, ] > 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)
cargo_abundance <- cargo_abundance[order(-cargo_abundance$mean_raw_counts), ]

cat("=== TRUE TOP-15 CARGO miRNAs (by mean raw count, NA-as-zero) ===\n")
print(head(cargo_abundance, 15), row.names = FALSE, digits = 1)

cat("\n=== BOTTOM-5 CARGO miRNAs (lowest mean raw count) ===\n")
print(tail(cargo_abundance, 5), row.names = FALSE, digits = 2)

# Save corrected abundance table
write.csv(cargo_abundance,
          file.path(PROJECT_ROOT, "tables/table_A1_cargo_abundance_corrected.csv"),
          row.names = FALSE)
cat("\nSaved: tables/table_A1_cargo_abundance_corrected.csv\n\n")

# -----------------------------------------------------------------------------
# 2. Regenerate Figure A3 with proper top-30
# -----------------------------------------------------------------------------
top30_corrected <- head(cargo_abundance$miRNA, 30)

# Build heatmap matrix on log2-transformed raw counts (NA -> 0 -> log2(1) = 0)
heatmap_mat <- raw_mat[top30_corrected, ]
heatmap_mat[is.na(heatmap_mat)] <- 0
heatmap_mat <- log2(heatmap_mat + 1)   # log2(x+1) handles zeros

col_annot <- data.frame(Source = sample_meta$Source,
                        row.names = sample_meta$Sample)
annot_colors <- list(
  Source = c(Adipose = "#E41A1C",
             BoneMarrow = "#377EB8",
             WhartonJelly = "#4DAF4A")
)

pdf(file.path(PROJECT_ROOT, "figures/fig_A3_top30_cargo_heatmap.pdf"),
    width = 9, height = 8)
pheatmap(heatmap_mat,
         annotation_col = col_annot,
         annotation_colors = annot_colors,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("white", "#FFF5A0", "#FB8E4F", "#C8102E"))(100),
         main = "Top-30 cargo miRNAs by mean raw count (log2+1 scale)",
         fontsize_row = 8,
         fontsize_col = 9,
         display_numbers = FALSE,
         show_colnames = TRUE)
dev.off()

cat("Saved: figures/fig_A3_top30_cargo_heatmap.pdf (CORRECTED)\n\n")

# -----------------------------------------------------------------------------
# 3. Verify cargo definition unchanged
# -----------------------------------------------------------------------------
cat("=== INTEGRITY CHECK ===\n")
cat("Cargo size (>= 6/12 detection):", length(cargo_mirnas),
    "  (should be 276)\n")
cat("Total miRNAs in data:", nrow(raw_mat),
    "  (should be 1869)\n")
cat("Duplicated rownames:", sum(duplicated(rownames(raw_mat))),
    "  (should be 0)\n")

# -----------------------------------------------------------------------------
# 4. Save updated cargo_data with abundance table
# -----------------------------------------------------------------------------
cargo_data$cargo_abundance <- cargo_abundance
saveRDS(cargo_data, file.path(PROJECT_ROOT, "data/processed/cargo_miRNAs.rds"))
cat("Saved: cargo_miRNAs.rds (with corrected abundance table)\n\n")

cat("============================================\n")
cat("ABUNDANCE CORRECTION COMPLETE\n")
cat("============================================\n")