# =============================================================================
# 07_OA_DEG_refinement.R
# -----------------------------------------------------------------------------
# Block B: Refine OA-dysregulated gene list for target intersection
#
# INPUT  : data/processed/de_results_annotated.csv
#
# OUTPUT : data/processed/oa_degs_filtered.rds
#          tables/table_B1_OA_DEG_summary.csv
#          tables/table_B2_DEG_sensitivity.csv
#          tables/table_B3_anchor_genes_DE.csv
#          figures/fig_B1_volcano.pdf
#          figures/fig_B2_anchor_forest.pdf
#          logs/07_OA_DEG_sessionInfo.txt
#
# NOTES  :
#   - Dual threshold: padj < 0.05 AND |log2FC| >= 0.585 (1.5-fold)
#   - GPX4 belongs to both Antioxidant and Ferroptosis pathway sets.
#     For the forest plot (display only), it is deduplicated under
#     Ferroptosis (its more current OA-relevant role). Full pathway
#     memberships are preserved in the saved `anchor_results` object.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("BLOCK B: OA DEG Refinement\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load DE results
# -----------------------------------------------------------------------------
de <- read.csv(file.path(PROJECT_ROOT, "data/processed/de_results_annotated.csv"),
               stringsAsFactors = FALSE)

cat("=== INPUT DE RESULTS ===\n")
cat("Genes loaded:", nrow(de), "\n")
cat("Columns:", paste(names(de), collapse = ", "), "\n")

required_cols <- c("Gene", "log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, names(de))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
cat("All required columns present.\n\n")

n_na_padj <- sum(is.na(de$padj))
cat("Genes with NA padj (low-count filtered by DESeq2):", n_na_padj, "\n")
cat("Genes with valid padj:", sum(!is.na(de$padj)), "\n\n")

# -----------------------------------------------------------------------------
# 2. Apply dual threshold
# -----------------------------------------------------------------------------
LOG2FC_THRESHOLD <- 0.585  # 1.5-fold change
PADJ_THRESHOLD   <- 0.05

oa_degs <- de %>%
  filter(!is.na(padj),
         padj < PADJ_THRESHOLD,
         abs(log2FoldChange) >= LOG2FC_THRESHOLD) %>%
  mutate(
    Direction    = ifelse(log2FoldChange > 0, "Up_in_OA", "Down_in_OA"),
    FC_magnitude = 2^abs(log2FoldChange)
  ) %>%
  arrange(padj)

oa_up   <- oa_degs %>% filter(Direction == "Up_in_OA")
oa_down <- oa_degs %>% filter(Direction == "Down_in_OA")

cat("=== PRIMARY OA DEG SET ===\n")
cat(sprintf("Threshold: padj < %.2f AND |log2FC| >= %.3f (>= %.1f-fold)\n",
            PADJ_THRESHOLD, LOG2FC_THRESHOLD, 2^LOG2FC_THRESHOLD))
cat(sprintf("OA DEGs total: %d\n",   nrow(oa_degs)))
cat(sprintf("  Up in OA:    %d\n",   nrow(oa_up)))
cat(sprintf("  Down in OA:  %d\n\n", nrow(oa_down)))

up_down_ratio <- nrow(oa_up) / nrow(oa_down)
cat(sprintf("Up/Down ratio: %.2f  (expect roughly balanced for OA)\n\n",
            up_down_ratio))

# -----------------------------------------------------------------------------
# 3. Anchor gene inspection
# -----------------------------------------------------------------------------
anchor_genes <- list(
  Antioxidant = c("ALDH2","SOD1","SOD2","SOD3","CAT","GPX1","GPX3","GPX4",
                  "NQO1","HMOX1","NFE2L2","KEAP1","TXNRD1","TXN","PRDX1",
                  "GSTA1","GSTP1","AKR1B1","PTGES3"),
  Ferroptosis = c("GPX4","SLC7A11","SLC3A2","ACSL4","LPCAT3","TFRC","FTH1",
                  "FTL","NCOA4","STEAP3","VDAC2","CISD1","NFS1","AIFM2"),
  miRNA_Biogenesis = c("DROSHA","DGCR8","DICER1","AGO1","AGO2","TARBP2","XRN1",
                       "DDX5","DDX17","TNRC6A","HNRNPA1","LIN28A","KHSRP","ADAR",
                       "QKI","SRSF1","SRSF3")
)

anchor_results <- data.frame()
cat("=== ANCHOR GENES IN PRIMARY DEG SET ===\n")
for (set_name in names(anchor_genes)) {
  genes_in_set <- anchor_genes[[set_name]]
  in_de        <- de %>% filter(Gene %in% genes_in_set)
  hits         <- oa_degs %>% filter(Gene %in% genes_in_set)
  cat(sprintf("\n%s: %d/%d in DE table, %d/%d pass threshold\n",
              set_name,
              nrow(in_de), length(genes_in_set),
              nrow(hits), length(genes_in_set)))
  
  if (nrow(hits) > 0) {
    display <- hits[, c("Gene", "log2FoldChange", "padj", "Direction")]
    display$log2FoldChange <- round(display$log2FoldChange, 3)
    display$padj <- formatC(display$padj, format = "e", digits = 2)
    print(display, row.names = FALSE)
    anchor_results <- rbind(anchor_results,
                            data.frame(Gene      = hits$Gene,
                                       Pathway   = set_name,
                                       log2FC    = hits$log2FoldChange,
                                       padj      = hits$padj,
                                       Direction = hits$Direction))
  }
  
  not_passing <- in_de %>% filter(!(Gene %in% hits$Gene))
  if (nrow(not_passing) > 0) {
    cat("  Anchor genes in DE table but not passing FC threshold:\n")
    display_not <- not_passing[, c("Gene", "log2FoldChange", "padj")]
    display_not$log2FoldChange <- round(display_not$log2FoldChange, 3)
    display_not$padj <- formatC(display_not$padj, format = "e", digits = 2)
    print(display_not, row.names = FALSE)
  }
}
cat("\n")

write.csv(anchor_results,
          file.path(PROJECT_ROOT, "tables/table_B3_anchor_genes_DE.csv"),
          row.names = FALSE)
cat("Saved: tables/table_B3_anchor_genes_DE.csv\n\n")

# -----------------------------------------------------------------------------
# 4. Sensitivity analysis
# -----------------------------------------------------------------------------
cat("=== SENSITIVITY TO log2FC CUTOFF ===\n")
fc_thresholds <- c(0, 0.263, 0.585, 1, 1.585)
fc_labels     <- c("No FC filter (1x)", "1.2-fold (0.263)",
                   "1.5-fold (0.585) [PRIMARY]", "2-fold (1.0)", "3-fold (1.585)")

sens_tbl <- data.frame(
  Threshold     = fc_labels,
  log2FC_cutoff = fc_thresholds,
  N_up   = sapply(fc_thresholds, function(k)
    sum(de$padj < 0.05 & de$log2FoldChange >=  k, na.rm = TRUE)),
  N_down = sapply(fc_thresholds, function(k)
    sum(de$padj < 0.05 & de$log2FoldChange <= -k, na.rm = TRUE))
)
sens_tbl$N_total <- sens_tbl$N_up + sens_tbl$N_down
print(sens_tbl, row.names = FALSE)
cat("\n")

write.csv(sens_tbl,
          file.path(PROJECT_ROOT, "tables/table_B2_DEG_sensitivity.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Summary table
# -----------------------------------------------------------------------------
summary_tbl <- data.frame(
  Category = c("Total_genes_in_DE",
               "Genes_with_valid_padj",
               "Significant_padj_only",
               "Significant_with_1.5FC",
               "Up_in_OA",
               "Down_in_OA",
               "Anchor_Antioxidant_DE",
               "Anchor_Ferroptosis_DE",
               "Anchor_miRNA_biogenesis_DE"),
  Count    = c(nrow(de),
               sum(!is.na(de$padj)),
               sum(de$padj < 0.05, na.rm = TRUE),
               nrow(oa_degs),
               nrow(oa_up),
               nrow(oa_down),
               sum(oa_degs$Gene %in% anchor_genes$Antioxidant),
               sum(oa_degs$Gene %in% anchor_genes$Ferroptosis),
               sum(oa_degs$Gene %in% anchor_genes$miRNA_Biogenesis))
)

cat("=== SUMMARY ===\n")
print(summary_tbl, row.names = FALSE)
cat("\n")

write.csv(summary_tbl,
          file.path(PROJECT_ROOT, "tables/table_B1_OA_DEG_summary.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Save refined DEG objects (MOVED BEFORE FIGURES so save never depends on
#    figure generation succeeding)
# -----------------------------------------------------------------------------
refined <- list(
  oa_degs          = oa_degs,
  oa_up            = oa_up,
  oa_down          = oa_down,
  anchor_results   = anchor_results,
  threshold_padj   = PADJ_THRESHOLD,
  threshold_log2FC = LOG2FC_THRESHOLD,
  anchor_gene_sets = anchor_genes
)
saveRDS(refined,
        file.path(PROJECT_ROOT, "data/processed/oa_degs_filtered.rds"))
cat("Saved: data/processed/oa_degs_filtered.rds\n\n")

# -----------------------------------------------------------------------------
# 7. Figure B1 — Volcano plot
# -----------------------------------------------------------------------------
de_plot <- de %>%
  filter(!is.na(padj)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    category = case_when(
      padj < 0.05 & log2FoldChange >=  LOG2FC_THRESHOLD ~ "Up in OA",
      padj < 0.05 & log2FoldChange <= -LOG2FC_THRESHOLD ~ "Down in OA",
      TRUE ~ "Not significant"
    )
  )

de_plot$neg_log10_padj_capped <- pmin(de_plot$neg_log10_padj, 20)
de_plot$log2FC_capped          <- pmax(pmin(de_plot$log2FoldChange, 6), -6)

all_anchors <- unique(unlist(anchor_genes))
de_plot$label <- ifelse(de_plot$Gene %in% all_anchors &
                          de_plot$category != "Not significant",
                        de_plot$Gene, NA)

p_volcano <- ggplot(de_plot, aes(x = log2FC_capped, y = neg_log10_padj_capped)) +
  geom_point(aes(color = category), size = 0.8, alpha = 0.6) +
  geom_point(data = de_plot %>% filter(!is.na(label)),
             color = "black", size = 2.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  box.padding = 0.5, segment.size = 0.3) +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Up in OA"        = "#E41A1C",
                                "Down in OA"      = "#377EB8",
                                "Not significant" = "grey70")) +
  labs(title    = "OA cartilage differentially expressed genes",
       subtitle = sprintf("padj < 0.05 and |log2FC| >= %.3f (>= %.1f-fold)",
                          LOG2FC_THRESHOLD, 2^LOG2FC_THRESHOLD),
       x        = "log2 fold change (OA vs Control)",
       y        = "-log10 adjusted p-value",
       color    = "") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(PROJECT_ROOT, "figures/fig_B1_volcano.pdf"),
       p_volcano, width = 9, height = 7)
cat("Saved: figures/fig_B1_volcano.pdf\n")

# -----------------------------------------------------------------------------
# 8. Figure B2 — Anchor gene forest plot (DEDUPLICATED)
# -----------------------------------------------------------------------------
# GPX4 is in both Antioxidant and Ferroptosis lists, so anchor_results
# has it twice. For display, deduplicate: assign GPX4 to Ferroptosis.
# The original anchor_results object (with both memberships) is preserved in
# the saved RDS; this dedup only affects the forest plot.

if (nrow(anchor_results) > 0) {
  
  forest_df <- anchor_results %>%
    mutate(pathway_priority = case_when(
      Gene == "GPX4" & Pathway == "Ferroptosis" ~ 1,
      Gene == "GPX4" & Pathway == "Antioxidant" ~ 2,
      TRUE                                      ~ 1
    )) %>%
    arrange(Gene, pathway_priority) %>%
    group_by(Gene) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-pathway_priority) %>%
    arrange(log2FC) %>%
    as.data.frame()
  
  dup_note <- ""
  if (any(duplicated(forest_df$Gene))) {
    stop("Forest plot: unexpected duplicated genes after dedup.")
  }
  if ("GPX4" %in% anchor_results$Gene[anchor_results$Pathway == "Antioxidant"] &&
      "GPX4" %in% anchor_results$Gene[anchor_results$Pathway == "Ferroptosis"]) {
    dup_note <- "\nNote: GPX4 is a member of both sets; displayed under Ferroptosis."
  }
  
  forest_df$Gene <- factor(forest_df$Gene, levels = forest_df$Gene)
  
  p_forest <- ggplot(forest_df, aes(x = log2FC, y = Gene, color = Pathway)) +
    geom_point(size = 3.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
               linetype = "dotted", color = "grey70") +
    scale_color_manual(values = c("Antioxidant"      = "#FF7F00",
                                  "Ferroptosis"      = "#377EB8",
                                  "miRNA_Biogenesis" = "#984EA3")) +
    labs(title    = "Anchor-pathway genes dysregulated in OA cartilage",
         subtitle = paste0("padj < 0.05 and |log2FC| >= 0.585 (1.5-fold)",
                           dup_note),
         x        = "log2 fold change (OA vs Control)",
         y        = "",
         color    = "Pathway") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  ggsave(file.path(PROJECT_ROOT, "figures/fig_B2_anchor_forest.pdf"),
         p_forest, width = 8, height = 6)
  cat("Saved: figures/fig_B2_anchor_forest.pdf\n")
  cat("  (GPX4 deduplicated for display; full pathway memberships preserved in",
      "anchor_results)\n")
}

# -----------------------------------------------------------------------------
# 9. File existence verification
# -----------------------------------------------------------------------------
cat("\n=== FILE EXISTENCE VERIFICATION ===\n")
files_to_check <- c(
  "data/processed/oa_degs_filtered.rds",
  "tables/table_B1_OA_DEG_summary.csv",
  "tables/table_B2_DEG_sensitivity.csv",
  "tables/table_B3_anchor_genes_DE.csv",
  "figures/fig_B1_volcano.pdf",
  "figures/fig_B2_anchor_forest.pdf"
)
for (f in files_to_check) {
  full <- file.path(PROJECT_ROOT, f)
  cat(sprintf("  %-50s  %s\n", f, ifelse(file.exists(full), "OK", "MISSING")))
}

# -----------------------------------------------------------------------------
# 10. Session info
# -----------------------------------------------------------------------------
sink(file.path(PROJECT_ROOT, "logs/07_OA_DEG_sessionInfo.txt"))
cat("Script: 07_OA_DEG_refinement.R\n")
cat("Run:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\n============================================\n")
cat("BLOCK B COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")