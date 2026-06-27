# =============================================================================
# 17_cross_tissue_comparison_figure.R
# -----------------------------------------------------------------------------
# Figure F1: Cross-tissue comparison (cartilage vs synovium)
#   Panel A: Cargo-OA enrichment observed vs expected
#   Panel B: Directional cargo targeting (% of OA-up vs OA-down targeted)
#
# INPUT  : data/processed/intersection_results.rds (cartilage)
#          data/processed/syn_intersection_results.rds (synovium)
# OUTPUT : figures/fig_F1_cross_tissue_comparison.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# -----------------------------------------------------------------------------
# Panel A — observed vs expected
# -----------------------------------------------------------------------------
panel_a_df <- data.frame(
  Tissue = rep(c("Cartilage\n(GSE114007)", "Synovium\n(GSE55235)"), each = 2),
  Type   = rep(c("Expected", "Observed"), 2),
  Value  = c(2968.1, 2918,
             1314.1, 1394),
  stringsAsFactors = FALSE
)
panel_a_df$Tissue <- factor(panel_a_df$Tissue,
                            levels = c("Cartilage\n(GSE114007)", "Synovium\n(GSE55235)"))
panel_a_df$Type   <- factor(panel_a_df$Type, levels = c("Expected", "Observed"))

p_a <- ggplot(panel_a_df, aes(x = Tissue, y = Value, fill = Type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.62,
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", Value)),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 3.3, fontface = "bold", color = "grey20") +
  annotate("text", x = 1, y = 3150, label = "p = 0.98 (ns)",
           size = 3.8, fontface = "italic", color = "grey30") +
  annotate("text", x = 2, y = 1540,
           label = "p = 1.0 x 10^-7",
           size = 3.8, fontface = "italic", color = "#8B0000") +
  scale_fill_manual(values = c("Expected" = "#BBBBBB",
                               "Observed" = "#2166AC"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A  Cargo-target / OA-DEG overlap",
       subtitle = "Observed vs expected under null",
       x = NULL, y = "Number of genes") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x   = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top"
  )

# -----------------------------------------------------------------------------
# Panel B — directional targeting
# -----------------------------------------------------------------------------
panel_b_df <- data.frame(
  Tissue = rep(c("Cartilage\n(GSE114007)", "Synovium\n(GSE55235)"), each = 2),
  Direction = rep(c("OA-up genes", "OA-down genes"), 2),
  Pct_targeted = c(70.6, 67.4,
                   78.6, 90.7),
  stringsAsFactors = FALSE
)
panel_b_df$Tissue    <- factor(panel_b_df$Tissue,
                               levels = c("Cartilage\n(GSE114007)", "Synovium\n(GSE55235)"))
panel_b_df$Direction <- factor(panel_b_df$Direction,
                               levels = c("OA-up genes", "OA-down genes"))

bg_cart <- 70.2
bg_syn  <- 78.6

p_b <- ggplot(panel_b_df, aes(x = Tissue, y = Pct_targeted, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.62,
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Pct_targeted)),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 3.3, fontface = "bold", color = "grey20") +
  annotate("segment", x = 0.55, xend = 1.45, y = bg_cart, yend = bg_cart,
           color = "grey45", linetype = "dashed", linewidth = 0.5) +
  annotate("segment", x = 1.55, xend = 2.45, y = bg_syn, yend = bg_syn,
           color = "grey45", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 1, y = bg_cart - 3,
           label = "70.2% (background)",
           size = 2.9, color = "grey40", fontface = "italic") + +
  annotate("text", x = 2, y = bg_syn - 3,
           label = "78.6% (background)",
           size = 2.9, color = "grey40", fontface = "italic") + +
  annotate("text", x = 1, y = 80, label = "p = 0.026",
           size = 3.5, fontface = "italic", color = "grey30") +
  annotate("text", x = 2, y = 100,
           label = "p = 3.1 x 10^-11",
           size = 3.5, fontface = "italic", color = "#8B0000") +
  scale_fill_manual(values = c("OA-up genes"   = "#E41A1C",
                               "OA-down genes" = "#377EB8"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "B  Directional cargo targeting",
       subtitle = "% of OA-up and OA-down genes that are cargo targets",
       x = NULL, y = "% of genes targeted") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x   = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top"
  )

# -----------------------------------------------------------------------------
# Combine and save
# -----------------------------------------------------------------------------
combined <- p_a + p_b +
  plot_layout(ncol = 2) +
  plot_annotation(
    title    = "Cross-tissue cargo-target enrichment: cartilage vs synovium",
    subtitle = "Cargo targets validated in miRTarBase (tier-1 strict); background: tissue-expressed genes that are miRNA-targetable",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, color = "grey30"))
  )

ggsave(file.path(PROJECT_ROOT, "figures/fig_F1_cross_tissue_comparison.pdf"),
       combined, width = 12, height = 5.5)

cat("Saved: figures/fig_F1_cross_tissue_comparison.pdf\n\n")

cat("=== FIGURE F1 SUMMARY ===\n")
cat("Panel A -- Overlap observed vs expected:\n")
cat("  Cartilage: expected 2968, observed 2918 (ratio 0.98, p = 0.98)\n")
cat("  Synovium:  expected 1314, observed 1394 (ratio 1.06, p = 1.0e-7)\n\n")
cat("Panel B -- Directional cargo targeting:\n")
cat("  Cartilage: OA-up 70.6%, OA-down 67.4% (Fisher OR 1.16, p = 0.026)\n")
cat("  Synovium:  OA-up 78.6%, OA-down 90.7% (Fisher OR 0.38, p = 3.1e-11)\n\n")
cat("Panel B interpretation:\n")
cat("  Cartilage marginally prefers OA-up targeting\n")
cat("  Synovium strongly prefers OA-down targeting (opposite direction)\n")