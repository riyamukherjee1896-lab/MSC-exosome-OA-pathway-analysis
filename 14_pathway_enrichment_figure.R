# =============================================================================
# 14_pathway_enrichment_figure.R
# -----------------------------------------------------------------------------
# Figure D3: Pathway-level cargo targeting coverage vs background rate
#
# Shows the three pre-specified anchor pathways: coverage of cargo targeting
# within each set, compared to the cartilage-expressed background rate.
#
# Rationale for coverage % over OR: two pathways have 100% coverage, giving
# OR = Inf which cannot be plotted. Coverage % is more interpretable and
# scientifically equivalent.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# -----------------------------------------------------------------------------
# Load Block D pathway results
# -----------------------------------------------------------------------------
pw <- read.csv(
  file.path(PROJECT_ROOT, "tables/table_D3_pathway_enrichment.csv"),
  stringsAsFactors = FALSE)

cat("=== PATHWAY DATA FROM BLOCK D ===\n")
print(pw, row.names = FALSE)
cat("\n")

# Background cargo-targeting rate from Block D:
# 10,844 cargo targets / 15,440 cartilage-expressed = 70.2%
bg_rate_pct <- 70.2

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
pw_plot <- pw %>%
  dplyr::mutate(
    Pathway = dplyr::recode(Pathway,
                            "Antioxidant"      = "Antioxidant",
                            "Ferroptosis"      = "Ferroptosis",
                            "miRNA_Biogenesis" = "miRNA biogenesis"),
    Pathway = factor(Pathway,
                     levels = c("Antioxidant", "Ferroptosis", "miRNA biogenesis")),
    # Build human-readable p-value label
    p_label = sprintf("p = %.3f", P_value),
    # Make a descriptive label: "X of Y genes (Z%)"
    count_label = sprintf("%d of %d genes (%.0f%%)",
                          Cargo_targeted, Pathway_size_in_bg, Pct_targeted)
  )

cat("=== PLOT-READY DATA ===\n")
print(pw_plot[, c("Pathway", "Pathway_size_in_bg", "Cargo_targeted",
                  "Pct_targeted", "P_value")], row.names = FALSE)
cat("\n")

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
pathway_colors <- c(
  "Antioxidant"      = "#FF7F00",
  "Ferroptosis"      = "#377EB8",
  "miRNA biogenesis" = "#984EA3"
)

p <- ggplot(pw_plot, aes(x = Pct_targeted, y = Pathway, fill = Pathway)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_vline(xintercept = bg_rate_pct,
             linetype = "dashed", color = "grey35", linewidth = 0.6) +
  annotate("text",
           x = bg_rate_pct, y = 3.55,
           label = sprintf("Background rate\n(%.1f%%)", bg_rate_pct),
           hjust = 1.05, size = 3.2, color = "grey30", fontface = "italic") +
  # Coverage label inside the bar (or just after if close to edge)
  geom_text(aes(label = count_label),
            hjust = 1.1, size = 3.3, fontface = "bold", color = "white") +
  # p-value label outside the bar
  geom_text(aes(label = p_label, x = Pct_targeted + 2),
            hjust = 0, size = 3.5, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = pathway_colors) +
  scale_x_continuous(
    limits = c(0, 118),
    breaks = c(0, 25, 50, 70.2, 100),
    labels = c("0", "25", "50", "70", "100"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title    = "Anchor-pathway enrichment in cargo-target set",
    subtitle = paste0("Cargo coverage of pre-specified OA anchor gene sets\n",
                      "vs. 70.2% background cargo-targeting rate ",
                      "(Fisher's exact test, one-sided)"),
    x        = "% of pathway genes that are cargo targets",
    y        = NULL
  ) +
  guides(fill = "none") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    axis.text.y      = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(size = 9),
    plot.margin      = margin(10, 20, 15, 10)
  )

ggsave(file.path(PROJECT_ROOT, "figures/fig_D3_pathway_enrichment.pdf"),
       p, width = 9, height = 4.5)

cat("Saved: figures/fig_D3_pathway_enrichment.pdf\n\n")

cat("=== INTERPRETATION ===\n")
cat(sprintf("All three pathways show coverage above the %.1f%% background rate.\n",
            bg_rate_pct))
cat("Statistical significance (Fisher's exact):\n")
for (i in seq_len(nrow(pw_plot))) {
  cat(sprintf("  %-18s %.0f%% coverage, p = %.3f\n",
              pw_plot$Pathway[i],
              pw_plot$Pct_targeted[i],
              pw_plot$P_value[i]))
}