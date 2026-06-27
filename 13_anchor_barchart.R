# =============================================================================
# 13_anchor_barchart.R
# -----------------------------------------------------------------------------
# Figure D2: Anchor gene cargo-targeting bar chart
#
# Shows number of cargo mature miRNAs targeting each anchor OA gene,
# colored by pre-specified pathway (antioxidant / ferroptosis / miRNA biogenesis).
#
# GPX4 is a member of BOTH antioxidant and ferroptosis sets. For display,
# assigned to Ferroptosis (consistent with Block B forest plot).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
targets <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/mirna_targets_validated.rds"))
oa <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/oa_degs_filtered.rds"))

# Use mature-miRNA-level target frequency if available (from Block C audit fix)
# Otherwise recompute
if (!is.null(targets$target_frequency_mature)) {
  tf <- targets$target_frequency_mature
} else {
  tf <- targets$cargo_edges %>%
    dplyr::group_by(target_symbol) %>%
    dplyr::summarise(
      n_mature_miRNAs_targeting = dplyr::n_distinct(mature_mirna_id),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# Define anchor-pathway assignments (GPX4 -> Ferroptosis for display)
# -----------------------------------------------------------------------------
anchor_df <- data.frame(
  Gene = c(
    # Antioxidant
    "ALDH2", "SOD2", "GPX3", "NQO1", "HMOX1", "NFE2L2", "PRDX1",
    # Ferroptosis
    "GPX4", "ACSL4", "SLC3A2", "TFRC",
    # miRNA biogenesis
    "AGO2", "DICER1", "DDX17", "DGCR8", "DROSHA", "SRSF3"
  ),
  Pathway = c(
    rep("Antioxidant",      7),
    rep("Ferroptosis",       4),
    rep("miRNA biogenesis",  6)
  ),
  stringsAsFactors = FALSE
)

# Join with target counts
anchor_plot <- anchor_df %>%
  dplyr::left_join(tf, by = c("Gene" = "target_symbol")) %>%
  dplyr::mutate(
    n_mature_miRNAs_targeting = ifelse(is.na(n_mature_miRNAs_targeting), 0,
                                       n_mature_miRNAs_targeting),
    Pathway = factor(Pathway,
                     levels = c("Antioxidant", "Ferroptosis", "miRNA biogenesis"))
  ) %>%
  dplyr::arrange(Pathway, n_mature_miRNAs_targeting) %>%
  dplyr::mutate(Gene = factor(Gene, levels = Gene))

cat("=== ANCHOR BAR CHART DATA ===\n")
print(anchor_plot, row.names = FALSE)
cat("\n")

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
pathway_colors <- c(
  "Antioxidant"      = "#FF7F00",
  "Ferroptosis"      = "#377EB8",
  "miRNA biogenesis" = "#984EA3"
)

p <- ggplot(anchor_plot,
            aes(x = Gene, y = n_mature_miRNAs_targeting, fill = Pathway)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = n_mature_miRNAs_targeting),
            hjust = -0.25, size = 3.3, fontface = "bold",
            color = "grey20") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pathway_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = c(0, 10, 20, 30, 40, 50, 60)) +
  labs(
    title    = "Anchor OA-DEG genes are validated cargo targets",
    subtitle = paste0("Number of distinct cargo mature miRNAs targeting each anchor gene\n",
                      "(miRTarBase validated, Non-Functional MTI excluded)"),
    x        = NULL,
    y        = "Number of cargo mature miRNAs (validated targeting)",
    fill     = "Anchor pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    axis.text.y      = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(size = 9),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    plot.margin      = margin(10, 20, 10, 10)
  )

ggsave(file.path(PROJECT_ROOT, "figures/fig_D2_anchor_bar.pdf"),
       p, width = 9, height = 6.5)

cat("Saved: figures/fig_D2_anchor_bar.pdf\n\n")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("=== SUMMARY ===\n")
summary_tbl <- anchor_plot %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(
    n_anchors_in_set       = dplyr::n(),
    n_targeted             = sum(n_mature_miRNAs_targeting > 0),
    total_miRNA_edges      = sum(n_mature_miRNAs_targeting),
    max_miRNAs_per_anchor  = max(n_mature_miRNAs_targeting),
    .groups = "drop"
  )
print(summary_tbl, row.names = FALSE)

zero_hit <- anchor_plot %>% dplyr::filter(n_mature_miRNAs_targeting == 0)
if (nrow(zero_hit) > 0) {
  cat("\nAnchors with ZERO cargo miRNAs targeting:\n")
  cat(paste(zero_hit$Gene, collapse = ", "), "\n")
}