# =============================================================================
# 20_cross_tissue_pathway_figure.R
# -----------------------------------------------------------------------------
# Figure G2: Cross-tissue shared pathway heatmap (cartilage vs synovium)
#
# Panel A: Top 20 GO BP terms enriched in BOTH tissues
# Panel B: Shared KEGG pathways (3 total)
#
# INPUT  : data/processed/enrichment_results.rds (cartilage)
#          data/processed/syn_enrichment_results.rds (synovium)
# OUTPUT : figures/fig_G2_cross_tissue_pathways.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
cart_enr <- readRDS(file.path(PROJECT_ROOT, "data/processed/enrichment_results.rds"))
syn_enr  <- readRDS(file.path(PROJECT_ROOT, "data/processed/syn_enrichment_results.rds"))

cart_bp <- cart_enr$go_bp@result %>% dplyr::filter(p.adjust < 0.05)
syn_bp  <- syn_enr$go_bp@result %>% dplyr::filter(p.adjust < 0.05)

cart_kegg <- cart_enr$kegg@result %>% dplyr::filter(p.adjust < 0.05)
syn_kegg  <- syn_enr$kegg@result %>% dplyr::filter(p.adjust < 0.05)

# -----------------------------------------------------------------------------
# Find shared terms
# -----------------------------------------------------------------------------
shared_bp_ids <- intersect(cart_bp$ID, syn_bp$ID)
shared_kegg_ids <- intersect(cart_kegg$ID, syn_kegg$ID)

cat("Shared GO BP terms:", length(shared_bp_ids), "\n")
cat("Shared KEGG pathways:", length(shared_kegg_ids), "\n\n")

# -----------------------------------------------------------------------------
# Build top-20 shared GO BP heatmap data
# -----------------------------------------------------------------------------
shared_bp <- data.frame(
  ID = shared_bp_ids,
  Description = cart_bp$Description[match(shared_bp_ids, cart_bp$ID)],
  Cartilage_FDR = cart_bp$p.adjust[match(shared_bp_ids, cart_bp$ID)],
  Synovium_FDR  = syn_bp$p.adjust[match(shared_bp_ids, syn_bp$ID)],
  stringsAsFactors = FALSE
)

# Combined ranking by geometric mean of -log10(FDR)
shared_bp$combined_score <- (-log10(shared_bp$Cartilage_FDR) *
                               -log10(shared_bp$Synovium_FDR))^0.5

shared_bp_top <- shared_bp %>%
  dplyr::arrange(dplyr::desc(combined_score)) %>%
  head(20)

# Long format for ggplot
shared_bp_long <- shared_bp_top %>%
  dplyr::mutate(
    Cartilage = -log10(Cartilage_FDR),
    Synovium  = -log10(Synovium_FDR),
    Description = ifelse(nchar(Description) > 50,
                         paste0(substr(Description, 1, 47), "..."),
                         Description)
  ) %>%
  dplyr::select(Description, Cartilage, Synovium, combined_score) %>%
  tidyr::pivot_longer(cols = c(Cartilage, Synovium),
                      names_to = "Tissue", values_to = "NegLog10FDR")

# Preserve order: most significant at top
shared_bp_long$Description <- factor(shared_bp_long$Description,
                                     levels = rev(unique(shared_bp_top %>%
                                                           dplyr::arrange(dplyr::desc(combined_score)) %>%
                                                           dplyr::mutate(Description = ifelse(
                                                             nchar(Description) > 50,
                                                             paste0(substr(Description, 1, 47), "..."),
                                                             Description)) %>%
                                                           dplyr::pull(Description))))
shared_bp_long$Tissue <- factor(shared_bp_long$Tissue,
                                levels = c("Cartilage", "Synovium"))

cat("=== TOP 20 SHARED GO BP TERMS (by combined -log10 FDR) ===\n")
print(shared_bp_top[, c("Description", "Cartilage_FDR", "Synovium_FDR")],
      row.names = FALSE, digits = 3)
cat("\n")

# -----------------------------------------------------------------------------
# Build shared KEGG heatmap data
# -----------------------------------------------------------------------------
shared_kegg <- data.frame(
  ID = shared_kegg_ids,
  Description = cart_kegg$Description[match(shared_kegg_ids, cart_kegg$ID)],
  Cartilage_FDR = cart_kegg$p.adjust[match(shared_kegg_ids, cart_kegg$ID)],
  Synovium_FDR  = syn_kegg$p.adjust[match(shared_kegg_ids, syn_kegg$ID)],
  stringsAsFactors = FALSE
)
shared_kegg$combined_score <- (-log10(shared_kegg$Cartilage_FDR) *
                                 -log10(shared_kegg$Synovium_FDR))^0.5
shared_kegg <- shared_kegg %>% dplyr::arrange(dplyr::desc(combined_score))

shared_kegg_long <- shared_kegg %>%
  dplyr::mutate(
    Cartilage = -log10(Cartilage_FDR),
    Synovium  = -log10(Synovium_FDR)
  ) %>%
  dplyr::select(Description, Cartilage, Synovium, combined_score) %>%
  tidyr::pivot_longer(cols = c(Cartilage, Synovium),
                      names_to = "Tissue", values_to = "NegLog10FDR")
shared_kegg_long$Description <- factor(shared_kegg_long$Description,
                                       levels = rev(shared_kegg$Description))
shared_kegg_long$Tissue <- factor(shared_kegg_long$Tissue,
                                  levels = c("Cartilage", "Synovium"))

cat("=== SHARED KEGG PATHWAYS ===\n")
print(shared_kegg[, c("Description", "Cartilage_FDR", "Synovium_FDR")],
      row.names = FALSE, digits = 3)
cat("\n")

# -----------------------------------------------------------------------------
# Color scale: shared across both panels for direct comparison
# -----------------------------------------------------------------------------
all_neglog <- c(shared_bp_long$NegLog10FDR, shared_kegg_long$NegLog10FDR)
color_min <- 1.3   # corresponds to FDR = 0.05 (just barely significant)
color_max <- max(all_neglog) * 1.05

# -----------------------------------------------------------------------------
# Panel A — GO BP
# -----------------------------------------------------------------------------
p_go <- ggplot(shared_bp_long,
               aes(x = Tissue, y = Description, fill = NegLog10FDR)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", NegLog10FDR)),
            color = ifelse(shared_bp_long$NegLog10FDR > 6, "white", "grey20"),
            size = 3, fontface = "bold") +
  scale_fill_gradient(
    low = "#FFF5E0", high = "#7B1F1F",
    limits = c(color_min, color_max),
    name = expression(-log[10]*"(FDR)")
  ) +
  scale_x_discrete(position = "top") +
  labs(title    = "A  Top 20 GO BP terms enriched in BOTH cartilage and synovium",
       subtitle = "155 GO biological process terms reach FDR < 0.05 in both tissues",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    axis.text.x        = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 8.5),
    panel.grid         = element_blank(),
    legend.position    = "right"
  )

# -----------------------------------------------------------------------------
# Panel B — KEGG
# -----------------------------------------------------------------------------
p_kegg <- ggplot(shared_kegg_long,
                 aes(x = Tissue, y = Description, fill = NegLog10FDR)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", NegLog10FDR)),
            color = ifelse(shared_kegg_long$NegLog10FDR > 6, "white", "grey20"),
            size = 3.2, fontface = "bold") +
  scale_fill_gradient(
    low = "#FFF5E0", high = "#7B1F1F",
    limits = c(color_min, color_max),
    name = expression(-log[10]*"(FDR)")
  ) +
  scale_x_discrete(position = "top") +
  labs(title    = "B  KEGG pathways enriched in BOTH tissues",
       subtitle = "3 KEGG pathways reach FDR < 0.05 in both cartilage and synovium",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    axis.text.x        = element_text(face = "bold", size = 11),
    axis.text.y        = element_text(size = 9),
    panel.grid         = element_blank(),
    legend.position    = "right"
  )

# -----------------------------------------------------------------------------
# Combine
# -----------------------------------------------------------------------------
combined <- p_go / p_kegg +
  plot_layout(heights = c(20, 4)) +
  plot_annotation(
    title    = "Cross-tissue convergent pathway enrichment",
    subtitle = "Pathways enriched (FDR < 0.05) in cargo-target / OA-DEG intersection of BOTH cartilage (GSE114007) and synovium (GSE55235)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, color = "grey30"))
  )

ggsave(file.path(PROJECT_ROOT, "figures/fig_G2_cross_tissue_pathways.pdf"),
       combined, width = 11, height = 11)

cat("Saved: figures/fig_G2_cross_tissue_pathways.pdf\n")