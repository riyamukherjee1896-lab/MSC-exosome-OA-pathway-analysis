# =============================================================================
# 19_synovium_enrichment_figure.R
# -----------------------------------------------------------------------------
# Figure G1: Synovium GO/KEGG enrichment dotplot (2-panel)
# Parallels fig_E1 for direct cross-tissue visual comparison
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

enr <- readRDS(file.path(PROJECT_ROOT, "data/processed/syn_enrichment_results.rds"))
go_bp_df <- enr$go_bp@result
kegg_df  <- enr$kegg@result

cat("=== INPUT ===\n")
cat("GO BP at FDR < 0.05:", sum(go_bp_df$p.adjust < 0.05), "\n")
cat("KEGG at FDR < 0.05:", sum(kegg_df$p.adjust < 0.05), "\n\n")

parse_ratio <- function(ratio_str) {
  parts <- strsplit(ratio_str, "/")
  sapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]))
}

# Panel A: Top 15 GO BP
go_top <- go_bp_df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  head(15) %>%
  dplyr::mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    NegLog10FDR   = -log10(p.adjust),
    Description   = ifelse(nchar(Description) > 55,
                           paste0(substr(Description, 1, 52), "..."),
                           Description),
    Description   = factor(Description, levels = rev(Description))
  )

p_go <- ggplot(go_top, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = NegLog10FDR)) +
  scale_color_gradient(low = "#9DC8E6", high = "#0D2642",
                       name = expression(-log[10]*"(FDR)")) +
  scale_size_continuous(name = "Gene count", range = c(3, 8)) +
  labs(title    = "A  GO Biological Process",
       subtitle = paste0("Top 15 enriched terms (synovium intersection, n = ",
                         enr$query_size, ")"),
       x = "Gene ratio", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    axis.text.y        = element_text(size = 8.5),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "right",
    legend.box         = "vertical"
  )

# Panel B: Top 15 KEGG
kegg_top <- kegg_df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  head(15) %>%
  dplyr::mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    NegLog10FDR   = -log10(p.adjust),
    Description   = ifelse(nchar(Description) > 55,
                           paste0(substr(Description, 1, 52), "..."),
                           Description),
    Description   = factor(Description, levels = rev(Description))
  )

p_kegg <- ggplot(kegg_top, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = NegLog10FDR)) +
  scale_color_gradient(low = "#F5B4B4", high = "#8B0000",
                       name = expression(-log[10]*"(FDR)")) +
  scale_size_continuous(name = "Gene count", range = c(3, 8)) +
  labs(title    = "B  KEGG Pathways",
       subtitle = paste0("Top 15 enriched pathways (",
                         sum(kegg_df$p.adjust < 0.05), " significant total)"),
       x = "Gene ratio", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    axis.text.y        = element_text(size = 8.5),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position    = "right",
    legend.box         = "vertical"
  )

combined <- p_go / p_kegg +
  plot_annotation(
    title    = "Pathway enrichment of synovium cargo-target / OA-DEG intersection",
    subtitle = paste0("1,394 intersection genes; background: ",
                      enr$bg_size, " synovium-expressed and miRNA-targetable genes"),
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, color = "grey30"))
  )

ggsave(file.path(PROJECT_ROOT, "figures/fig_G1_synovium_GO_KEGG_dotplot.pdf"),
       combined, width = 11, height = 10)

cat("Saved: figures/fig_G1_synovium_GO_KEGG_dotplot.pdf\n\n")

cat("=== GO BP TERMS IN FIGURE ===\n")
print(go_top[, c("Description", "Count", "p.adjust")], row.names = FALSE)
cat("\n=== KEGG PATHWAYS IN FIGURE ===\n")
print(kegg_top[, c("Description", "Count", "p.adjust")], row.names = FALSE)