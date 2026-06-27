# =============================================================================
# 18_synovium_GO_KEGG_enrichment.R
# -----------------------------------------------------------------------------
# Block G: GO/KEGG enrichment for synovium cargo-target ∩ OA-DEG intersection
#
# Mirrors Block E (cartilage enrichment) for direct cross-tissue comparison.
#
# INPUT  : data/processed/syn_intersection_results.rds
#          data/processed/syn_de_results.rds
#          data/processed/mirna_targets_validated.rds
# OUTPUT : data/processed/syn_enrichment_results.rds
#          tables/table_G1_synovium_GO_BP_enrichment.csv
#          tables/table_G2_synovium_KEGG_enrichment.csv
# =============================================================================

# Packages should already be installed from Block E
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("BLOCK G: SYNOVIUM GO/KEGG ENRICHMENT\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
syn_inter <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/syn_intersection_results.rds"))
syn_de <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/syn_de_results.rds"))
targets <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/mirna_targets_validated.rds"))

# Gene sets
intersection_genes <- syn_inter$intersection            # 1394 genes
cargo_targetable    <- unique(targets$cargo_edges$target_symbol)
synovium_expressed  <- syn_de$Gene                       # 13039 genes

# Background: synovium-expressed AND cargo-targetable (matches Block E logic)
background <- intersect(synovium_expressed, cargo_targetable)

cat("=== GENE SETS ===\n")
cat(sprintf("Intersection (query):       %d\n", length(intersection_genes)))
cat(sprintf("Synovium-expressed:         %d\n", length(synovium_expressed)))
cat(sprintf("Cargo-targetable (any):     %d\n", length(cargo_targetable)))
cat(sprintf("Background (intersection):  %d\n\n", length(background)))

# -----------------------------------------------------------------------------
# Convert SYMBOL -> ENTREZ
# -----------------------------------------------------------------------------
cat("Converting gene symbols to ENTREZ IDs...\n")

query_entrez <- bitr(intersection_genes,
                     fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
bg_entrez    <- bitr(background,
                     fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

cat(sprintf("  Query: %d / %d symbols mapped (%.1f%%)\n",
            nrow(query_entrez), length(intersection_genes),
            100 * nrow(query_entrez) / length(intersection_genes)))
cat(sprintf("  Background: %d / %d symbols mapped (%.1f%%)\n\n",
            nrow(bg_entrez), length(background),
            100 * nrow(bg_entrez) / length(background)))

# -----------------------------------------------------------------------------
# GO BP enrichment
# -----------------------------------------------------------------------------
cat("=== GO BIOLOGICAL PROCESS ENRICHMENT ===\n")
t0 <- Sys.time()
go_bp <- enrichGO(
  gene          = query_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
cat(sprintf("Completed in %.1f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  bp_signif <- go_bp@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust)
  
  cat(sprintf("Total GO BP terms: %d\n", nrow(go_bp@result)))
  cat(sprintf("Significant at FDR < 0.05: %d\n\n", nrow(bp_signif)))
  
  cat("Top 20 enriched GO BP terms:\n")
  print(head(bp_signif[, c("ID", "Description", "GeneRatio",
                           "BgRatio", "pvalue", "p.adjust", "Count")], 20),
        row.names = FALSE)
  
  write.csv(go_bp@result,
            file.path(PROJECT_ROOT, "tables/table_G1_synovium_GO_BP_enrichment.csv"),
            row.names = FALSE)
  cat("\nSaved: tables/table_G1_synovium_GO_BP_enrichment.csv\n\n")
} else {
  cat("No GO BP enrichment at FDR < 0.05\n\n")
}

# -----------------------------------------------------------------------------
# KEGG pathway enrichment
# -----------------------------------------------------------------------------
cat("=== KEGG PATHWAY ENRICHMENT ===\n")
t0 <- Sys.time()
kegg <- enrichKEGG(
  gene          = query_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
cat(sprintf("Completed in %.1f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  kegg_signif <- kegg@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust)
  
  cat(sprintf("Total KEGG pathways: %d\n", nrow(kegg@result)))
  cat(sprintf("Significant at FDR < 0.05: %d\n\n", nrow(kegg_signif)))
  
  cat("Top 20 enriched KEGG pathways:\n")
  print(head(kegg_signif[, c("ID", "Description", "GeneRatio",
                             "BgRatio", "pvalue", "p.adjust", "Count")], 20),
        row.names = FALSE)
  
  write.csv(kegg@result,
            file.path(PROJECT_ROOT, "tables/table_G2_synovium_KEGG_enrichment.csv"),
            row.names = FALSE)
  cat("\nSaved: tables/table_G2_synovium_KEGG_enrichment.csv\n\n")
} else {
  cat("No KEGG enrichment at FDR < 0.05\n\n")
}

# -----------------------------------------------------------------------------
# Cross-tissue comparison: which terms appear in BOTH?
# -----------------------------------------------------------------------------
cat("=== CROSS-TISSUE PATHWAY COMPARISON ===\n")

cart_enr <- readRDS(file.path(PROJECT_ROOT, "data/processed/enrichment_results.rds"))
cart_bp <- cart_enr$go_bp@result %>% dplyr::filter(p.adjust < 0.05)
cart_kegg <- cart_enr$kegg@result %>% dplyr::filter(p.adjust < 0.05)

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  syn_bp <- go_bp@result %>% dplyr::filter(p.adjust < 0.05)
  shared_bp <- intersect(cart_bp$ID, syn_bp$ID)
  cart_only_bp <- setdiff(cart_bp$ID, syn_bp$ID)
  syn_only_bp <- setdiff(syn_bp$ID, cart_bp$ID)
  
  cat(sprintf("GO BP terms — Cartilage signif: %d, Synovium signif: %d\n",
              nrow(cart_bp), nrow(syn_bp)))
  cat(sprintf("  Shared between tissues: %d\n", length(shared_bp)))
  cat(sprintf("  Cartilage-only:         %d\n", length(cart_only_bp)))
  cat(sprintf("  Synovium-only:          %d\n", length(syn_only_bp)))
  
  if (length(shared_bp) > 0) {
    shared_terms <- cart_bp %>%
      dplyr::filter(ID %in% shared_bp) %>%
      dplyr::arrange(p.adjust) %>%
      head(15)
    cat("\nTop 15 GO BP terms enriched in BOTH tissues (sorted by cartilage FDR):\n")
    print(shared_terms[, c("ID", "Description")], row.names = FALSE)
  }
  cat("\n")
}

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  syn_kegg <- kegg@result %>% dplyr::filter(p.adjust < 0.05)
  shared_kegg <- intersect(cart_kegg$ID, syn_kegg$ID)
  cart_only_kegg <- setdiff(cart_kegg$ID, syn_kegg$ID)
  syn_only_kegg <- setdiff(syn_kegg$ID, cart_kegg$ID)
  
  cat(sprintf("KEGG pathways — Cartilage signif: %d, Synovium signif: %d\n",
              nrow(cart_kegg), nrow(syn_kegg)))
  cat(sprintf("  Shared between tissues: %d\n", length(shared_kegg)))
  cat(sprintf("  Cartilage-only:         %d\n", length(cart_only_kegg)))
  cat(sprintf("  Synovium-only:          %d\n", length(syn_only_kegg)))
  
  if (length(shared_kegg) > 0) {
    shared_kegg_terms <- cart_kegg %>%
      dplyr::filter(ID %in% shared_kegg) %>%
      dplyr::arrange(p.adjust)
    cat("\nKEGG pathways enriched in BOTH tissues:\n")
    print(shared_kegg_terms[, c("ID", "Description")], row.names = FALSE)
  }
  cat("\n")
}

# -----------------------------------------------------------------------------
# Save consolidated results
# -----------------------------------------------------------------------------
syn_enrichment <- list(
  go_bp      = go_bp,
  kegg       = kegg,
  query_size = length(intersection_genes),
  bg_size    = length(background)
)
saveRDS(syn_enrichment,
        file.path(PROJECT_ROOT, "data/processed/syn_enrichment_results.rds"))
cat("Saved: data/processed/syn_enrichment_results.rds\n\n")

cat("============================================\n")
cat("BLOCK G COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")