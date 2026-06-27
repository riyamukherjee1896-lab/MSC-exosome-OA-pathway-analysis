# =============================================================================
# 10_GO_KEGG_enrichment.R
# -----------------------------------------------------------------------------
# Block E: GO biological process + KEGG pathway enrichment on
#          cargo-target / OA-DEG intersection
#
# BACKGROUND: custom, cartilage-expressed AND cargo-targetable
# GENE SET:   2,918 intersection genes from Block D
# FDR:        Benjamini-Hochberg, threshold 0.05
# =============================================================================

# Install if needed
required <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("BLOCK E: GO/KEGG Enrichment\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
intersection_data <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/intersection_results.rds"))
oa_data <- readRDS(
  file.path(PROJECT_ROOT, "data/processed/oa_degs_filtered.rds"))
de_full <- read.csv(
  file.path(PROJECT_ROOT, "data/processed/de_results_annotated.csv"),
  stringsAsFactors = FALSE)

intersection_genes <- intersection_data$intersection
cargo_target_strict <- intersection_data$cargo_target_strict

# Background: cartilage-expressed AND cargo-targetable
background <- intersect(de_full$Gene, cargo_target_strict)

cat("=== GENE SETS ===\n")
cat("Intersection (query):", length(intersection_genes), "\n")
cat("Background:", length(background), "\n")
cat("(Background = cartilage-expressed AND cargo-targetable)\n\n")

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

cat(sprintf("  Query: %d / %d symbols mapped\n",
            nrow(query_entrez), length(intersection_genes)))
cat(sprintf("  Background: %d / %d symbols mapped\n\n",
            nrow(bg_entrez), length(background)))

# -----------------------------------------------------------------------------
# GO BP enrichment
# -----------------------------------------------------------------------------
cat("=== GO BIOLOGICAL PROCESS ENRICHMENT ===\n")
cat("Running enrichGO (BP)...\n")

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
cat(sprintf("  Completed in %.1f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  # Top 20 enriched BP terms
  bp_top <- go_bp@result %>%
    dplyr::arrange(p.adjust) %>%
    head(20) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count)
  
  cat("\nTop 20 enriched GO BP terms:\n")
  print(bp_top, row.names = FALSE)
  
  # Save full table
  write.csv(go_bp@result,
            file.path(PROJECT_ROOT, "tables/table_E1_GO_BP_enrichment.csv"),
            row.names = FALSE)
  cat("\nSaved: tables/table_E1_GO_BP_enrichment.csv (",
      nrow(go_bp@result), "rows total)\n")
} else {
  cat("No GO BP enrichment at FDR < 0.05\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# KEGG pathway enrichment
# -----------------------------------------------------------------------------
cat("=== KEGG PATHWAY ENRICHMENT ===\n")
cat("Running enrichKEGG...\n")

t0 <- Sys.time()
kegg <- enrichKEGG(
  gene          = query_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
cat(sprintf("  Completed in %.1f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  kegg_top <- kegg@result %>%
    dplyr::arrange(p.adjust) %>%
    head(20) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count)
  
  cat("\nTop 20 enriched KEGG pathways:\n")
  print(kegg_top, row.names = FALSE)
  
  write.csv(kegg@result,
            file.path(PROJECT_ROOT, "tables/table_E2_KEGG_enrichment.csv"),
            row.names = FALSE)
  cat("\nSaved: tables/table_E2_KEGG_enrichment.csv (",
      nrow(kegg@result), "rows total)\n")
} else {
  cat("No KEGG enrichment at FDR < 0.05\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# Targeted check: do our anchor pathway terms appear?
# -----------------------------------------------------------------------------
cat("=== TARGETED ANCHOR-RELEVANT TERMS ===\n")

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  anchor_keywords <- c("oxidative", "redox", "reactive oxygen", "glutathione",
                       "ferroptosis", "iron", "lipid peroxidation",
                       "miRNA", "microRNA", "RNA interference", "gene silencing",
                       "NADPH", "cellular response to oxidative")
  
  hits_in_bp <- go_bp@result %>%
    dplyr::filter(grepl(paste(anchor_keywords, collapse = "|"),
                        Description, ignore.case = TRUE)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(ID, Description, GeneRatio, p.adjust, Count)
  
  if (nrow(hits_in_bp) > 0) {
    cat("Anchor-relevant GO BP terms:\n")
    print(hits_in_bp, row.names = FALSE)
  } else {
    cat("No anchor-relevant GO BP terms in significant enrichment\n")
  }
}

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  kegg_keywords <- c("Ferroptosis", "Glutathione", "NRF2", "oxidative",
                     "RNA", "Ribosome")
  
  hits_in_kegg <- kegg@result %>%
    dplyr::filter(grepl(paste(kegg_keywords, collapse = "|"),
                        Description, ignore.case = TRUE)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(ID, Description, GeneRatio, p.adjust, Count)
  
  if (nrow(hits_in_kegg) > 0) {
    cat("\nAnchor-relevant KEGG pathways:\n")
    print(hits_in_kegg, row.names = FALSE)
  } else {
    cat("\nNo anchor-relevant KEGG pathways in significant enrichment\n")
  }
}
cat("\n")

# -----------------------------------------------------------------------------
# Save everything
# -----------------------------------------------------------------------------
enrichment_results <- list(
  go_bp      = go_bp,
  kegg       = kegg,
  query_size = length(intersection_genes),
  bg_size    = length(background)
)
saveRDS(enrichment_results,
        file.path(PROJECT_ROOT, "data/processed/enrichment_results.rds"))
cat("Saved: data/processed/enrichment_results.rds\n\n")

cat("============================================\n")
cat("BLOCK E COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")