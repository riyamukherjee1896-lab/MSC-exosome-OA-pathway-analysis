# =============================================================================
# 09_intersection_test.R
# -----------------------------------------------------------------------------
# Block D: Statistical test of cargo-target / OA-DEG convergence
#
# ANALYSES:
#   1. Hypergeometric test (multiple backgrounds for robustness)
#   2. Directional analysis (cargo targeting of OA-down vs OA-up genes)
#   3. Pathway-specific Fisher's exact tests (anchor gene sets)
#
# INPUT  : data/processed/cargo_miRNAs.rds
#          data/processed/oa_degs_filtered.rds
#          data/processed/mirna_targets_validated.rds
#          data/processed/de_results_annotated.csv (for cartilage background)
#
# OUTPUT : data/processed/intersection_results.rds
#          tables/table_D1_hypergeometric_tests.csv
#          tables/table_D2_directional_test.csv
#          tables/table_D3_pathway_enrichment.csv
#          tables/table_D4_intersection_genes.csv
#          logs/09_intersection_sessionInfo.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# Reload
oa_data <- readRDS(file.path(PROJECT_ROOT, "data/processed/oa_degs_filtered.rds"))
targets <- readRDS(file.path(PROJECT_ROOT, "data/processed/mirna_targets_validated.rds"))
de_full <- read.csv(file.path(PROJECT_ROOT, "data/processed/de_results_annotated.csv"),
                    stringsAsFactors = FALSE)

cat("============================================\n")
cat("BLOCK D REVISED: Tier-1 strict filter\n")
cat("============================================\n\n")

# Cargo-target universe at tier-1 strict (Functional MTI only)
cargo_edges_strict <- targets$cargo_edges %>%
  dplyr::filter(tier1_strict == TRUE)

cargo_target_strict <- unique(cargo_edges_strict$target_symbol)

oa_degs         <- oa_data$oa_degs$Gene
oa_up_genes     <- oa_data$oa_up$Gene
oa_down_genes   <- oa_data$oa_down$Gene
cartilage_expressed <- de_full$Gene

cat("=== INPUT SETS (TIER-1 STRICT) ===\n")
cat("Cargo-target universe (tier-1 strict):", length(cargo_target_strict), "\n")
cat("OA DEGs:", length(oa_degs), "\n")
cat("OA up:", length(oa_up_genes), "\n")
cat("OA down:", length(oa_down_genes), "\n")
cat("Cartilage expressed:", length(cartilage_expressed), "\n\n")

# Observed intersection
intersection <- intersect(cargo_target_strict, oa_degs)
intersection_up   <- intersect(cargo_target_strict, oa_up_genes)
intersection_down <- intersect(cargo_target_strict, oa_down_genes)

cat("=== OBSERVED INTERSECTION (TIER-1 STRICT) ===\n")
cat("Cargo-targets AND OA-DEG:  ", length(intersection), "\n")
cat("Cargo-targets AND OA-up:   ", length(intersection_up), "\n")
cat("Cargo-targets AND OA-down: ", length(intersection_down), "\n\n")

# =============================================================================
# Hypergeometric test with background = cartilage-expressed
# =============================================================================
bg <- unique(cartilage_expressed)

oa_in_bg    <- intersect(oa_degs, bg)
cargo_in_bg <- intersect(cargo_target_strict, bg)
obs         <- intersect(oa_in_bg, cargo_in_bg)

N <- length(bg)
K <- length(oa_in_bg)
n <- length(cargo_in_bg)
k <- length(obs)

expected <- (K * n) / N
p_val    <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

# Fisher for OR
tbl <- matrix(c(
  k,               # in cargo AND in OA
  n - k,           # in cargo NOT in OA
  K - k,           # NOT in cargo AND in OA
  N - n - K + k    # NOT in cargo NOT in OA
), nrow = 2)

fisher_result <- fisher.test(tbl, alternative = "greater")

cat("=== HYPERGEOMETRIC TEST ===\n")
cat("Background: cartilage-expressed (n =", N, ")\n")
cat(sprintf("  OA-DEGs in background:     %d\n", K))
cat(sprintf("  Cargo-targets in background: %d\n", n))
cat(sprintf("  Observed overlap:            %d\n", k))
cat(sprintf("  Expected under null:         %.1f\n", expected))
cat(sprintf("  Enrichment ratio:            %.2f\n", k / expected))
cat(sprintf("  Odds ratio (Fisher):         %.2f (95%% CI: %.2f - %.2f)\n",
            fisher_result$estimate,
            fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("  Hypergeometric p-value:      %.3g\n", p_val))
cat(sprintf("  -log10(p):                   %.1f\n\n", -log10(p_val)))

hypergeom_df <- data.frame(
  Background = "cartilage-expressed",
  N = N, K = K, n = n,
  Observed = k,
  Expected = round(expected, 1),
  Enrichment = round(k / expected, 3),
  OR = round(fisher_result$estimate, 3),
  OR_CI_lower = round(fisher_result$conf.int[1], 3),
  OR_CI_upper = round(fisher_result$conf.int[2], 3),
  P_value = signif(p_val, 3),
  Neg_log10_p = round(-log10(p_val), 2)
)
print(hypergeom_df, row.names = FALSE)
cat("\n")

# =============================================================================
# Directional: up vs down enrichment
# =============================================================================
up_in_bg    <- intersect(oa_up_genes, bg)
down_in_bg  <- intersect(oa_down_genes, bg)

# cargo targeting rate among up vs down
up_cargo    <- length(intersect(cargo_in_bg, up_in_bg))
up_total    <- length(up_in_bg)
down_cargo  <- length(intersect(cargo_in_bg, down_in_bg))
down_total  <- length(down_in_bg)

cat("=== DIRECTIONAL ANALYSIS ===\n")
cat(sprintf("OA-up in background:   %d, cargo-targeted: %d (%.1f%%)\n",
            up_total, up_cargo, 100 * up_cargo / up_total))
cat(sprintf("OA-down in background: %d, cargo-targeted: %d (%.1f%%)\n",
            down_total, down_cargo, 100 * down_cargo / down_total))

# 2x2: (in_cargo_target) x (OA_up vs OA_down)
tbl_dir <- matrix(c(
  up_cargo, up_total - up_cargo,
  down_cargo, down_total - down_cargo
), nrow = 2,
dimnames = list(c("Cargo_target", "Not_cargo_target"),
                c("OA_up", "OA_down")))
cat("\n2x2 contingency:\n")
print(tbl_dir)

fisher_dir <- fisher.test(tbl_dir)
cat(sprintf("\nFisher's exact test (up vs down):\n"))
cat(sprintf("  Odds ratio: %.3f (95%% CI: %.3f - %.3f)\n",
            fisher_dir$estimate,
            fisher_dir$conf.int[1], fisher_dir$conf.int[2]))
cat(sprintf("  p-value: %.3g\n", fisher_dir$p.value))
if (fisher_dir$estimate > 1) {
  cat("  Cargo preferentially targets OA-UP genes\n")
} else {
  cat("  Cargo preferentially targets OA-DOWN genes\n")
}
cat("\n")

# =============================================================================
# Pathway-specific enrichment (anchor sets, tier-1 strict)
# =============================================================================
anchor_sets <- oa_data$anchor_gene_sets

pathway_results <- data.frame()

for (set_name in names(anchor_sets)) {
  gene_set       <- anchor_sets[[set_name]]
  gene_set_in_bg <- intersect(gene_set, bg)
  if (length(gene_set_in_bg) == 0) next
  
  in_set_in_cargo   <- length(intersect(gene_set_in_bg, cargo_in_bg))
  in_set_not_cargo  <- length(gene_set_in_bg) - in_set_in_cargo
  not_set_in_cargo  <- length(cargo_in_bg) - in_set_in_cargo
  not_set_not_cargo <- length(bg) - length(gene_set_in_bg) -
    (length(cargo_in_bg) - in_set_in_cargo)
  
  tbl_path <- matrix(c(
    in_set_in_cargo,  in_set_not_cargo,
    not_set_in_cargo, not_set_not_cargo
  ), nrow = 2)
  
  fisher_path <- fisher.test(tbl_path, alternative = "greater")
  
  pathway_results <- rbind(pathway_results, data.frame(
    Pathway            = set_name,
    Pathway_size_in_bg = length(gene_set_in_bg),
    Cargo_targeted     = in_set_in_cargo,
    Pct_targeted       = round(100 * in_set_in_cargo / length(gene_set_in_bg), 1),
    Odds_ratio         = round(fisher_path$estimate, 2),
    OR_CI_lower        = round(fisher_path$conf.int[1], 2),
    OR_CI_upper        = round(fisher_path$conf.int[2], 2),
    P_value            = signif(fisher_path$p.value, 3),
    stringsAsFactors   = FALSE
  ))
}

cat("=== PATHWAY-SPECIFIC ENRICHMENT (tier-1 strict) ===\n")
print(pathway_results, row.names = FALSE)
cat("\n")

# Save
intersection_genes <- data.frame(Gene = intersection, stringsAsFactors = FALSE) %>%
  dplyr::left_join(
    oa_data$oa_degs[, c("Gene", "log2FoldChange", "padj", "Direction")],
    by = "Gene"
  ) %>%
  dplyr::arrange(padj)

write.csv(intersection_genes,
          file.path(PROJECT_ROOT, "tables/table_D4_intersection_genes.csv"),
          row.names = FALSE)
write.csv(hypergeom_df,
          file.path(PROJECT_ROOT, "tables/table_D1_hypergeometric_tests.csv"),
          row.names = FALSE)
write.csv(pathway_results,
          file.path(PROJECT_ROOT, "tables/table_D3_pathway_enrichment.csv"),
          row.names = FALSE)

saveRDS(list(
  cargo_target_strict = cargo_target_strict,
  intersection = intersection,
  intersection_genes = intersection_genes,
  hypergeom = hypergeom_df,
  directional = list(
    up_cargo = up_cargo, up_total = up_total,
    down_cargo = down_cargo, down_total = down_total,
    fisher_result = fisher_dir
  ),
  pathway_results = pathway_results,
  evidence_filter = "tier-1 strict (Functional MTI + Functional MTI Weak)"
), file.path(PROJECT_ROOT, "data/processed/intersection_results.rds"))
cat("Saved: intersection_results.rds\n\n")

cat("============================================\n")
cat("BLOCK D REVISED COMPLETE\n")
cat("============================================\n")