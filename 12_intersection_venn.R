# =============================================================================
# 12_intersection_venn.R
# -----------------------------------------------------------------------------
# Figure D1: Intersection Venn diagram (cargo targets vs OA DEGs)
#
# INPUT  : data/processed/intersection_results.rds
#          data/processed/oa_degs_filtered.rds
#          data/processed/de_results_annotated.csv
# OUTPUT : figures/fig_D1_intersection_venn.pdf
# =============================================================================

# Install eulerr if needed
if (!requireNamespace("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
}
suppressPackageStartupMessages({
  library(eulerr)
  library(grid)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
inter <- readRDS(file.path(PROJECT_ROOT, "data/processed/intersection_results.rds"))
oa    <- readRDS(file.path(PROJECT_ROOT, "data/processed/oa_degs_filtered.rds"))
de_full <- read.csv(file.path(PROJECT_ROOT, "data/processed/de_results_annotated.csv"),
                    stringsAsFactors = FALSE)

cargo_targets <- inter$cargo_target_strict
oa_degs       <- oa$oa_degs$Gene
cartilage     <- de_full$Gene

cat("=== SET SIZES ===\n")
cat("Cartilage-expressed:", length(cartilage), "\n")
cat("Cargo-target universe (tier-1 strict):", length(cargo_targets), "\n")
cat("OA DEGs:", length(oa_degs), "\n\n")

# -----------------------------------------------------------------------------
# Compute Euler set sizes
# -----------------------------------------------------------------------------
# Restrict cargo and OA to cartilage-expressed for fair comparison
cargo_in_cart <- intersect(cargo_targets, cartilage)
oa_in_cart    <- intersect(oa_degs, cartilage)

# Three sets for Euler: Cartilage (all), Cargo (within cartilage), OA (within cartilage)
# We want region counts:
#   Cartilage only (not cargo, not OA):       A
#   Cargo only (within cart, not OA):          B
#   OA only (within cart, not cargo):          C
#   Cargo AND OA (within cart):                D
#
#   Total cartilage = A + B + C + D = 15,440 (approximately)

overlap_cargo_oa  <- length(intersect(cargo_in_cart, oa_in_cart))
cargo_only        <- length(setdiff(cargo_in_cart, oa_in_cart))
oa_only           <- length(setdiff(oa_in_cart, cargo_in_cart))
cart_only         <- length(cartilage) - cargo_only - oa_only - overlap_cargo_oa

cat("=== REGION COUNTS (in cartilage-expressed universe) ===\n")
cat("Cartilage only (not cargo, not OA):", cart_only, "\n")
cat("Cargo-target only (not OA):         ", cargo_only, "\n")
cat("OA DEG only (not cargo-target):     ", oa_only, "\n")
cat("Cargo-target AND OA DEG:            ", overlap_cargo_oa, "\n")
cat("Sum check:",
    cart_only + cargo_only + oa_only + overlap_cargo_oa,
    "(should be ~", length(cartilage), ")\n\n")

# -----------------------------------------------------------------------------
# Build Euler diagram
# -----------------------------------------------------------------------------
fit <- euler(c(
  "Cartilage-expressed"                                       = cart_only,
  "Cargo-target"                                               = 0,
  "OA DEG"                                                     = 0,
  "Cartilage-expressed&Cargo-target"                           = cargo_only,
  "Cartilage-expressed&OA DEG"                                 = oa_only,
  "Cargo-target&OA DEG"                                        = 0,
  "Cartilage-expressed&Cargo-target&OA DEG"                    = overlap_cargo_oa
))

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
pdf(file.path(PROJECT_ROOT, "figures/fig_D1_intersection_venn.pdf"),
    width = 8, height = 7)

plot(fit,
     fills      = list(fill = c("#C0C0C0", "#E41A1C", "#377EB8"),
                       alpha = c(0.20, 0.45, 0.45)),
     edges      = list(col = c("grey50", "#B22222", "#1F3A5F"), lwd = 2),
     labels     = list(font = 2, cex = 1.0,
                       col = c("grey30", "#8B0000", "#0D2642")),
     quantities = list(cex = 1.1, font = 2),
     main       = list(label = "Cargo-target and OA-DEG overlap within cartilage-expressed transcriptome",
                       cex = 1.0, font = 2))

# Add informational caption below
grid.text(paste0("Total cartilage-expressed: ", length(cartilage),
                 "  |  Cargo-target coverage: ",
                 round(100 * length(cargo_in_cart) / length(cartilage), 1),
                 "%  |  Cargo ∩ OA DEG: ", overlap_cargo_oa, " genes"),
          x = 0.5, y = 0.03, gp = gpar(fontsize = 9, col = "grey25"))

dev.off()

cat("Saved: figures/fig_D1_intersection_venn.pdf\n\n")

# -----------------------------------------------------------------------------
# Report summary
# -----------------------------------------------------------------------------
cat("=== FIGURE SUMMARY ===\n")
cat(sprintf("%-45s %6d\n", "Cartilage-expressed (total):",     length(cartilage)))
cat(sprintf("%-45s %6d  (%.1f%%)\n",
            "  Cargo-target coverage:", length(cargo_in_cart),
            100 * length(cargo_in_cart) / length(cartilage)))
cat(sprintf("%-45s %6d  (%.1f%%)\n",
            "  OA-DEG coverage:", length(oa_in_cart),
            100 * length(oa_in_cart) / length(cartilage)))
cat(sprintf("%-45s %6d\n", "  Cargo-target AND OA-DEG:", overlap_cargo_oa))
cat(sprintf("%-45s %.2f\n", "Observed/expected overlap ratio:",
            overlap_cargo_oa /
              ((length(cargo_in_cart) * length(oa_in_cart)) / length(cartilage))))