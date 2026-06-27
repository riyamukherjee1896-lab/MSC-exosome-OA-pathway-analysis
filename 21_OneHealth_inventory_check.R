# =============================================================================
# 21_OneHealth_inventory_check.R
# -----------------------------------------------------------------------------
# Stage 1: Verify what mammalian species miRNA data is accessible.
#
# Tests miRBaseConverter and miRBase API availability before committing to
# full cross-species analysis.
# =============================================================================

# Install miRBaseConverter if needed
if (!requireNamespace("miRBaseConverter", quietly = TRUE)) {
  cat("Installing miRBaseConverter from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("miRBaseConverter", update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(miRBaseConverter)
})

cat("============================================\n")
cat("ONE HEALTH INVENTORY CHECK\n")
cat("============================================\n\n")

# -----------------------------------------------------------------------------
# 1. Check miRBaseConverter version and supported species
# -----------------------------------------------------------------------------
cat("=== miRBaseConverter version ===\n")
cat("Version:", as.character(packageVersion("miRBaseConverter")), "\n\n")

# -----------------------------------------------------------------------------
# 2. Test what species are queryable
# -----------------------------------------------------------------------------
cat("=== TESTING SPECIES AVAILABILITY ===\n")

# Get miRBase v22 dataset
data("miRNATable")
cat("Total entries in miRNATable:", nrow(miRNATable), "\n")
cat("Columns:", paste(colnames(miRNATable), collapse = ", "), "\n\n")

# Show species breakdown
cat("Species in miRNATable (top 20 by count):\n")
species_counts <- sort(table(miRNATable$Species), decreasing = TRUE)
print(head(species_counts, 20))
cat("\n")

# Specifically check our target species
target_species <- c("hsa", "mmu", "rno", "cfa", "eca", "bta")
cat("Target species check:\n")
for (sp in target_species) {
  in_data <- sp %in% miRNATable$Species
  count <- if (in_data) sum(miRNATable$Species == sp) else 0
  species_full <- switch(sp,
                         "hsa" = "Human (Homo sapiens)",
                         "mmu" = "Mouse (Mus musculus)",
                         "rno" = "Rat (Rattus norvegicus)",
                         "cfa" = "Dog (Canis familiaris)",
                         "eca" = "Horse (Equus caballus)",
                         "bta" = "Cattle (Bos taurus)",
                         sp)
  cat(sprintf("  %s (%-5s)  : %s  (%d miRNAs)\n",
              ifelse(in_data, "OK    ", "MISSING"),
              sp, species_full, count))
}
cat("\n")

# -----------------------------------------------------------------------------
# 3. Inspect what we can actually look up
# -----------------------------------------------------------------------------
cat("=== SAMPLE LOOKUP: hsa-miR-21-5p ===\n")

# Find human miR-21-5p
test_query <- miRNATable[miRNATable$Species == "hsa" &
                           grepl("miR-21$|miR-21-5p$", miRNATable$Mature1_ID,
                                 ignore.case = TRUE), ]
cat("Human miR-21 entries:\n")
print(test_query, row.names = FALSE)
cat("\n")

# Same query for mouse
mouse_mir21 <- miRNATable[miRNATable$Species == "mmu" &
                            grepl("miR-21$|miR-21-5p$", miRNATable$Mature1_ID,
                                  ignore.case = TRUE), ]
cat("Mouse miR-21 entries:\n")
print(mouse_mir21, row.names = FALSE)
cat("\n")

cat("============================================\n")
cat("INVENTORY COMPLETE\n")
cat("============================================\n")