# =============================================================================
# 15_synovium_download_and_inspect.R
# -----------------------------------------------------------------------------
# Stage 1: Download GSE55235 from GEO, inspect sample metadata
#
# GSE55235 is Affymetrix HG-U133A, contains normal/OA/RA samples from
# synovial membrane. We filter to OA-vs-Normal synovium only.
#
# OUTPUT : data/raw/GSE55235_expr.rds  (expression matrix)
#          data/raw/GSE55235_pheno.csv (sample metadata)
# =============================================================================

# Install GEOquery if needed
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  cat("Installing GEOquery...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
})

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("STAGE 1: GSE55235 DOWNLOAD & INSPECT\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# Download
# -----------------------------------------------------------------------------
cat("=== DOWNLOADING GSE55235 ===\n")
cat("This can take 1-3 minutes depending on connection...\n")

dest_dir <- file.path(PROJECT_ROOT, "data/raw/GSE55235_raw")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
gse <- tryCatch({
  getGEO("GSE55235", GSEMatrix = TRUE, destdir = dest_dir,
         AnnotGPL = TRUE)
}, error = function(e) {
  cat("ERROR downloading:", conditionMessage(e), "\n")
  NULL
})
t_end <- Sys.time()

if (is.null(gse)) {
  cat("Download failed. Check network. BAILOUT.\n")
  stop("Cannot proceed without data")
}

cat(sprintf("Download completed in %.1f seconds\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))
cat("ExpressionSet(s) returned:", length(gse), "\n")

# GSE55235 has typically 1 platform (HG-U133A)
eset <- gse[[1]]

cat("Expression matrix dims:", nrow(exprs(eset)), "x", ncol(exprs(eset)), "\n")
cat("Platform:", annotation(eset), "\n\n")

# -----------------------------------------------------------------------------
# Inspect phenotype data
# -----------------------------------------------------------------------------
pheno <- pData(eset)

cat("=== PHENOTYPE COLUMNS AVAILABLE ===\n")
cat(paste(colnames(pheno), collapse = "\n"), "\n\n")

# The critical columns are usually title, source_name_ch1, or characteristics_ch1.*
cat("=== SAMPLE TITLES ===\n")
print(pheno$title)
cat("\n")

cat("=== SOURCE NAMES ===\n")
print(pheno$source_name_ch1)
cat("\n")

# Look for tissue and disease info in characteristics columns
char_cols <- grep("characteristics", colnames(pheno), value = TRUE)
if (length(char_cols) > 0) {
  cat("=== CHARACTERISTICS COLUMNS ===\n")
  for (col in char_cols) {
    cat(sprintf("-- %s --\n", col))
    print(unique(pheno[[col]]))
    cat("\n")
  }
}

# -----------------------------------------------------------------------------
# Attempt to parse condition from title (first heuristic)
# -----------------------------------------------------------------------------
cat("=== CONDITION INFERENCE FROM TITLE ===\n")
titles_lower <- tolower(as.character(pheno$title))

# Tag each sample
pheno$Condition_inferred <- dplyr::case_when(
  grepl("normal|healthy|control", titles_lower) ~ "Normal",
  grepl("osteoarth|^oa|[^a-z]oa[^a-z]", titles_lower) ~ "OA",
  grepl("rheumatoid|[^a-z]ra[^a-z]", titles_lower) ~ "RA",
  TRUE ~ "UNCLEAR"
)

cat("Condition breakdown:\n")
print(table(pheno$Condition_inferred))
cat("\n")

# -----------------------------------------------------------------------------
# Save inspected data for next stage
# -----------------------------------------------------------------------------
saveRDS(eset, file.path(PROJECT_ROOT, "data/raw/GSE55235_eset.rds"))
write.csv(pheno,
          file.path(PROJECT_ROOT, "data/raw/GSE55235_pheno.csv"),
          row.names = TRUE)

cat("Saved: data/raw/GSE55235_eset.rds\n")
cat("Saved: data/raw/GSE55235_pheno.csv\n\n")

cat("============================================\n")
cat("STAGE 1 COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")