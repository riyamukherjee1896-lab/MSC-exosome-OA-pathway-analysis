# =============================================================================
# 08_miRNA_target_mapping.R
# -----------------------------------------------------------------------------
# Block C: Map cargo miRNAs to experimentally validated target genes
#
# STRATEGY:
#   Query multiMiR (Bioconductor) using MATURE miRNA names (hsa-miR-X-5p/3p).
#   Primary-transcript queries return incomplete results for many miRNAs due
#   to multiMiR's internal name mapping; mature names resolve this.
#   Filter to validated interactions, excluding Non-Functional MTI (validated
#   non-interactions). Flag Tier-1 strict subset (Functional MTI only) for
#   sensitivity analysis.
#
# DATABASE:
#   miRTarBase only (TarBase server was intermittently unavailable; miRTarBase
#   alone provides 240K+ validated interactions and is the field-standard gold
#   standard for validated miRNA-target data).
#
# INPUT  : data/processed/cargo_miRNAs.rds
#
# OUTPUT : data/processed/mirna_targets_validated.rds
#          data/processed/multimir_requery_raw.rds
#          tables/table_C1_cargo_coverage.csv
#          tables/table_C2_top_targeted_genes.csv
#          tables/table_C3_top_targeting_miRNAs.csv
#          logs/08_target_mapping_sessionInfo.txt
# =============================================================================

# Namespace-safe: explicit dplyr:: calls throughout to avoid AnnotationDbi /
# clusterProfiler / MASS conflicts with select() and filter()
suppressPackageStartupMessages({
  library(dplyr)
})

# Install multiMiR if needed
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  cat("multiMiR not found. Installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("multiMiR", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages(library(multiMiR))

PROJECT_ROOT <- "D:/Research Work/Final Version Paper 3/new paper V1/Paper3_clean"

cat("============================================\n")
cat("BLOCK C: miRNA Target Mapping\n")
cat("============================================\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# STEP 1: Load cargo
# =============================================================================
cargo_data <- readRDS(file.path(PROJECT_ROOT, "data/processed/cargo_miRNAs.rds"))
cargo_mirnas <- cargo_data$cargo_mirnas

cat("Cargo miRNAs loaded:", length(cargo_mirnas), "\n")
cat("First 5:", paste(head(cargo_mirnas, 5), collapse = ", "), "\n\n")

# =============================================================================
# STEP 2: Generate candidate mature names for each cargo primary
#
# Mapping rules:
#   hsa-mir-21        -> hsa-miR-21-5p, hsa-miR-21-3p, hsa-miR-21 (no-arm)
#   hsa-let-7b        -> hsa-let-7b-5p, hsa-let-7b-3p, hsa-let-7b (no-arm)
#   hsa-mir-199a-1    -> hsa-miR-199a-5p, hsa-miR-199a-3p (locus suffix DROPPED)
#   hsa-mir-486-2     -> hsa-miR-486-5p, hsa-miR-486-3p (locus suffix DROPPED)
# =============================================================================

generate_mature_forms <- function(primary) {
  is_let <- grepl("^hsa-let-", primary)
  prefix <- if (is_let) "hsa-let-" else "hsa-miR-"
  
  # Strip the "hsa-mir-" / "hsa-let-" prefix
  base <- sub("^hsa-(mir|let)-", "", primary)
  
  # Drop trailing locus suffix (-1, -2, -3) which belongs to primary only
  mature_base <- sub("-[0-9]+$", "", base)
  
  # Generate candidate mature names
  forms <- c(
    paste0(prefix, mature_base),           # no-arm form (e.g., hsa-miR-484)
    paste0(prefix, mature_base, "-5p"),    # 5p arm
    paste0(prefix, mature_base, "-3p"),    # 3p arm
    paste0(prefix, base)                    # also try case-changed primary
  )
  unique(forms)
}

mature_forms_list <- lapply(cargo_mirnas, generate_mature_forms)
names(mature_forms_list) <- cargo_mirnas
all_mature_forms <- unique(unlist(mature_forms_list))

cat("=== MATURE NAME GENERATION ===\n")
cat("Generated", length(all_mature_forms), "unique candidate mature forms\n")
cat("Example mappings:\n")
for (ex in c("hsa-mir-21", "hsa-let-7b", "hsa-mir-199a-1", "hsa-mir-486")) {
  if (ex %in% names(mature_forms_list)) {
    cat(sprintf("  %-20s -> %s\n", ex,
                paste(mature_forms_list[[ex]], collapse = ", ")))
  }
}
cat("\n")

# =============================================================================
# STEP 3: Query multiMiR (miRTarBase only)
#
# Note: this script queries miRTarBase exclusively because:
#   (a) TarBase via multiMiR returned server errors during development
#   (b) miRTarBase alone provides ~240K validated interactions, sufficient
#   (c) miRTarBase is the field-standard for validated miRNA-target data
#
# If TarBase data is needed later, an additional query can be added.
# =============================================================================

cat("=== QUERYING miRTarBase (validated) ===\n")
cat("Querying", length(all_mature_forms), "mature names against miRTarBase...\n")
cat("Expected runtime: 60-120 seconds\n\n")

t_start <- Sys.time()
targets_raw <- tryCatch({
  get_multimir(
    org   = "hsa",
    mirna = all_mature_forms,
    table = "mirtarbase"
  )
}, error = function(e) {
  cat("ERROR in multiMiR query:", conditionMessage(e), "\n")
  cat("Try running the query again in 10-15 minutes (server may be temporarily unavailable).\n")
  NULL
})
t_end <- Sys.time()
cat(sprintf("Query completed in %.1f seconds\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))

if (is.null(targets_raw)) {
  stop("Cannot proceed without multiMiR data.")
}

targets_df <- targets_raw@data

cat("\n=== RAW QUERY RESULTS ===\n")
cat("Total miRTarBase interactions returned:", nrow(targets_df), "\n")
cat("Unique mature miRNAs:", length(unique(targets_df$mature_mirna_id)), "\n")
cat("Unique target genes:", length(unique(targets_df$target_symbol)), "\n\n")

cat("Support types in raw data:\n")
print(table(targets_df$support_type))
cat("\n")

# Save raw return for reproducibility
saveRDS(list(
  targets_df        = targets_df,
  mature_forms_list = mature_forms_list,
  all_mature_forms  = all_mature_forms,
  query_time        = t_end - t_start,
  source_database   = "mirtarbase_only"
),
file.path(PROJECT_ROOT, "data/processed/multimir_requery_raw.rds"))

# =============================================================================
# STEP 4: Filter
#
# EXCLUDE: Non-Functional MTI entries (these are validated NON-interactions,
#          not evidence of targeting).
# INCLUDE: Everything else returned by table = "validated".
# FLAG:    Functional MTI + Functional MTI (Weak) as Tier-1 strict
#          (for sensitivity analysis, not the primary filter).
# =============================================================================

validated <- targets_df %>%
  dplyr::filter(!support_type %in% c("Non-Functional MTI",
                                     "Non-Functional MTI (Weak)")) %>%
  dplyr::filter(!is.na(target_symbol), target_symbol != "") %>%
  dplyr::select(database, mature_mirna_id, mature_mirna_acc,
                target_symbol, target_entrez,
                experiment, support_type, pubmed_id) %>%
  dplyr::distinct()

# Flag tier-1 strict for sensitivity
validated$tier1_strict <- validated$support_type %in%
  c("Functional MTI", "Functional MTI (Weak)")

cat("=== FILTERED EVIDENCE ===\n")
cat("After removing Non-Functional MTI and empty targets:\n")
cat("  Total validated interactions:", nrow(validated), "\n")
cat("  Unique mature miRNAs:", length(unique(validated$mature_mirna_id)), "\n")
cat("  Unique target genes:", length(unique(validated$target_symbol)), "\n")
cat("  Tier-1 strict subset (Functional MTI only):",
    sum(validated$tier1_strict), "\n\n")

# =============================================================================
# STEP 5: Map cargo primary names to mature names in validated data
# =============================================================================

validated_mature_ids <- unique(validated$mature_mirna_id)
validated_mature_lc  <- tolower(validated_mature_ids)

cargo_coverage <- data.frame(
  cargo_miRNA              = cargo_mirnas,
  mature_forms_in_evidence = NA_character_,
  n_validated_targets      = 0L,
  n_tier1_strict_targets   = 0L,
  stringsAsFactors = FALSE
)

for (i in seq_along(cargo_mirnas)) {
  cm    <- cargo_mirnas[i]
  cm_lc <- tolower(cm)
  
  # Build normalized target forms by rule:
  #   "hsa-mir-X" -> look for "hsa-mir-<mature_base>" and "hsa-mir-<mature_base>-..."
  is_let <- startsWith(cm_lc, "hsa-let-")
  prefix <- if (is_let) "hsa-let-" else "hsa-mir-"
  base        <- sub("^hsa-(mir|let)-", "", cm_lc)
  mature_base <- sub("-[0-9]+$", "", base)   # drop locus suffix
  
  target_exact  <- paste0(prefix, mature_base)          # exact (no arm)
  target_prefix <- paste0(prefix, mature_base, "-")     # prefix before arm
  
  # Match: exact OR starts with prefix+hyphen (hyphen prevents miR-21 -> miR-21b)
  matches_exact  <- which(validated_mature_lc == target_exact)
  matches_prefix <- which(startsWith(validated_mature_lc, target_prefix))
  matched_idx    <- unique(c(matches_exact, matches_prefix))
  matched_mature <- validated_mature_ids[matched_idx]
  
  cargo_coverage$mature_forms_in_evidence[i] <- paste(matched_mature,
                                                      collapse = "|")
  cargo_coverage$n_validated_targets[i] <-
    sum(validated$mature_mirna_id %in% matched_mature)
  cargo_coverage$n_tier1_strict_targets[i] <-
    sum(validated$mature_mirna_id %in% matched_mature & validated$tier1_strict)
}

n_covered     <- sum(cargo_coverage$n_validated_targets > 0)
n_not_covered <- sum(cargo_coverage$n_validated_targets == 0)

cat("=== CARGO COVERAGE ===\n")
cat(sprintf("Cargo miRNAs with >= 1 validated target: %d / %d (%.1f%%)\n",
            n_covered, length(cargo_mirnas),
            100 * n_covered / length(cargo_mirnas)))
cat(sprintf("Cargo miRNAs with NO validated target: %d\n\n", n_not_covered))

cat("Distribution of validated targets per covered cargo miRNA:\n")
print(summary(cargo_coverage$n_validated_targets[cargo_coverage$n_validated_targets > 0]))

cat("\nTop 20 cargo miRNAs by number of validated targets:\n")
top20 <- cargo_coverage %>%
  dplyr::arrange(dplyr::desc(n_validated_targets)) %>%
  head(20)
print(top20[, c("cargo_miRNA", "n_validated_targets",
                "n_tier1_strict_targets")], row.names = FALSE)

write.csv(cargo_coverage,
          file.path(PROJECT_ROOT, "tables/table_C1_cargo_coverage.csv"),
          row.names = FALSE)
cat("\nSaved: tables/table_C1_cargo_coverage.csv\n\n")

# =============================================================================
# STEP 6: Canonical MSC-exosome miRNA sanity check
# =============================================================================

cat("=== CANONICAL MSC miRNA SANITY CHECK ===\n")
canonical <- c("hsa-mir-21", "hsa-mir-100", "hsa-mir-143", "hsa-mir-145",
               "hsa-let-7b", "hsa-let-7c", "hsa-mir-181a-1",
               "hsa-mir-26a-1", "hsa-mir-27a", "hsa-mir-16-1")
for (m in canonical) {
  row <- cargo_coverage %>% dplyr::filter(cargo_miRNA == m)
  if (nrow(row) == 1) {
    cat(sprintf("  %-20s  validated: %5d   tier1-strict: %5d\n",
                m, row$n_validated_targets, row$n_tier1_strict_targets))
  }
}
cat("\n(Expected: canonical MSC miRNAs should have hundreds to thousands of targets.\n")
cat(" If any show 0, the mature-name mapping is broken.)\n\n")

# =============================================================================
# STEP 7: Build cargo-target edge list
# =============================================================================

cargo_mature_set <- unique(unlist(
  strsplit(cargo_coverage$mature_forms_in_evidence, "\\|")
))
cargo_mature_set <- cargo_mature_set[cargo_mature_set != ""]

cargo_edges <- validated %>%
  dplyr::filter(mature_mirna_id %in% cargo_mature_set)

# Build mature -> primary lookup
mature_to_primary <- data.frame(mature  = character(),
                                primary = character(),
                                stringsAsFactors = FALSE)
for (i in seq_along(cargo_mirnas)) {
  forms <- strsplit(cargo_coverage$mature_forms_in_evidence[i], "\\|")[[1]]
  forms <- forms[forms != ""]
  if (length(forms) > 0) {
    mature_to_primary <- rbind(mature_to_primary,
                               data.frame(mature  = forms,
                                          primary = cargo_mirnas[i],
                                          stringsAsFactors = FALSE))
  }
}

cargo_edges <- cargo_edges %>%
  dplyr::left_join(mature_to_primary,
                   by = c("mature_mirna_id" = "mature"))

cat("=== CARGO-TARGET EDGE LIST ===\n")
cat("Total cargo-target edges:", nrow(cargo_edges), "\n")
cat("Unique cargo miRNAs represented:", length(unique(cargo_edges$primary)), "\n")
cat("Unique target genes:", length(unique(cargo_edges$target_symbol)), "\n\n")

# =============================================================================
# STEP 8: Top-targeted genes
# =============================================================================

target_frequency <- cargo_edges %>%
  dplyr::group_by(target_symbol) %>%
  dplyr::summarise(
    n_mature_miRNAs_targeting = dplyr::n_distinct(mature_mirna_id),
    n_primary_cargo_targeting = dplyr::n_distinct(primary),
    n_tier1_strict_evidence   = sum(tier1_strict),
    mature_miRNAs  = paste(sort(unique(mature_mirna_id)), collapse = "; "),
    primary_cargo  = paste(sort(unique(primary)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_mature_miRNAs_targeting))

cat("=== TOP 25 MOST-TARGETED GENES ===\n")
print(head(target_frequency[, c("target_symbol",
                                "n_cargo_miRNAs_targeting",
                                "n_tier1_strict_evidence")], 25),
      row.names = FALSE)

write.csv(target_frequency,
          file.path(PROJECT_ROOT, "tables/table_C2_top_targeted_genes.csv"),
          row.names = FALSE)
cat("\nSaved: tables/table_C2_top_targeted_genes.csv\n\n")

# =============================================================================
# STEP 9: Anchor gene coverage (paper-critical result)
# =============================================================================

anchor_check <- c("ALDH2", "SOD2", "GPX3", "GPX4", "NQO1", "HMOX1", "NFE2L2",
                  "PRDX1", "ACSL4", "SLC3A2", "TFRC", "AGO2", "SRSF3",
                  "DDX17", "DGCR8", "DROSHA", "DICER1")

anchor_hits <- target_frequency %>%
  dplyr::filter(target_symbol %in% anchor_check) %>%
  dplyr::select(target_symbol, n_cargo_miRNAs_targeting,
                n_tier1_strict_evidence, targeting_miRNAs)

cat("=== ANCHOR GENE COVERAGE IN CARGO TARGETS ===\n")
cat("Anchor genes present as cargo targets:\n")
print(anchor_hits[, c("target_symbol", "n_cargo_miRNAs_targeting",
                      "n_tier1_strict_evidence")],
      row.names = FALSE)

not_hit <- setdiff(anchor_check, anchor_hits$target_symbol)
cat("\nAnchor genes NOT targeted by any cargo miRNA:\n")
if (length(not_hit) > 0) {
  cat(" ", paste(not_hit, collapse=", "), "\n")
} else {
  cat("  (none)\n")
}
cat("\n")

# =============================================================================
# STEP 10: Top targeting miRNAs table
# =============================================================================

mirna_targeting <- cargo_edges %>%
  dplyr::group_by(primary) %>%
  dplyr::summarise(
    n_targets = dplyr::n_distinct(target_symbol),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_targets))

write.csv(mirna_targeting,
          file.path(PROJECT_ROOT, "tables/table_C3_top_targeting_miRNAs.csv"),
          row.names = FALSE)

# =============================================================================
# STEP 11: Save consolidated results
# =============================================================================

target_results <- list(
  cargo_coverage   = cargo_coverage,
  cargo_edges      = cargo_edges,
  target_frequency = target_frequency,
  mirna_targeting  = mirna_targeting,
  anchor_hits      = anchor_hits,
  n_cargo_query    = length(cargo_mirnas),
  n_cargo_covered  = n_covered,
  total_edges      = nrow(cargo_edges),
  unique_targets   = length(unique(cargo_edges$target_symbol)),
  mature_forms_list = mature_forms_list,
  evidence_filter  = "multiMiR miRTarBase validated (mature-name query, Non-Functional MTI excluded)",
  tier1_filter     = "Functional MTI + Functional MTI (Weak) — sensitivity subset"
)

saveRDS(target_results,
        file.path(PROJECT_ROOT, "data/processed/mirna_targets_validated.rds"))
cat("Saved: data/processed/mirna_targets_validated.rds\n\n")

# =============================================================================
# STEP 12: File verification
# =============================================================================

cat("=== FILE EXISTENCE VERIFICATION ===\n")
files_to_check <- c(
  "data/processed/mirna_targets_validated.rds",
  "data/processed/multimir_requery_raw.rds",
  "tables/table_C1_cargo_coverage.csv",
  "tables/table_C2_top_targeted_genes.csv",
  "tables/table_C3_top_targeting_miRNAs.csv"
)
for (f in files_to_check) {
  full <- file.path(PROJECT_ROOT, f)
  exists_flag <- file.exists(full)
  size_mb <- if (exists_flag) round(file.info(full)$size / 1024^2, 2) else NA
  cat(sprintf("  %-55s  %s  %s MB\n", f,
              ifelse(exists_flag, "OK", "MISSING"),
              ifelse(is.na(size_mb), "-", as.character(size_mb))))
}

# =============================================================================
# Session info
# =============================================================================
sink(file.path(PROJECT_ROOT, "logs/08_target_mapping_sessionInfo.txt"))
cat("Script: 08_miRNA_target_mapping.R\n")
cat("Run:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\n============================================\n")
cat("BLOCK C COMPLETE\n")
cat("============================================\n")
cat("Finished:", format(Sys.time()), "\n")