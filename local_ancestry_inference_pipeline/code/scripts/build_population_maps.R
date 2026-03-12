#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("kgp", quietly = TRUE)) {
    install.packages("kgp", repos = "https://cloud.r-project.org")
  }
  library(kgp)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  script_path <- normalizePath(file.path(getwd(), "code", "scripts", "build_population_maps.R"), mustWork = FALSE)
}
root_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
meta_dir <- file.path(root_dir, "code", "meta")

if (!dir.exists(meta_dir)) {
  stop("Missing directory: ", meta_dir, "\nRun build_sample_sets.sh first.")
}

all_samples_file <- file.path(meta_dir, "samples.all_3202.tsv")
ref_strict_file <- file.path(meta_dir, "reference.strict.tsv")
ref_child_file <- file.path(meta_dir, "reference.child_only.tsv")

for (p in c(all_samples_file, ref_strict_file, ref_child_file)) {
  if (!file.exists(p)) stop("Missing required file: ", p)
}

load_kgp_metadata <- function() {
  ns <- asNamespace("kgp")

  # Preferred objects observed across kgp versions.
  for (obj in c("kgpe", "kgp")) {
    if (exists(obj, envir = ns, inherits = FALSE)) {
      return(as.data.frame(get(obj, envir = ns), stringsAsFactors = FALSE))
    }
  }

  # Fallback: try data() loading into an isolated env.
  e <- new.env(parent = emptyenv())
  for (obj in c("kgpe", "kgp")) {
    suppressWarnings(try(data(list = obj, package = "kgp", envir = e), silent = TRUE))
    if (exists(obj, envir = e, inherits = FALSE)) {
      return(as.data.frame(get(obj, envir = e), stringsAsFactors = FALSE))
    }
  }

  stop("Could not find a metadata table in the kgp package (tried objects: kgpe, kgp).")
}

normalize_superpopulation <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA", "NAN", "NULL")] <- NA_character_

  alias <- c(
    "AFRICA" = "AFR",
    "AFRICAN" = "AFR",
    "AMERICA" = "AMR",
    "AMERICAS" = "AMR",
    "AMERICAN" = "AMR",
    "ADMIXED AMERICAN" = "AMR",
    "ADMIXED_AMERICAN" = "AMR",
    "LATIN AMERICAN" = "AMR",
    "LATIN_AMERICAN" = "AMR",
    "EAST ASIAN" = "EAS",
    "EAST_ASIAN" = "EAS",
    "EUROPE" = "EUR",
    "EUROPEAN" = "EUR",
    "SOUTH ASIAN" = "SAS",
    "SOUTH_ASIAN" = "SAS"
  )

  idx <- y %in% names(alias)
  y[idx] <- alias[y[idx]]

  allowed <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  y[!is.na(y) & !(y %in% allowed)] <- NA_character_
  y
}

superpop_from_population <- function(pop) {
  p <- toupper(trimws(as.character(pop)))
  pop_map <- c(
    "ACB" = "AFR", "ASW" = "AFR", "ESN" = "AFR", "GWD" = "AFR", "LWK" = "AFR", "MSL" = "AFR", "YRI" = "AFR",
    "CDX" = "EAS", "CHB" = "EAS", "CHS" = "EAS", "JPT" = "EAS", "KHV" = "EAS",
    "CEU" = "EUR", "FIN" = "EUR", "GBR" = "EUR", "IBS" = "EUR", "TSI" = "EUR",
    "BEB" = "SAS", "GIH" = "SAS", "ITU" = "SAS", "PJL" = "SAS", "STU" = "SAS",
    "CLM" = "AMR", "MXL" = "AMR", "PEL" = "AMR", "PUR" = "AMR"
  )
  out <- unname(pop_map[p])
  out[is.na(out)] <- NA_character_
  out
}

meta <- load_kgp_metadata()
nm <- names(meta)
nm_l <- tolower(nm)

pick_col <- function(candidates) {
  idx <- match(candidates, nm_l)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

id_idx <- pick_col(c("id", "sampleid", "sample_id", "sample", "individual", "iid"))
pop_idx <- pick_col(c("pop", "population", "population_code"))
sup_idx <- pick_col(c("superpop", "super_pop", "super.population", "superpopulation", "super.population_code", "reg", "region"))

if (is.na(id_idx) || is.na(pop_idx)) {
  stop(
    "Could not detect required sample/population columns in kgp metadata.\n",
    "Available columns: ", paste(nm, collapse = ", ")
  )
}

sup_raw <- if (!is.na(sup_idx)) as.character(meta[[sup_idx]]) else rep(NA_character_, nrow(meta))
sup_norm <- normalize_superpopulation(sup_raw)
sup_fill <- superpop_from_population(meta[[pop_idx]])
needs_fill <- is.na(sup_norm) | !nzchar(sup_norm)
sup_norm[needs_fill] <- sup_fill[needs_fill]

meta_clean <- data.frame(
  sample_id = as.character(meta[[id_idx]]),
  population = as.character(meta[[pop_idx]]),
  superpopulation = as.character(sup_norm),
  stringsAsFactors = FALSE
)

meta_clean <- meta_clean[!is.na(meta_clean$sample_id) & nzchar(meta_clean$sample_id), ]
meta_clean <- meta_clean[!duplicated(meta_clean$sample_id), ]

all_samples <- read.table(all_samples_file, sep = "\t", stringsAsFactors = FALSE, col.names = "sample_id")
ref_strict <- read.table(ref_strict_file, sep = "\t", stringsAsFactors = FALSE, col.names = "sample_id")
ref_child <- read.table(ref_child_file, sep = "\t", stringsAsFactors = FALSE, col.names = "sample_id")

idx_all <- match(all_samples$sample_id, meta_clean$sample_id)
map_all <- data.frame(
  sample_id = all_samples$sample_id,
  population = meta_clean$population[idx_all],
  superpopulation = meta_clean$superpopulation[idx_all],
  stringsAsFactors = FALSE
)
missing_n <- sum(is.na(map_all$population) | is.na(map_all$superpopulation))
if (missing_n > 0) {
  missing_ids <- map_all$sample_id[is.na(map_all$population) | is.na(map_all$superpopulation)]
  writeLines(missing_ids, con = file.path(meta_dir, "missing_population_labels.txt"))
  stop(
    "Missing population/superpopulation labels for ", missing_n, " samples.\n",
    "See: ", file.path(meta_dir, "missing_population_labels.txt")
  )
}

# Canonical full maps for all 3202 samples.
write.table(
  map_all[, c("sample_id", "population")],
  file = file.path(meta_dir, "sample_to_population.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  map_all[, c("sample_id", "superpopulation")],
  file = file.path(meta_dir, "sample_to_superpopulation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  map_all[, c("sample_id", "population", "superpopulation")],
  file = file.path(meta_dir, "sample_to_pop_superpop.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  map_all,
  file = file.path(meta_dir, "kgpe_3202_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

subset_map <- function(ref_df, label_col) {
  idx <- match(ref_df$sample_id, map_all$sample_id)
  out <- data.frame(
    sample_id = ref_df$sample_id,
    label = map_all[[label_col]][idx],
    stringsAsFactors = FALSE
  )
  if (any(is.na(out$label) | !nzchar(out$label))) {
    bad <- out$sample_id[is.na(out$label) | !nzchar(out$label)]
    writeLines(bad, con = file.path(meta_dir, paste0("missing_labels_in_", label_col, ".txt")))
    stop("Missing labels detected for ", length(bad), " reference samples in ", label_col)
  }
  out
}

write_flare <- function(df, out_file) {
  # FLARE ref-panel uses whitespace-delimited: <sample_id> <panel_name>
  write.table(df, file = out_file, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_tsv_map <- function(df, out_file) {
  # RFMix sample maps are tab-delimited: <sample_id>\t<label>
  write.table(df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Strict reference maps.
strict_pop <- subset_map(ref_strict, "population")
strict_sup <- subset_map(ref_strict, "superpopulation")

# Child-only reference maps.
child_pop <- subset_map(ref_child, "population")
child_sup <- subset_map(ref_child, "superpopulation")

# FLARE
write_flare(strict_pop, file.path(meta_dir, "flare_ref_panel.strict.population.txt"))
write_flare(strict_sup, file.path(meta_dir, "flare_ref_panel.strict.superpopulation.txt"))
write_flare(child_pop, file.path(meta_dir, "flare_ref_panel.child_only.population.txt"))
write_flare(child_sup, file.path(meta_dir, "flare_ref_panel.child_only.superpopulation.txt"))

# RFMix
write_tsv_map(strict_pop, file.path(meta_dir, "rfmix_sample_map.strict.population.tsv"))
write_tsv_map(strict_sup, file.path(meta_dir, "rfmix_sample_map.strict.superpopulation.tsv"))
write_tsv_map(child_pop, file.path(meta_dir, "rfmix_sample_map.child_only.population.tsv"))
write_tsv_map(child_sup, file.path(meta_dir, "rfmix_sample_map.child_only.superpopulation.tsv"))

summary_df <- data.frame(
  metric = c(
    "n_all_samples",
    "n_ref_strict",
    "n_ref_child_only",
    "n_unique_population_labels",
    "n_unique_superpopulation_labels"
  ),
  value = c(
    nrow(map_all),
    nrow(ref_strict),
    nrow(ref_child),
    length(unique(map_all$population)),
    length(unique(map_all$superpopulation))
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = file.path(meta_dir, "summary.population_maps.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Wrote population/superpopulation maps and FLARE/RFMix panel maps to:\n", meta_dir, "\n", sep = "")
