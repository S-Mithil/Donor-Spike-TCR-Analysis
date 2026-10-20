# ==============================================================================
# Pt1006_TCR_analysis_complete.R
# ==============================================================================
# Manuscript analysis code for longitudinal TCRbeta repertoire analyses in a
# liver-transplant recipient with SARS-CoV-2 vaccine-associated TCMR.
#
# This single script contains the custom R analyses used for the principal
# clonotype-level analyses in Figure 2. Each analysis is separated into a
# clearly labeled section and annotated with the corresponding figure panel(s).
#
# ANALYSIS MAP
#   Section 1  - Project setup, packages, input files and thresholds
#   Section 2  - Helper functions
#   Section 3  - Prepare cleaned TCRbeta repertoire matrix
#   Section 4  - Define donor-reactive (DRTCC) and spike-reactive (SRTCC) clones
#                Relevant to Fig. 2a-e and Fig. 2h
#   Section 5  - 12 MLR x spike pairwise overlap comparisons
#                Fig. 2b,d
#   Section 6  - Cumulative productive frequency of shared clonotypes
#                Fig. 2c,e
#   Section 7  - Dominant circulating CD8 clonotypes / alluvial analysis
#                Fig. 2f
#   Section 8  - Biopsy representation by clonotype class
#                Fig. 2h
#   Section 9  - Productive frequency and abundance-rank analysis
#                Fig. 2i,j
#   Section 10 - Save package/version information for reproducibility
#
# INPUT DATA
# The script expects the following files in data/raw/:
#   1. Pt1006_All_Rearrangements_Ambiguous_removed_plus_30M_54M.csv
#   2. CD4_Direct_PreDRTCCs.csv
#   3. CD8_Direct_PreDRTCCs.csv
#   4. Pt1006_Spike_Reactive_Clones_30M_33-7M_54M.xlsx
#
# The complete TCRbeta sequencing dataset can be deposited separately in
# ImmuneACCESS. This script can be distributed without patient-identifiable
# raw clinical data; it operates on the processed repertoire-level inputs.
#
# IMPORTANT ANALYTIC DEFINITIONS
# - Post-transplant DRTCC: productive frequency >=1e-5 in CFSE-low MLR cells
#   and >2-fold enrichment relative to the designated unstimulated repertoire.
# - SRTCC assignments are imported from the curated timepoint-specific spike-
#   reactive workbook generated using the study's predefined enrichment rules.
# - Fig. 2f: union of the top 20 circulating CD8 clonotypes at each of five
#   longitudinal blood time points; retain clones with mean frequency >=0.75%
#   across the five samples. Selection is independent of DRTCC/SRTCC status.
# - Fig. 2h CRTCC hierarchy: any shared clone with pre-transplant donor-reactive
#   assignment is classified as Pre-CRTCC. Post-CRTCC therefore contains only
#   shared clones lacking pre-transplant donor-reactive assignment.
# - Fig. 2j: rank 1 is the most abundant detected clonotype; identical
#   productive frequencies share the same rank (competition ranking); clones
#   with zero frequency are ND and are not ranked.
#
# Frequencies in the repertoire input are fractions (0-1). Values reported in
# figures/source-data tables are converted to percentages where indicated.
# ==============================================================================

# ==============================================================================
# SECTION 1. PROJECT SETUP, PACKAGES, INPUT FILES AND THRESHOLDS
# Figures: all downstream analyses
# ==============================================================================

required_packages <- c(
  "readr", "readxl", "dplyr", "tidyr", "ggplot2", "ggalluvial",
  "RColorBrewer", "openxlsx", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(RColorBrewer)
  library(openxlsx)
  library(scales)
})

# Set TCMR_PROJECT_DIR to the repository root, or run the script from that root.
# Example:
# Sys.setenv(TCMR_PROJECT_DIR = "/path/to/donor-spike-tcr-repertoire-analysis")
project_dir <- Sys.getenv("TCMR_PROJECT_DIR", unset = getwd())

raw_dir       <- file.path(project_dir, "data", "raw")
processed_dir <- file.path(project_dir, "data", "processed")
results_dir   <- file.path(project_dir, "results")
figures_dir   <- file.path(results_dir, "figures")
tables_dir    <- file.path(results_dir, "tables")

for (d in c(processed_dir, results_dir, figures_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

repertoire_file <- file.path(
  raw_dir, "Pt1006_All_Rearrangements_Ambiguous_removed_plus_30M_54M.csv"
)
pre_cd4_drtcc_file <- file.path(raw_dir, "CD4_Direct_PreDRTCCs.csv")
pre_cd8_drtcc_file <- file.path(raw_dir, "CD8_Direct_PreDRTCCs.csv")
spike_reactive_file <- file.path(
  raw_dir, "Pt1006_Spike_Reactive_Clones_30M_33-7M_54M.xlsx"
)

# Functional-clonotype thresholds
min_functional_frequency <- 1e-5
fold_enrichment_cutoff <- 2

# Fig. 2h comparator-clone definition
comparator_fold_expansion <- 5
comparator_min_templates <- 8
comparator_min_frequency <- 1e-5

# Productive template counts from sequencing summary/source data
cd8_30m_blood_total_templates <- 43614
cd8_fcb_total_templates <- 13765

# Fig. 2f selection criteria
alluvial_top_n_each_timepoint <- 20
alluvial_mean_percent_cutoff <- 0.75

# Fig. 2i,j representative CRTCCs
representative_crtccs <- c(
  "CASSLGRLRSPLHF",
  "CASWSNRGQETQYF",
  "CASSEFAPRGSGDYEQYF",
  "CASSQSTGDWDTQYF",
  "CASSSGGLEGEQYV"
)

# Fig. 2f highlighted donor-spike shared clonotypes
special_clone_colors <- c(
  "CASSLGRLRSPLHF" = "#E31A1C",
  "CASWSNRGQETQYF" = "#1F78B4",
  "CASSSGGLEGEQYV" = "#6A3D9A"
)

# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# Figures: all downstream analyses
# ==============================================================================

clean_cdr3 <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

as_numeric_zero <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

assert_columns <- function(df, cols, object_name = "data frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

read_cdr3_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (!"CDR3" %in% names(x)) names(x)[1] <- "CDR3"
  unique(stats::na.omit(clean_cdr3(x$CDR3)))
}

truthy <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  if (is.numeric(x)) return(replace(x == 1, is.na(x), FALSE))
  z <- toupper(trimws(as.character(x)))
  replace(z %in% c("TRUE", "T", "1", "YES", "Y"), is.na(z), FALSE)
}

make_post_drtcc <- function(df, cfse_col, unstim_col,
                            min_freq = 1e-5, fold = 2) {
  assert_columns(df, c("CDR3", cfse_col, unstim_col), "repertoire matrix")
  df |>
    dplyr::filter(
      .data[[cfse_col]] >= min_freq,
      .data[[unstim_col]] == 0 |
        .data[[cfse_col]] > fold * .data[[unstim_col]]
    ) |>
    dplyr::pull(CDR3) |>
    unique()
}

competition_rank <- function(values) {
  # Rank 1 = most abundant detected clonotype.
  # Ties share the same rank; zero/absent values are ND (NA internally).
  out <- rep(NA_integer_, length(values))
  keep <- !is.na(values) & values > 0
  out[keep] <- rank(-values[keep], ties.method = "min")
  out
}

safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*:/\\\\?]", "_", x)
  substr(x, 1, 31)
}

write_csv_clean <- function(df, path) {
  readr::write_csv(df, path, na = "")
  message("Wrote: ", path)
}

# ==============================================================================
# SECTION 3. PREPARE CLEANED TCRbeta REPERTOIRE MATRIX
# Figures: prerequisite for Fig. 2 analyses
# ==============================================================================
# The source repertoire matrix has already undergone the laboratory-contaminant
# and CD4/CD8 lineage-ambiguity filtering described in the Methods. This section
# removes blank/header rows, converts frequency columns to numeric values and
# collapses any exact duplicate CDR3 rows that remain.

if (!file.exists(repertoire_file)) stop("Input not found: ", repertoire_file)

rep_raw <- read_csv(repertoire_file, show_col_types = FALSE)
if (!"CDR3" %in% names(rep_raw)) names(rep_raw)[1] <- "CDR3"
rep_raw$CDR3 <- clean_cdr3(rep_raw$CDR3)

freq_cols <- names(rep_raw)[grepl("TCRB$", names(rep_raw))]
if (length(freq_cols) == 0) stop("No repertoire columns ending in 'TCRB' found.")

rep_raw <- rep_raw |>
  filter(!is.na(CDR3), CDR3 != "CDR3")
rep_raw[freq_cols] <- lapply(rep_raw[freq_cols], as_numeric_zero)

rep <- rep_raw |>
  group_by(CDR3) |>
  summarise(across(all_of(freq_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop")

write_csv_clean(rep, file.path(processed_dir, "repertoire_matrix_clean.csv"))
saveRDS(rep, file.path(processed_dir, "repertoire_matrix_clean.rds"))

message("Unique CDR3 sequences in cleaned matrix: ", nrow(rep))

# ==============================================================================
# SECTION 4. DEFINE DRTCC AND SRTCC CLONOTYPE SETS
# Figures: Fig. 2a-e and Fig. 2h
# ==============================================================================
# Donor-reactive T-cell clonotypes (DRTCCs):
#   - Pre-transplant DRTCCs are imported from the curated pre-transplant lists.
#   - Post-transplant DRTCCs require stimulated productive frequency >=1e-5
#     and >2-fold enrichment over the designated unstimulated blood repertoire.
#
# Spike-reactive T-cell clonotypes (SRTCCs):
#   - Timepoint-specific assignments are imported from the curated spike-
#     reactive workbook generated from the spike-reactive culture analysis.

pre_cd4 <- read_cdr3_csv(pre_cd4_drtcc_file)
pre_cd8 <- read_cdr3_csv(pre_cd8_drtcc_file)

drtcc <- list(
  CD4 = list(
    PreTx = pre_cd4,
    `30M` = make_post_drtcc(
      rep, "30M_CFSE_Low_CD4_TCRB", "30M_Unstim_CD4_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    ),
    `33.7M` = make_post_drtcc(
      rep, "33-7M_CFSE_low_CD4_TCRB", "30M_Unstim_CD4_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    ),
    `54M` = make_post_drtcc(
      rep, "54M_CFSE_Low_CD4_TCRB", "54M_Unstim_CD4_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    )
  ),
  CD8 = list(
    PreTx = pre_cd8,
    `30M` = make_post_drtcc(
      rep, "30M_CFSE_Low_CD8_TCRB", "30M_Unstim_CD8_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    ),
    `33.7M` = make_post_drtcc(
      rep, "33-7M_CFSE_low_CD8_TCRB", "30M_Unstim_CD8_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    ),
    `54M` = make_post_drtcc(
      rep, "54M_CFSE_Low_CD8_TCRB", "54M_Unstim_CD8_TCRB",
      min_functional_frequency, fold_enrichment_cutoff
    )
  )
)

if (!file.exists(spike_reactive_file)) stop("Input not found: ", spike_reactive_file)
spike_cd4 <- read_excel(spike_reactive_file, sheet = "CD4_spike_reactive")
spike_cd8 <- read_excel(spike_reactive_file, sheet = "CD8_spike_reactive")

assert_columns(
  spike_cd4,
  c("CDR3", "CD4_hit_30M_spike", "CD4_hit_33_7M_spike", "CD4_hit_54M_spike"),
  "CD4 spike workbook"
)
assert_columns(
  spike_cd8,
  c("CDR3", "CD8_hit_30M_spike", "CD8_hit_33_7M_spike", "CD8_hit_54M_spike"),
  "CD8 spike workbook"
)

spike_cd4$CDR3 <- clean_cdr3(spike_cd4$CDR3)
spike_cd8$CDR3 <- clean_cdr3(spike_cd8$CDR3)

get_spike_set <- function(df, hit_col) {
  unique(df$CDR3[!is.na(df$CDR3) & truthy(df[[hit_col]])])
}

srtcc <- list(
  CD4 = list(
    `30M` = get_spike_set(spike_cd4, "CD4_hit_30M_spike"),
    `33.7M` = get_spike_set(spike_cd4, "CD4_hit_33_7M_spike"),
    `54M` = get_spike_set(spike_cd4, "CD4_hit_54M_spike")
  ),
  CD8 = list(
    `30M` = get_spike_set(spike_cd8, "CD8_hit_30M_spike"),
    `33.7M` = get_spike_set(spike_cd8, "CD8_hit_33_7M_spike"),
    `54M` = get_spike_set(spike_cd8, "CD8_hit_54M_spike")
  )
)

saveRDS(
  list(DRTCC = drtcc, SRTCC = srtcc),
  file.path(processed_dir, "functional_clone_sets.rds")
)

functional_counts <- data.frame(
  Lineage = rep(c("CD4", "CD8"), each = 7),
  Population = rep(
    c("DRTCC_PreTx", "DRTCC_30M", "DRTCC_33.7M", "DRTCC_54M",
      "SRTCC_30M", "SRTCC_33.7M", "SRTCC_54M"),
    2
  ),
  N = c(
    length(drtcc$CD4$PreTx), length(drtcc$CD4$`30M`),
    length(drtcc$CD4$`33.7M`), length(drtcc$CD4$`54M`),
    length(srtcc$CD4$`30M`), length(srtcc$CD4$`33.7M`),
    length(srtcc$CD4$`54M`),
    length(drtcc$CD8$PreTx), length(drtcc$CD8$`30M`),
    length(drtcc$CD8$`33.7M`), length(drtcc$CD8$`54M`),
    length(srtcc$CD8$`30M`), length(srtcc$CD8$`33.7M`),
    length(srtcc$CD8$`54M`)
  )
)
write_csv_clean(functional_counts, file.path(tables_dir, "functional_clone_counts.csv"))

# ==============================================================================
# SECTION 5. 12 MLR x SPIKE PAIRWISE OVERLAP COMPARISONS
# Figures: Fig. 2b (CD4) and Fig. 2d (CD8)
# ==============================================================================
# Four donor-reactive MLR repertoires (PreTx, 30M, 33.7M, 54M) are compared
# with three spike-reactive repertoires (30M, 33.7M, 54M), producing 12
# pairwise comparisons per lineage. Exact shared CDR3 sequences are counted.

mlr_freq_cols <- list(
  PreTx = c(CD4 = "UPITT_1006_CFSElow_CD4_TCRB", CD8 = "UPITT_1006_CFSElow_CD8_TCRB"),
  `30M` = c(CD4 = "30M_CFSE_Low_CD4_TCRB", CD8 = "30M_CFSE_Low_CD8_TCRB"),
  `33.7M` = c(CD4 = "33-7M_CFSE_low_CD4_TCRB", CD8 = "33-7M_CFSE_low_CD8_TCRB"),
  `54M` = c(CD4 = "54M_CFSE_Low_CD4_TCRB", CD8 = "54M_CFSE_Low_CD8_TCRB")
)

spike_freq_cols <- list(
  `30M` = c(CD4 = "30M_SCoV2_S_CD4_TCRB", CD8 = "30M_SCoV2_S_CD8_TCRB"),
  `33.7M` = c(CD4 = "33-7M_SCoV2_S_CD4_TCRB", CD8 = "33-7M_SCoV2_S_CD8_TCRB"),
  `54M` = c(CD4 = "54M_SCoV2_S_CD4_TCRB", CD8 = "54M_SCoV2_S_CD8_TCRB")
)

run_pair <- function(lineage, mlr_tp, spike_tp) {
  mlr_set <- drtcc[[lineage]][[mlr_tp]]
  spike_set <- srtcc[[lineage]][[spike_tp]]
  overlap <- intersect(mlr_set, spike_set)

  mlr_col <- unname(mlr_freq_cols[[mlr_tp]][lineage])
  spike_col <- unname(spike_freq_cols[[spike_tp]][lineage])

  clone_table <- data.frame(CDR3 = overlap) |>
    left_join(rep |> select(CDR3, all_of(c(mlr_col, spike_col))), by = "CDR3") |>
    transmute(
      CDR3,
      MLR_productive_frequency = .data[[mlr_col]],
      MLR_productive_frequency_percent = 100 * .data[[mlr_col]],
      Spike_productive_frequency = .data[[spike_col]],
      Spike_productive_frequency_percent = 100 * .data[[spike_col]]
    ) |>
    arrange(desc(MLR_productive_frequency), desc(Spike_productive_frequency))

  summary <- data.frame(
    Lineage = lineage,
    MLR_timepoint = mlr_tp,
    Spike_timepoint = spike_tp,
    MLR_DRTCC_count = length(mlr_set),
    Spike_SRTCC_count = length(spike_set),
    Shared_count = length(overlap),
    Percent_MLR_shared = ifelse(length(mlr_set) > 0,
                                100 * length(overlap) / length(mlr_set), NA_real_),
    Percent_Spike_shared = ifelse(length(spike_set) > 0,
                                  100 * length(overlap) / length(spike_set), NA_real_),
    Cumulative_MLR_frequency_percent = sum(clone_table$MLR_productive_frequency_percent, na.rm = TRUE),
    Cumulative_Spike_frequency_percent = sum(clone_table$Spike_productive_frequency_percent, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(summary = summary, clones = clone_table)
}

pair_results <- list()
for (lineage in c("CD4", "CD8")) {
  for (m in c("PreTx", "30M", "33.7M", "54M")) {
    for (s in c("30M", "33.7M", "54M")) {
      key <- paste(lineage, m, s, sep = "__")
      pair_results[[key]] <- run_pair(lineage, m, s)
    }
  }
}

pair_summary <- bind_rows(lapply(pair_results, `[[`, "summary"))
write_csv_clean(pair_summary, file.path(tables_dir, "Fig2BD_pairwise_overlap_summary.csv"))

# Union of all CD8 donor-spike shared clonotypes across the 12 comparisons.
cd8_shared <- sort(unique(unlist(lapply(
  pair_results[grepl("^CD8__", names(pair_results))],
  function(x) x$clones$CDR3
))))

pairwise_recurrence <- bind_rows(lapply(names(pair_results), function(k) {
  if (!grepl("^CD8__", k)) return(NULL)
  x <- pair_results[[k]]$clones
  if (nrow(x) == 0) return(NULL)
  data.frame(CDR3 = x$CDR3, Comparison = k)
})) |>
  distinct() |>
  count(CDR3, name = "Pairwise_overlap_count")

cd8_shared_catalog <- data.frame(CDR3 = cd8_shared) |>
  left_join(pairwise_recurrence, by = "CDR3") |>
  arrange(desc(Pairwise_overlap_count), CDR3)

write_csv_clean(
  cd8_shared_catalog,
  file.path(processed_dir, "CD8_shared_clone_catalog.csv")
)
message("Unique CD8 donor-spike shared clonotypes: ", nrow(cd8_shared_catalog))

# ==============================================================================
# SECTION 6. CUMULATIVE PRODUCTIVE FREQUENCY OF SHARED CLONOTYPES
# Figures: Fig. 2c (CD4) and Fig. 2e (CD8)
# ==============================================================================
# For each of the 12 MLR x spike comparisons, the productive frequencies of the
# exact shared-clonotype set are summed separately in the donor-reactive MLR and
# spike-reactive repertoires. Values are reported as percentages.

fig2bd_counts <- pair_summary |>
  select(Lineage, MLR_timepoint, Spike_timepoint, Shared_count,
         MLR_DRTCC_count, Spike_SRTCC_count)
write_csv_clean(fig2bd_counts, file.path(tables_dir, "Fig2BD_overlap_counts.csv"))

fig2ce_frequency <- bind_rows(
  pair_summary |>
    transmute(
      Lineage, MLR_timepoint, Spike_timepoint,
      Repertoire = "MLR",
      Cumulative_productive_frequency_percent = Cumulative_MLR_frequency_percent
    ),
  pair_summary |>
    transmute(
      Lineage, MLR_timepoint, Spike_timepoint,
      Repertoire = "Spike",
      Cumulative_productive_frequency_percent = Cumulative_Spike_frequency_percent
    )
)
write_csv_clean(
  fig2ce_frequency,
  file.path(tables_dir, "Fig2CE_cumulative_shared_frequency.csv")
)

# Optional reproducibility plot. Final publication formatting may be performed
# separately using the exported source-data table.
p_ce <- ggplot(
  fig2ce_frequency,
  aes(
    x = MLR_timepoint,
    y = Cumulative_productive_frequency_percent,
    fill = Repertoire
  )
) +
  geom_col(
    position = position_dodge(width = 0.75), width = 0.65,
    colour = "black", linewidth = 0.2
  ) +
  facet_grid(Lineage ~ Spike_timepoint, scales = "free_y") +
  labs(x = NULL, y = "Cumulative productive frequency (%)") +
  theme_classic(base_size = 10)

ggsave(
  file.path(figures_dir, "Fig2CE_cumulative_shared_frequency.pdf"),
  p_ce, width = 9, height = 4.5, units = "in"
)

# ==============================================================================
# SECTION 7. DOMINANT CIRCULATING CD8 CLONOTYPES / ALLUVIAL ANALYSIS
# Figure: Fig. 2f
# ==============================================================================
# Selection is deliberately independent of donor- or spike-reactive assignment.
# Step 1: identify the top 20 circulating CD8 clonotypes at each of five blood
#         time points (PreTx, 6M, 12M, 30M, 54M).
# Step 2: take the union of those five top-20 lists.
# Step 3: calculate mean productive frequency across all five blood samples,
#         assigning 0 to undetected clones because the matrix stores absence as 0.
# Step 4: retain clonotypes with mean productive frequency >=0.75%.
# In the finalized analysis this yields 19 clonotypes. Three donor-spike shared
# clonotypes among these 19 are highlighted after abundance-based selection.

blood_cols <- c(
  PreTx = "UPITT_1006_UNSTIM_CD8_TCRB",
  `6M` = "UPITT_1006_6mth_CD8_TCRB",
  `12M` = "UPITT_1006_12mth_CD8_TCRB",
  `30M` = "30M_Unstim_CD8_TCRB",
  `54M` = "54M_Unstim_CD8_TCRB"
)
assert_columns(rep, c("CDR3", unname(blood_cols)), "repertoire matrix")

top_union <- unique(unlist(lapply(blood_cols, function(col) {
  rep |>
    arrange(desc(.data[[col]])) |>
    slice_head(n = alluvial_top_n_each_timepoint) |>
    pull(CDR3)
})))

message("Unique clonotypes in union of five top-20 lists: ", length(top_union))

fig2f_wide <- rep |>
  filter(CDR3 %in% top_union) |>
  transmute(
    CDR3,
    PreTx = .data[[blood_cols["PreTx"]]],
    `6M` = .data[[blood_cols["6M"]]],
    `12M` = .data[[blood_cols["12M"]]],
    `30M` = .data[[blood_cols["30M"]]],
    `54M` = .data[[blood_cols["54M"]]]
  ) |>
  mutate(
    Mean_percent = 100 * rowMeans(across(c(PreTx, `6M`, `12M`, `30M`, `54M`))),
    Max_frequency = pmax(PreTx, `6M`, `12M`, `30M`, `54M`),
    Shared_DRTCC_SRTCC = CDR3 %in% cd8_shared
  ) |>
  filter(Mean_percent >= alluvial_mean_percent_cutoff) |>
  arrange(desc(Max_frequency))

message("Clonotypes retained for Fig. 2f: ", nrow(fig2f_wide))
if (nrow(fig2f_wide) != 19) {
  warning(
    "The finalized manuscript analysis contains 19 clonotypes; this run produced ",
    nrow(fig2f_wide), ". Confirm input files and threshold definitions."
  )
}

clone_order <- fig2f_wide$CDR3
fig2f_long <- fig2f_wide |>
  select(CDR3, PreTx, `6M`, `12M`, `30M`, `54M`) |>
  pivot_longer(-CDR3, names_to = "Timepoint", values_to = "Frequency") |>
  mutate(
    Timepoint = factor(Timepoint, levels = c("PreTx", "6M", "12M", "30M", "54M")),
    CDR3 = factor(CDR3, levels = rev(clone_order)),
    Productive_frequency_percent = 100 * Frequency
  )

candidate_colors <- unique(c(
  brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"),
  brewer.pal(8, "Accent"), brewer.pal(9, "Set1")
))
candidate_colors <- candidate_colors[
  !toupper(candidate_colors) %in% toupper(unname(special_clone_colors))
]
ordinary_clones <- setdiff(clone_order, names(special_clone_colors))
clone_colors <- setNames(rep(NA_character_, length(clone_order)), clone_order)
clone_colors[ordinary_clones] <- candidate_colors[seq_along(ordinary_clones)]
for (nm in intersect(names(special_clone_colors), clone_order)) {
  clone_colors[nm] <- special_clone_colors[nm]
}

p_f <- ggplot(
  fig2f_long,
  aes(
    x = Timepoint, stratum = CDR3, alluvium = CDR3,
    y = Frequency, fill = CDR3
  )
) +
  geom_alluvium(
    width = 0.25, knot.pos = 0.5, curve_type = "xspline",
    colour = "black", linewidth = 0.25
  ) +
  geom_stratum(
    width = 0.25, colour = "black", linewidth = 0.25, alpha = 0.35
  ) +
  scale_fill_manual(values = clone_colors, guide = "none") +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Cumulative productive frequency (%)") +
  theme_classic(base_size = 11)

ggsave(
  file.path(figures_dir, "Fig2F_dominant_CD8_alluvial.pdf"),
  p_f, width = 8.5, height = 4.6, units = "in"
)
ggsave(
  file.path(figures_dir, "Fig2F_dominant_CD8_alluvial.png"),
  p_f, width = 8.5, height = 4.6, units = "in", dpi = 600
)

fig2f_wide |>
  mutate(Plot_color = clone_colors[CDR3]) |>
  write_csv_clean(file.path(tables_dir, "Fig2F_alluvial_source_data.csv"))

# ==============================================================================
# SECTION 8. BIOPSY REPRESENTATION BY CLONOTYPE CLASS
# Figure: Fig. 2h
# ==============================================================================
# This analysis compares representation of predefined CD8 clonotype classes in
# the rejecting for-cause biopsy (FcB):
#   - Pre-DRTCC only
#   - Post-DRTCC only
#   - SRTCC only
#   - Pre-CRTCC
#   - Post-CRTCC
#   - nDRTCC-nSpike comparator clones
#
# CRTCC classes are mutually exclusive. If a shared clonotype has both pre- and
# post-transplant donor-reactive assignments, it is classified as Pre-CRTCC.
# Post-CRTCCs therefore lack a pre-transplant donor-reactive assignment.

pre_drtcc <- drtcc$CD8$PreTx
post_drtcc <- unique(c(drtcc$CD8$`30M`, drtcc$CD8$`33.7M`, drtcc$CD8$`54M`))
all_srtcc <- unique(c(srtcc$CD8$`30M`, srtcc$CD8$`33.7M`, srtcc$CD8$`54M`))
all_drtcc <- union(pre_drtcc, post_drtcc)

pre_crtcc <- intersect(pre_drtcc, all_srtcc)
post_crtcc <- setdiff(intersect(post_drtcc, all_srtcc), pre_crtcc)

pre_drtcc_only <- setdiff(pre_drtcc, all_srtcc)
post_drtcc_only <- setdiff(post_drtcc, all_srtcc)
srtcc_only <- setdiff(all_srtcc, all_drtcc)

assert_columns(
  rep,
  c("UPITT_1006_12mth_CD8_TCRB", "30M_Unstim_CD8_TCRB"),
  "repertoire matrix"
)

rep_comparator <- rep |>
  mutate(
    Blood12 = UPITT_1006_12mth_CD8_TCRB,
    Blood30 = `30M_Unstim_CD8_TCRB`,
    Blood30_templates = round(Blood30 * cd8_30m_blood_total_templates),
    Fold_12_to_30 = case_when(
      Blood12 > 0 ~ Blood30 / Blood12,
      Blood12 == 0 & Blood30 > 0 ~ Inf,
      TRUE ~ 0
    )
  )

comparator <- rep_comparator |>
  filter(
    Blood30 >= comparator_min_frequency,
    Blood30_templates >= comparator_min_templates,
    Fold_12_to_30 >= comparator_fold_expansion,
    !CDR3 %in% union(all_drtcc, all_srtcc)
  ) |>
  pull(CDR3)

fig2h_categories <- list(
  `Pre-DRTCC only` = pre_drtcc_only,
  `Post-DRTCC only` = post_drtcc_only,
  `SRTCC only` = srtcc_only,
  `Pre-CRTCC` = pre_crtcc,
  `Post-CRTCC` = post_crtcc,
  `nDRTCC-nSpike` = comparator
)

fcb_col <- "Upitt_1006_FcB9M_post_wean_TCRB"
assert_columns(rep, c("CDR3", fcb_col), "repertoire matrix")

fig2h_summary <- bind_rows(lapply(names(fig2h_categories), function(cat) {
  ids <- unique(fig2h_categories[[cat]])

  d <- rep |>
    filter(CDR3 %in% ids) |>
    mutate(
      FcB_frequency = .data[[fcb_col]],
      FcB_templates = round(FcB_frequency * cd8_fcb_total_templates),
      Detected = FcB_frequency > 0
    )

  n_total <- length(ids)
  n_detected <- sum(d$Detected)
  total_templates <- sum(d$FcB_templates, na.rm = TRUE)

  data.frame(
    Category = cat,
    Total_clonotypes = n_total,
    Detected_in_FcB = n_detected,
    Clonal_detection_rate_percent = ifelse(
      n_total > 0, 100 * n_detected / n_total, NA_real_
    ),
    Cumulative_productive_frequency_in_FcB_percent =
      100 * sum(d$FcB_frequency, na.rm = TRUE),
    Total_FcB_templates = total_templates,
    Biopsy_templates_per_detected_clonotype = ifelse(
      n_detected > 0, total_templates / n_detected, NA_real_
    )
  )
}))

write_csv_clean(fig2h_summary, file.path(tables_dir, "Fig2H_biopsy_enrichment_summary.csv"))

fig2h_membership <- bind_rows(lapply(names(fig2h_categories), function(cat) {
  data.frame(Category = cat, CDR3 = unique(fig2h_categories[[cat]]))
})) |>
  left_join(
    rep |>
      transmute(
        CDR3,
        FcB_productive_frequency_percent = 100 * .data[[fcb_col]],
        FcB_templates = round(.data[[fcb_col]] * cd8_fcb_total_templates)
      ),
    by = "CDR3"
  ) |>
  mutate(Detected_in_FcB = FcB_productive_frequency_percent > 0)

write_csv_clean(fig2h_membership, file.path(tables_dir, "Fig2H_clone_membership.csv"))

# ==============================================================================
# SECTION 9. PRODUCTIVE FREQUENCY AND ABUNDANCE-RANK ANALYSIS
# Figures: Fig. 2i,j
# ==============================================================================
# Five representative CRTCCs are tracked across longitudinal blood, donor-
# reactive MLRs, spike-reactive cultures and liver biopsies. Fig. 2i uses
# productive frequency; Fig. 2j uses within-repertoire abundance rank.
# Competition ranking is used: ties share the same rank, and undetected clones
# are reported as ND rather than being assigned an artificial terminal rank.

sample_cols <- c(
  `Blood 0M` = "UPITT_1006_UNSTIM_CD8_TCRB",
  `Blood 6M` = "UPITT_1006_6mth_CD8_TCRB",
  `Blood 12M` = "UPITT_1006_12mth_CD8_TCRB",
  `Blood 30M` = "30M_Unstim_CD8_TCRB",
  `Blood 54M` = "54M_Unstim_CD8_TCRB",
  `MLR 0M` = "UPITT_1006_CFSElow_CD8_TCRB",
  `MLR 30M` = "30M_CFSE_Low_CD8_TCRB",
  `MLR 33.7M` = "33-7M_CFSE_low_CD8_TCRB",
  `MLR 54M` = "54M_CFSE_Low_CD8_TCRB",
  `Spike 30M` = "30M_SCoV2_S_CD8_TCRB",
  `Spike 33.7M` = "33-7M_SCoV2_S_CD8_TCRB",
  `Spike 54M` = "54M_SCoV2_S_CD8_TCRB",
  `Liver 12M` = "UPITT_1006_12mth_LiverBx_TCRB",
  `Liver 33.7M` = "Upitt_1006_FcB9M_post_wean_TCRB"
)
assert_columns(rep, c("CDR3", unname(sample_cols)), "repertoire matrix")

rank_table <- rep |> select(CDR3)
for (nm in names(sample_cols)) {
  col <- unname(sample_cols[nm])
  rank_table[[nm]] <- competition_rank(rep[[col]])
}

freq_table <- rep |>
  filter(CDR3 %in% representative_crtccs) |>
  select(CDR3, all_of(unname(sample_cols)))
names(freq_table)[-1] <- names(sample_cols)
freq_table[-1] <- lapply(freq_table[-1], function(x) 100 * x)

rank_display <- rank_table |>
  filter(CDR3 %in% representative_crtccs)

freq_table$CDR3 <- factor(freq_table$CDR3, levels = representative_crtccs)
rank_display$CDR3 <- factor(rank_display$CDR3, levels = representative_crtccs)
freq_table <- arrange(freq_table, CDR3)
rank_display <- arrange(rank_display, CDR3)
freq_table$CDR3 <- as.character(freq_table$CDR3)
rank_display$CDR3 <- as.character(rank_display$CDR3)

fig2ij_long <- freq_table |>
  pivot_longer(-CDR3, names_to = "Sample", values_to = "Productive_frequency_percent") |>
  left_join(
    rank_display |>
      pivot_longer(-CDR3, names_to = "Sample", values_to = "Abundance_rank"),
    by = c("CDR3", "Sample")
  ) |>
  mutate(
    Rank_display = ifelse(is.na(Abundance_rank), "ND", as.character(Abundance_rank))
  )

write_csv_clean(
  fig2ij_long,
  file.path(tables_dir, "Fig2IJ_representative_CRTCC_frequency_rank.csv")
)

wb_ij <- createWorkbook()
addWorksheet(wb_ij, "Productive_frequency_percent")
writeData(wb_ij, "Productive_frequency_percent", freq_table)
addWorksheet(wb_ij, "Abundance_rank")
writeData(wb_ij, "Abundance_rank", rank_display)
addWorksheet(wb_ij, "Long_source_data")
writeData(wb_ij, "Long_source_data", fig2ij_long)
addWorksheet(wb_ij, "Notes")
writeData(
  wb_ij,
  "Notes",
  data.frame(Note = c(
    "Rank 1 denotes the most abundant detected clonotype within a repertoire.",
    "Competition ranking is used: identical productive frequencies share the same rank.",
    "Clonotypes with zero productive frequency are reported as ND and are not ranked."
  ))
)
saveWorkbook(
  wb_ij,
  file.path(tables_dir, "Fig2IJ_representative_CRTCC_frequency_rank.xlsx"),
  overwrite = TRUE
)

# ==============================================================================
# SECTION 10. PACKAGE AND VERSION INFORMATION
# Code availability / reproducibility
# ==============================================================================
# This writes the exact R and package versions from the environment in which
# the script is executed. Include the resulting file with the archived code.

session_info_file <- file.path(results_dir, "sessionInfo.txt")
sink(session_info_file)
cat("Pt1006 TCR repertoire analysis - R session information\n")
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
print(sessionInfo())
sink()
message("Wrote package/version information: ", session_info_file)

# ==============================================================================
# END OF ANALYSIS
# ==============================================================================
