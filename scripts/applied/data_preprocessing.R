# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
source("scripts/utils.R")

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "MatchIt", "fastDummies"
)
new_packages <- list_of_packages[
  !(list_of_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(
    new_packages,
    repos = "http://cran.us.r-project.org"
  )
}
invisible(lapply(list_of_packages, library, character.only = TRUE))

# Set working directory and paths
setwd("~/projects/alphafold-impact/")
figures <- "data/05_model_output/figures/"
tables <- "data/05_model_output/tables/"
figs <- list()

# Assign commonly used dplyr functions
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

# Read and clean the dataset
sb_data_qtly <- read_parquet(
  "data/04_outputs/applied_labs/staggered/outputs_quarterly.parquet"
)

# Clean column names
colnames(sb_data_qtly) <- gsub(",", "", colnames(sb_data_qtly))


# Complete missing rows and fill NA values
sb_data_qtly <- sb_data_qtly %>%
  # group_by(pi_id, parent_time) %>%
  # complete(time = 1:25) %>%
  # ungroup() %>%
  group_by(pi_id) %>%
  fill(
    seed, strong, institution_country_code,
    starts_with("institution_"),
    .direction = "downup"
  ) %>%
  replace_na(list(
    num_publications = 0, cited_by_count = 0, ct0 = 0, ct1 = 0,
    ca_count = 0, patent_count = 0, pdb_share = 0, # resolution = 0, R_free = 0,
    covid_share_2020 = 0,
    institution_country_code = "OTH", institution_type = "other"
  )) %>%
  ungroup()

# Factorize relevant columns
sb_data_qtly <- sb_data_qtly %>%
  mutate(
    strong = as.factor(strong),
    time_qtly = as.factor(time),
    country = as.factor(institution_country_code),
    institution_type = as.factor(institution_type)
  )

# Replace NA in columns starting with "mesh_" with 0
sb_data_qtly <- sb_data_qtly %>%
  mutate_at(vars(starts_with("mesh_")), ~ replace_na(., 0))

# Rename and create new variables
sb_data_qtly <- sb_data_qtly %>%
  rename(g = parent_time) %>%
  mutate(
    g = ifelse(seed == "other", NA, g),
    g_af = ifelse(str_detect(seed, "af"), g, NA),
    g_ct = ifelse(str_detect(seed, "ct"), g, NA),
    g_ct_ai = ifelse(str_detect(seed, "ct_ai"), g, NA),
    g_ct_noai = ifelse(str_detect(seed, "ct_noai"), g, NA),
    treatment_af_dyn = ifelse(time > g_af, 1, 0),
    treatment_af = ifelse(time > 15 & !is.na(g_af), 1, 0),
    treatment_ct = ifelse(time > g_ct, 1, 0),
    experimental_share = experimental / num_publications,
    protein_share = protein_concept / num_publications
  ) %>%
  replace_na(list(
    treatment_af_dyn = 0, treatment_af = 0, treatment_ct = 0,
    treatment_ct_ai = 0, treatment_ct_noai = 0, pdb_share = 0,
    protein_share = 0, experimental_share = 0
  )) %>%
  select(c(g, g_af, g_ct, g_ct_ai, g_ct_noai, time), everything()) %>%
  arrange(pi_id, time)

# Backfill pdb_share, resolution, R_free
sb_data_qtly <- sb_data_qtly %>%
  select(-institution_country_code) %>%
  group_by(pi_id) %>%
  fill(pdb_share, resolution, R_free) %>%
  ungroup()

# Log-transform some columns
sb_data_qtly <- sb_data_qtly %>%
  mutate(
    num_publications = log1p(num_publications),
    ca_count = log1p(ca_count),
    patent_count = log1p(patent_count),
    ct0 = log1p(ct0),
    ct1 = log1p(ct1),
    cited_by_count = log1p(cited_by_count)
  )

# Standardise publication-related counts
sb_data_qtly <- sb_data_qtly %>%
  group_by(time) %>%
  mutate(
    num_publications_std = scale(num_publications, center = TRUE, scale = TRUE),
    patent_count_std = scale(patent_count, center = TRUE, scale = TRUE),
    ca_count_std = scale(ca_count, center = TRUE, scale = TRUE),
    cited_by_count_std = scale(cited_by_count, center = TRUE, scale = TRUE)
  ) %>%
  ungroup()

sb_data_qtly <- sb_data_qtly %>%
  mutate_at(vars(starts_with("field_"), starts_with("mesh_")), ~ replace_na(., 0))
# ------------------------------------------------------------------------------
# CREATING DUMMIES FOR EVENT STUDY
# ------------------------------------------------------------------------------

sb_data_qtly <- sb_data_qtly %>%
  mutate(
    rel_treat_af = time - g_af,
    rel_treat_af_strong = ifelse(strong == "TRUE", time - g_af, NA),
    rel_treat_ct = time - g_ct,
    rel_treat_ct_ai = time - g_ct_ai,
    rel_treat_ct_noai = time - g_ct_noai
  )

sb_data_qtly <- sb_data_qtly %>%
  dummy_cols(
    select_columns = c(
      "rel_treat_af", "rel_treat_af_strong", "rel_treat_ct",
      "rel_treat_ct_ai", "rel_treat_ct_noai"
    ),
    remove_first_dummy = FALSE,
    remove_selected_columns = TRUE
  )

# Fill NAs in dummy columns with 0
sb_data_qtly <- sb_data_qtly %>%
  mutate_at(vars(starts_with("rel_treat_")), as.numeric) %>%
  mutate_at(vars(starts_with("rel_treat_")), ~ replace_na(., 0))

# ------------------------------------------------------------------------------
# SUBSET CREATION
# ------------------------------------------------------------------------------

# PDB subset: high-pdb researchers
sb_data_qtly_2018_2020 <- sb_data_qtly %>% filter(time %in% 1:13)
avg_pdb_share <- sb_data_qtly_2018_2020 %>%
  group_by(pi_id) %>%
  summarise(avg_pdb_share = mean(pdb_share, na.rm = TRUE))
percentile_75_pdb <- quantile(avg_pdb_share$avg_pdb_share, 0.75, na.rm = TRUE)
high_pdb <- avg_pdb_share %>%
  filter(avg_pdb_share > percentile_75_pdb) %>%
  pull(pi_id) # nolint
sb_data_qtly$high_pdb <- ifelse(sb_data_qtly$pi_id %in% high_pdb, 1, 0)

### SUBSETS ###
sub_samples <- list(
  "all" = sb_data_qtly,
  "af_ct" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | str_detect(seed, "ct")),
  "af_ct_ai" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | str_detect(seed, "ct_ai")),
  "af_ct_noai" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | str_detect(seed, "ct_noai")),
  "af_ct_w_high_pdb" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | str_detect(seed, "ct")) %>%
    semi_join(
      avg_pdb_share %>%
        filter(avg_pdb_share > percentile_75_pdb),
      by = "pi_id"
    )
)

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching)
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
cols <- c(
  "num_publications", "cited_by_count", "ct0", "ct1", "pdb_share",
  "protein_share", "experimental_share", "patent_count",
  "field_biochemistry_genetics_and_molecular_biology",
  "mesh_C"
)

# Filter and prepare data for collapsing
sb_data_qtly_y0_cem <- sb_data_qtly %>%
  filter(time %in% 1:5) %>%
  mutate(treatment_af_ct = ifelse(!is.na(g_af) & is.na(g_ct), 1, 0)) %>%
  filter(complete.cases(field_biochemistry_genetics_and_molecular_biology, mesh_C)) # nolint

# Collapse the data to ensure one observation per pi_id
sb_data_qtly_y0_cem_collapsed <- sb_data_qtly_y0_cem %>%
  group_by(pi_id) %>%
  summarize(
    treatment_af_ct = max(treatment_af_ct),  # Ensure treatment remains binary
    across(everything(), \(x) mean(x, na.rm = TRUE))
  )

# CEM matching on the collapsed data
match_out_af <- matchit(
  as.formula(paste0("treatment_af_ct ~ ", paste0(cols, collapse = " + "))),
  data = sb_data_qtly_y0_cem_collapsed, method = "cem", k2k = TRUE
)

# Store matched data
cem_data <- match.data(match_out_af)

# Sample based on the matched group
sb_data_qtly_cem <- sb_data_qtly %>% semi_join(cem_data, by = "pi_id")
sub_samples[["af_ct_cem"]] <- sb_data_qtly_cem

# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, "data/05_model_output/applied_sub_samples.rds")
