# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 100)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "aws.s3", "yaml", "zoo", "MatchIt"
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

# Set working directory
setwd("~/projects/alphafold-impact/")
pathdir <- "data/05_model_output/authors/ecr/quarterly/data/"

# Create directories if they do not exist
if (!dir.exists(pathdir)) {
  dir.create(pathdir, recursive = TRUE)
}

# Assign commonly used dplyr functions
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

credentials <- suppressWarnings(yaml.load_file("conf/base/credentials.yml"))

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key, # nolint
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret, # nolint
  "AWS_DEFAULT_REGION" = "eu-west-2" # nolint
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
path <- "oct/03_primary/ecr/publications_quarterly.parquet" # nolint

# Fetch the data from the S3 bucket
ecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# ------------------------------------------------------------------------------
# Strong Data Prep
# ------------------------------------------------------------------------------

# Calculate intent ratios for each author across all periods
author_intent_ratios <- ecr_data %>%
  group_by(author) %>%
  summarise(
    af_intent_ratio = sum(af_with_intent, na.rm = TRUE) /
      (sum(af_with_intent, na.rm = TRUE) + sum(af_unknown, na.rm = TRUE)),
    ct_ai_intent_ratio = sum(ct_ai_with_intent, na.rm = TRUE) /
      (sum(ct_ai_with_intent, na.rm = TRUE) + sum(ct_ai_unknown, na.rm = TRUE)),
    ct_noai_intent_ratio = sum(ct_noai_with_intent, na.rm = TRUE) /
      (sum(ct_noai_with_intent, na.rm = TRUE) + sum(ct_noai_unknown, na.rm = TRUE)) # nolint
  )

# Join the ratios back to main data and create the strong0/strong1 columns
ecr_data <- ecr_data %>%
  left_join(author_intent_ratios, by = "author") %>%
  mutate(
    # Calculate row mean of intent ratios, handling NaN cases properly
    all_intents_high = {
      means <- rowMeans(
        cbind(
          af_intent_ratio,
          ct_ai_intent_ratio,
          ct_noai_intent_ratio
        ),
        na.rm = TRUE
      )
      # If all values were NA (resulting in NaN), or mean <= 0.5, return FALSE
      !is.nan(means) & means > 0.5
    },

    # Only pass values if intent ratio mean is high enough
    "af_strong0" = if_else(all_intents_high, af_weak, NA_real_),
    "af_strong1" = if_else(all_intents_high, af_strong, NA_real_),
    "ct_ai_strong0" = if_else(all_intents_high, ct_ai_weak, NA_real_),
    "ct_ai_strong1" = if_else(all_intents_high, ct_ai_strong, NA_real_),
    "ct_noai_strong0" = if_else(all_intents_high, ct_noai_weak, NA_real_),
    "ct_noai_strong1" = if_else(all_intents_high, ct_noai_strong, NA_real_)
  )

# If strong1 exists for any type, set corresponding strong0 to 0
ecr_data <- ecr_data %>%
  group_by(author) %>%
  mutate(
    "af_strong0" = case_when(
      is.na(`af_strong0`) ~ NA_real_,
      any(`af_strong1` > 0, na.rm = TRUE) ~ 0,
      TRUE ~ `af_strong0`
    ),
    "ct_ai_strong0" = case_when(
      is.na(`ct_ai_strong0`) ~ NA_real_,
      any(`ct_ai_strong1` > 0, na.rm = TRUE) ~ 0,
      TRUE ~ `ct_ai_strong0`
    ),
    "ct_noai_strong0" = case_when(
      is.na(`ct_noai_strong0`) ~ NA_real_,
      any(`ct_noai_strong1` > 0, na.rm = TRUE) ~ 0,
      TRUE ~ `ct_noai_strong0`
    )
  ) %>%
  ungroup()

# create interactions between the strong variables
ecr_data <- ecr_data %>%
  mutate(
    "af_ct_ai_strong0" = `af_strong0` * `ct_ai_strong0`,
    "af_ct_ai_strong1" = `af_strong1` * `ct_ai_strong1`,
    "af_ct_noai_strong0" = `af_strong0` * `ct_noai_strong0`,
    "af_ct_noai_strong1" = `af_strong1` * `ct_noai_strong1`,
    "ct_ai_ct_noai_strong0" = `ct_ai_strong0` * `ct_noai_strong0`,
    "ct_ai_ct_noai_strong1" = `ct_ai_strong1` * `ct_noai_strong1`
  )

# ------------------------------------------------------------------------------
# Create Extensive Margin Variables
# ------------------------------------------------------------------------------

ecr_data <- ecr_data %>%
  mutate(
    # Basic variables - convert NAs to 0 and then to binary
    af = ifelse(af > 0, 1, 0),
    ct_ai = ifelse(ct_ai > 0, 1, 0),
    ct_noai = ifelse(ct_noai > 0, 1, 0),
    other = ifelse(other > 0, 1, 0),

    # Strong variables - keep NAs as NAs
    "af_strong0" = case_when(
      is.na(`af_strong0`) ~ NA_real_,
      `af_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "af_strong1" = case_when(
      is.na(`af_strong1`) ~ NA_real_,
      `af_strong1` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_ai_strong0" = case_when(
      is.na(`ct_ai_strong0`) ~ NA_real_,
      `ct_ai_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_ai_strong1" = case_when(
      is.na(`ct_ai_strong1`) ~ NA_real_,
      `ct_ai_strong1` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_noai_strong0" = case_when(
      is.na(`ct_noai_strong0`) ~ NA_real_,
      `ct_noai_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_noai_strong1" = case_when(
      is.na(`ct_noai_strong1`) ~ NA_real_,
      `ct_noai_strong1` > 0 ~ 1,
      TRUE ~ 0
    ),

    # Interaction variables - keep NAs as NAs
    "af_ct_ai_strong0" = case_when(
      is.na(`af_ct_ai_strong0`) ~ NA_real_,
      `af_ct_ai_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "af_ct_ai_strong1" = case_when(
      is.na(`af_ct_ai_strong1`) ~ NA_real_,
      `af_ct_ai_strong1` > 0 ~ 1,
      TRUE ~ 0
    ),
    "af_ct_noai_strong0" = case_when(
      is.na(`af_ct_noai_strong0`) ~ NA_real_,
      `af_ct_noai_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "af_ct_noai_strong1" = case_when(
      is.na(`af_ct_noai_strong1`) ~ NA_real_,
      `af_ct_noai_strong1` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_ai_ct_noai_strong0" = case_when(
      is.na(`ct_ai_ct_noai_strong0`) ~ NA_real_,
      `ct_ai_ct_noai_strong0` > 0 ~ 1,
      TRUE ~ 0
    ),
    "ct_ai_ct_noai_strong1" = case_when(
      is.na(`ct_ai_ct_noai_strong1`) ~ NA_real_,
      `ct_ai_ct_noai_strong1` > 0 ~ 1,
      TRUE ~ 0
    )
  )

# ------------------------------------------------------------------------------
# Data Prep
# ------------------------------------------------------------------------------

# creating quantiles for institution controls
ecr_data <- ecr_data %>%
  mutate(
    institution_2yr_mean_citedness = factor(
      ntile(institution_2yr_mean_citedness, 4)
    ),
    institution_h_index = factor(ntile(institution_h_index, 4)),
    institution_i10_index = factor(ntile(institution_i10_index, 4)),
    institution_cited_by_count = factor(
      ntile(institution_cited_by_count, 4)
    ),
    organism_rarity_mean_quantile = factor(ntile(organism_rarity_mean, 4)),
    organism_rarity_max_quantile = factor(ntile(organism_rarity_max, 4)),
    mean_tmscore_quantile = factor(ntile(mean_tmscore, 4)),
    max_tmscore_quantile = factor(ntile(max_tmscore, 4))
  )

# fill with nan
ecr_data <- ecr_data %>%
  mutate(
    institution_type = ifelse(
      is.na(institution_type), "unknown", institution_type
    ),
    institution_country_code = ifelse(
      is.na(institution_country_code), "unknown", institution_country_code
    ),
    institution_2yr_mean_citedness = ifelse(
      is.na(institution_2yr_mean_citedness), "unknown",
      institution_2yr_mean_citedness
    ),
    institution_h_index = ifelse(
      is.na(institution_h_index), "unknown", institution_h_index
    ),
    institution_i10_index = ifelse(
      is.na(institution_i10_index), "unknown", institution_i10_index
    ),
    institution_cited_by_count = ifelse(
      is.na(institution_cited_by_count), "unknown",
      institution_cited_by_count
    ),
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    covid_share_2020 = ifelse(is.na(covid_share_2020), 0, covid_share_2020),
    num_uniprot_structures = ifelse(
      is.na(num_uniprot_structures), 0, num_uniprot_structures
    ),
    num_pdb_submissions = ifelse(is.na(pdb_submission), 0, pdb_submission),
    num_pdb_ids = ifelse(is.na(num_pdb_ids), 0, num_pdb_ids),
    num_primary_submissions = ifelse(is.na(num_primary_submissions), 0, num_primary_submissions), # nolint
    across(starts_with("field_"), ~ ifelse(is.na(.), 0, .)),
    across(starts_with("mesh_"), ~ ifelse(is.na(.), 0, .))
  )

# create translational variables
ecr_data <- ecr_data %>%
  mutate(
    num_uniprot_structures_w_disease = ifelse(
      num_uniprot_structures > 0 & num_diseases > 0, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_disease = ifelse(
      num_primary_submissions > 0 & num_diseases > 0, num_primary_submissions, 0 # nolint
    ),
    num_uniprot_structures_w_rare_organisms = ifelse(
      num_uniprot_structures > 0 & organism_rarity_mean_quantile == 4, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_rare_organisms = ifelse(
      num_primary_submissions > 0 & organism_rarity_mean_quantile == 4, num_primary_submissions, 0 # nolint
    ),
    num_uniprot_structures_w_low_similarity = ifelse(
      num_uniprot_structures > 0 & mean_tmscore_quantile == 1, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_low_similarity = ifelse(
      num_primary_submissions > 0 & mean_tmscore_quantile == 1, num_primary_submissions, 0 # nolint
    )
  )
# create factors, log transforms
ecr_data <- ecr_data %>%
  mutate(
    author = as.factor(author),
    author_position = as.factor(author_position),
    institution = as.factor(institution),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    ln1p_cit_norm_perc = log1p(percentile_value),
    logit_cit_norm_perc = log(
      percentile_value /
        (1 - percentile_value)
    ),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    ln1p_resolution = log1p(as.numeric(resolution_mean)),
    ln1p_R_free = log1p(as.numeric(R_free_mean)),
    ln1p_score = log1p(as.numeric(score_mean))
  )
# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
ecr_data$primary_field <- recode(ecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching) - Define matching groups
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
coarse_cols <- c(
  "ln1p_cited_by_count", "num_publications"
)

exact_cols <- c("institution_country_code")

mode_function <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Filter and prepare data for collapsing
quarterly_cem <- ecr_data %>%
  group_by(author) %>%
  mutate(af = max(af)) %>%
  ungroup() %>%
  filter(complete.cases(across(coarse_cols))) %>%
  filter(quarter %in% 200:204) %>% # (2020-2021)
  select(af, author, all_of(coarse_cols), all_of(exact_cols)) %>%
  group_by(author) %>%
  summarise(
    af = max(af),
    across(coarse_cols, \(x) mean(x, na.rm = TRUE)),
    across(exact_cols, mode_function)
  )

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching) - Match
# ------------------------------------------------------------------------------

match_out_af_coarse <- matchit(
  as.formula(paste0(
    "af ~ ",
    paste(coarse_cols, collapse = " + ")
  )),
  data = quarterly_cem, method = "cem", k2k = FALSE
)

# Store matched data for coarse_cols
cem_data_coarse <- match.data(match_out_af_coarse)

# Exact matching on the collapsed data using exact_cols
match_out_af_exact <- matchit(
  as.formula(paste0("af ~ ", paste0(exact_cols, collapse = " + "))),
  data = quarterly_cem, method = "exact"
)

# Store matched data for exact_cols
cem_data_exact <- match.data(match_out_af_exact)

# Combine the matched results
combined_cem_data <- intersect(
  cem_data_coarse$author,
  cem_data_exact$author
)

matched_data <- ecr_data %>%
  filter(author %in% combined_cem_data)

# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------
# Define sub_samples as a list of samples

matched_data <- matched_data %>%
  select(
    # Sample-defining variables
    "all_intents_high",

    # Basic identifiers and time
    "quarter",
    "author",

    # Institution variables
    "institution",
    "institution_type",
    "institution_country_code",
    "institution_cited_by_count",
    "institution_2yr_mean_citedness",
    "institution_h_index",
    "institution_i10_index",

    # Publication metrics
    "num_publications",
    "ln1p_cited_by_count",
    "ln1p_cit_0",
    "ln1p_cit_1",
    "ln1p_fwci",
    "patent_count",
    "patent_citation",
    "ca_count",

    # Covid
    "covid_share_2020",

    # Structure quality metrics
    "ln1p_resolution",
    "ln1p_R_free",
    "ln1p_score",

    # Basic usage
    "af",
    "ct_ai",
    "ct_noai",
    "other",

    # Strong usage
    "af_strong0",
    "af_strong1",
    "ct_ai_strong0",
    "ct_ai_strong1",
    "ct_noai_strong0",
    "ct_noai_strong1",

    # Strong usage interactions
    "af_ct_ai_strong0",
    "af_ct_ai_strong1",
    "af_ct_noai_strong0",
    "af_ct_noai_strong1",
    "ct_ai_ct_noai_strong0",
    "ct_ai_ct_noai_strong1",

    # Field and classification
    "primary_field",
    grep("^field_", names(matched_data), value = TRUE),
    grep("^mesh_", names(matched_data), value = TRUE),

    # Structure metrics
    "num_pdb_ids",
    "num_pdb_submissions",
    "num_primary_submissions",
    "num_uniprot_structures",

    # Disease and organism metrics
    "num_diseases",
    "organism_rarity_mean",
    "mean_tmscore",

    # Translational variables
    "num_uniprot_structures_w_disease",
    "num_primary_submissions_w_disease",
    "num_uniprot_structures_w_rare_organisms",
    "num_primary_submissions_w_rare_organisms",
    "num_uniprot_structures_w_low_similarity",
    "num_primary_submissions_w_low_similarity"
  )

colnames(matched_data) <- gsub(",", "", colnames(matched_data))

# Define sub_samples as a list of samples
sub_samples <- list()
sub_groups <- c("All PDB")
unique_scopes <- c("All", "Intent") # Changed to All/Intent distinction
unique_fields <- c(
  "All Fields",
  "Molecular Biology",
  "Medicine"
)

# Create subsets for all combinations of depth, field, and sub_group
for (scope_lvl in unique_scopes) {
  for (field in unique_fields) {
    for (sub_group in sub_groups) {
      sample_name <- paste0(
        "scope_", scope_lvl, "__field_",
        field, "__subgroup_", sub_group
      )
      message("Creating sample: ", sample_name)

      # Start with the full dataset
      sub_sample <- matched_data

      # # Apply depth filter
      if (scope_lvl == "Intent") {
        sub_sample <- subset(sub_sample, all_intents_high == TRUE)
      } else {
        # drop strong
        sub_sample <- sub_sample %>%
          select(-matches("strong[01]$|_strong"))
      }

      # Apply field filter
      if (field != "All Fields") {
        sub_sample <- subset(sub_sample, primary_field == field)
      }

      # Store the subset
      sub_samples[[sample_name]] <- sub_sample
    }
  }
}


# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
