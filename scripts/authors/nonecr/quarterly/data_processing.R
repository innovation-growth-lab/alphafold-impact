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
pathdir <- "data/05_model_output/authors/nonecr/quarterly/data/"

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
path <- "2025Q1/03_primary/nonecr/publications_quarterly.parquet" # nolint # mistake saved to oct

# Fetch the data from the S3 bucket
nonecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# Fix quarter
to_year_quarter <- function(code_vec) {
  year <- 1970 + (code_vec %/% 4)
  quarter <- (code_vec %% 4) + 1
  paste0(year, "Q", quarter)
}

# impute quarter
nonecr_data <- nonecr_data %>%
  mutate(quarter = to_year_quarter(quarter))

# drop obs after 2025Q1 (ie. 2025Q2)
nonecr_data <- nonecr_data %>%
  filter(quarter <= "2025Q1")

# ------------------------------------------------------------------------------
# Strong Data Prep
# ------------------------------------------------------------------------------

# replace core counts to correct cumul ones
nonecr_data <- nonecr_data %>%
  mutate(
    af = af_strong + af_weak + af_mixed + af_unknown,
    ct_ai = ct_ai_strong + ct_ai_weak + ct_ai_mixed + ct_ai_unknown,
    ct_pp = ct_pp_strong + ct_pp_weak + ct_pp_mixed + ct_pp_unknown,
    ct_sb = ct_sb_strong + ct_sb_weak + ct_sb_mixed + ct_sb_unknown,
    other = other_strong + other_weak + other_mixed + other_unknown
  )

# Calculate intent ratios for each author across all periods
author_intent_ratios <- nonecr_data %>%
  group_by(author) %>%
  summarise(
    af_intent_ratio = sum(af_with_intent, na.rm = TRUE) /
      (sum(af_with_intent, na.rm = TRUE) + sum(af_unknown, na.rm = TRUE)),
    ct_ai_intent_ratio = sum(ct_ai_with_intent, na.rm = TRUE) /
      (sum(ct_ai_with_intent, na.rm = TRUE) + sum(ct_ai_unknown, na.rm = TRUE)),
    ct_pp_intent_ratio = sum(ct_pp_with_intent, na.rm = TRUE) /
      (sum(ct_pp_with_intent, na.rm = TRUE) + sum(ct_pp_unknown, na.rm = TRUE)),
    ct_sb_intent_ratio = sum(ct_sb_with_intent, na.rm = TRUE) /
      (sum(ct_sb_with_intent, na.rm = TRUE) + sum(ct_sb_unknown, na.rm = TRUE)) # nolint
  )

# ------------------------------------------------------------------------------
# Create Hierarchical Intent Variables (all untriggered = 0)
# ------------------------------------------------------------------------------

create_hierarchical_intent <- function(data, tech_prefix) {
  strong_var <- paste0(tech_prefix, "_strong")
  mixed_var <- paste0(tech_prefix, "_mixed")
  weak_var <- paste0(tech_prefix, "_weak")

  intent_strong_var <- paste0(tech_prefix, "_intent_strong")
  intent_mixed_var <- paste0(tech_prefix, "_intent_mixed")
  intent_weak_var <- paste0(tech_prefix, "_intent_weak")

  data %>% # nolint
    group_by(author) %>% # nolint
    mutate( # nolint
      has_strong = any(!!sym(strong_var) > 0, na.rm = TRUE), # nolint
      has_mixed = !has_strong & any(!!sym(mixed_var) > 0, na.rm = TRUE), # nolint
      has_weak = !has_strong & !has_mixed & any(!!sym(weak_var) > 0, na.rm = TRUE), # nolint
      !!intent_strong_var := ifelse(has_strong & !!sym(strong_var) > 0, 1, 0), # nolint
      !!intent_mixed_var := ifelse(has_strong, 0, ifelse(has_mixed & !!sym(mixed_var) > 0, 1, 0)), # nolint
      !!intent_weak_var := ifelse(has_strong | has_mixed, 0, ifelse(has_weak & !!sym(weak_var) > 0, 1, 0)) # nolint
    ) %>%
    ungroup() %>% # nolint
    select(-has_strong, -has_mixed, -has_weak)
}

nonecr_data <- create_hierarchical_intent(nonecr_data, "af")
nonecr_data <- create_hierarchical_intent(nonecr_data, "ct_ai")
nonecr_data <- create_hierarchical_intent(nonecr_data, "ct_pp")
nonecr_data <- create_hierarchical_intent(nonecr_data, "ct_sb")

# ------------------------------------------------------------------------------
# Mask intent variables for authors with too much unknown (less than 1/3 intent)
# ------------------------------------------------------------------------------

mask_intent_for_unknown <- function(data) {
  # Get all technology prefixes
  tech_prefixes <- c("af", "ct_ai", "ct_pp", "ct_sb")

  # Calculate mask per author based on any technology having too much unknown
  mask_df <- data %>% # nolint
    group_by(author) %>% # nolint
    summarise(
      # For each technology, calculate the ratio of intent to total citations
      across( # nolint
        all_of(tech_prefixes), # nolint
        ~ {
          unknown_sum <- sum(get(paste0(cur_column(), "_unknown")), na.rm = TRUE) # nolint
          intent_sum <- sum(get(paste0(cur_column(), "_strong")), na.rm = TRUE) + # nolint
            sum(get(paste0(cur_column(), "_mixed")), na.rm = TRUE) +
            sum(get(paste0(cur_column(), "_weak")), na.rm = TRUE)
          total_sum <- unknown_sum + intent_sum
          # If no citations, return 1 (don't mask)
          if (total_sum == 0) 1 else intent_sum / total_sum
        }
      ),
      .groups = "drop"
    ) %>%
    # If any technology has ratio < 1/3, mask all intent variables
    mutate( # nolint
      mask = if_else( # nolint
        if_any(all_of(tech_prefixes), ~ . < 1 / 4), # nolint
        TRUE,
        FALSE
      )
    ) %>%
    select(author, mask) # nolint

  # Join mask and set all intent variables to NA if mask is TRUE
  data <- data %>% # nolint
    left_join(mask_df, by = "author") %>% # nolint
    mutate( # nolint1
      across( # nolint
        c(
          ends_with(c("_intent_strong", "_intent_mixed", "_intent_weak")), # nolint
          ends_with("_with_intent") # nolint
        ),
        ~ if_else(mask, NA_real_, .)
      )
    ) %>%
    select(-mask) # nolint

  data
}

# Apply the masking function once for all technologies
nonecr_data <- mask_intent_for_unknown(nonecr_data)

# ------------------------------------------------------------------------------
# Create Extensive Margin Variables
# ------------------------------------------------------------------------------

# create is_applied variable
nonecr_data <- nonecr_data %>%
  mutate(is_applied = ifelse(depth == "applied", 1, 0))

# create is_applied variable
nonecr_data <- nonecr_data %>%
  mutate(
    af = ifelse(af > 0, 1, 0),
    ct_ai = ifelse(ct_ai > 0, 1, 0),
    ct_pp = ifelse(ct_pp > 0, 1, 0),
    ct_sb = ifelse(ct_sb > 0, 1, 0),
    other = ifelse(other > 0, 1, 0),
    af_with_intent = ifelse(af_with_intent > 0, 1, 0),
    ct_ai_with_intent = ifelse(ct_ai_with_intent > 0, 1, 0),
    ct_pp_with_intent = ifelse(ct_pp_with_intent > 0, 1, 0),
    ct_sb_with_intent = ifelse(ct_sb_with_intent > 0, 1, 0),
    other_with_intent = ifelse(other_with_intent > 0, 1, 0)
  )

# ------------------------------------------------------------------------------
# Data Prep
# ------------------------------------------------------------------------------

# creating quantiles for institution controls
nonecr_data <- nonecr_data %>%
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
    max_tmscore_quantile = factor(ntile(max_tmscore, 4)),
    max_score_quantile = factor(ntile(max_score, 4)),
    max_fident_quantile = factor(ntile(max_fident, 4))
  )

# fill with nan
nonecr_data <- nonecr_data %>%
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
nonecr_data <- nonecr_data %>%
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
      num_uniprot_structures > 0 & max_tmscore_quantile == 1, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_low_similarity = ifelse(
      num_primary_submissions > 0 & max_tmscore_quantile == 1, num_primary_submissions, 0 # nolint
    )
  )
# create factors, log transforms
nonecr_data <- nonecr_data %>%
  mutate(
    author = as.factor(author),
    institution = as.factor(institution),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    ln1p_resolution = log1p(as.numeric(resolution)),
    ln1p_R_free = log1p(as.numeric(R_free)),
    ln1p_organism_rarity_mean = log1p(as.numeric(organism_rarity_mean)),
    ln1p_organism_rarity_max = log1p(as.numeric(organism_rarity_max)),
    ln1p_max_tmscore = log1p(as.numeric(max_tmscore)),
    ln1p_max_fident = log1p(as.numeric(max_fident)),
    ln1p_max_score = log1p(as.numeric(max_score)),
    ln1p_mesh_C = log1p(as.numeric(mesh_C)),
    year = as.integer(str_sub(quarter, 1, 4))
  )


# Calculate pre-2021 PDB submissions per author and identify top decile authors
# First calculate the submissions per author
author_submissions <- nonecr_data %>%
  group_by(author) %>%
  summarise(pre_2021_pdb_submissions = sum(num_pdb_submissions[year < 2021]))

# Calculate the 90th percentile threshold across all authors
submission_threshold <- quantile(
  author_submissions$pre_2021_pdb_submissions, 0.99,
  na.rm = TRUE
)

# Join back and create the indicator
nonecr_data <- nonecr_data %>%
  left_join(author_submissions, by = "author") %>%
  mutate(
    high_pdb_pre2021 = pre_2021_pdb_submissions > submission_threshold
  )

# Print diagnostics
message(paste0(
  "90th percentile threshold of pre-2021 PDB submissions: ",
  round(submission_threshold, 2)
))

message(paste0(
  "Number of unique authors with high_pdb_pre2021 = TRUE: ",
  nonecr_data %>%
    filter(high_pdb_pre2021 == TRUE) %>%
    pull(author) %>%
    n_distinct()
))


# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
nonecr_data$primary_field <- recode(nonecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching) - Define matching groups
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
coarse_cols <- c(
  "cited_by_count", "ln1p_fwci", "num_publications", "num_pdb_submissions",
  "field_biochemist", "field_chemistry", "field_medicine",
  "covid_share_2020"
)

pdb_cols <- c(
  "ln1p_max_tmscore", "ln1p_max_fident", "ln1p_max_score"
)

exact_cols <- c(
  "institution_country_code",
  "institution_h_index"
)

mode_function <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Filter and prepare data for collapsing
quarterly_cem <- nonecr_data %>%
  group_by(author) %>%
  mutate(treatment = max(af, ct_ai, ct_pp, ct_sb), af = max(af)) %>%
  ungroup() %>%
  filter(complete.cases(across(coarse_cols))) %>%
  filter(year < 2021) %>% # (2015-2020)
  select(
    af, treatment, author, other, high_pdb_pre2021,
    all_of(coarse_cols), all_of(exact_cols), all_of(pdb_cols)
  ) %>%
  group_by(author) %>%
  summarise(
    af = max(af),
    other = max(other),
    treatment = max(treatment),
    high_pdb_pre2021 = max(high_pdb_pre2021),
    across(coarse_cols, \(x) mean(x, na.rm = TRUE)),
    across(exact_cols, mode_function),
    across(pdb_cols, \(x) mean(x, na.rm = TRUE))
  )


cat("Before pre-filtering:", nrow(quarterly_cem), "authors\n")

# Remove extreme outliers (1%-99%) and low-publication authors
quarterly_cem <- quarterly_cem %>%
  filter(
    cited_by_count >= quantile(cited_by_count, 0.01, na.rm = TRUE),
    cited_by_count <= quantile(cited_by_count, 0.99, na.rm = TRUE),
    num_publications > 1
  )

cat("After pre-filtering:", nrow(quarterly_cem), "authors\n")

# Downsample controls if ratio > 10:1 to reduce memory
n_treated <- sum(quarterly_cem$af)
n_control <- sum(1 - quarterly_cem$af)
control_ratio <- n_control / n_treated

if (control_ratio > 10) {
  cat(
    "Control ratio", round(control_ratio, 1),
    ":1 too high, downsampling controls\n"
  )

  # Stratified downsampling to preserve distribution
  treated_data <- quarterly_cem %>% filter(af == 1)
  control_data <- quarterly_cem %>%
    filter(af == 0) %>%
    slice_sample(n = n_treated * 8) # Keep 8:1 ratio

  quarterly_cem <- bind_rows(treated_data, control_data)

  cat("After ratio'd pre-filtering:", nrow(quarterly_cem), "authors\n")
}



# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching) - Match
# ------------------------------------------------------------------------------

# ---- Treatment ----
match_out_af_coarse_treatment <- matchit(
  as.formula(paste0(
    "treatment ~ ",
    paste(coarse_cols, collapse = " + ")
  )),
  data = quarterly_cem, method = "cem", k2k = TRUE
)
message("Matched data for coarse_cols")
# Store matched data for coarse_cols
cem_data_coarse_treatment <- match.data(match_out_af_coarse_treatment)

# ---- AlphaFold ----
match_out_af_coarse <- matchit(
  as.formula(paste0(
    "af ~ ",
    paste(coarse_cols, collapse = " + ")
  )),
  data = quarterly_cem, method = "cem", k2k = TRUE
)
message("Matched data for af")
cem_data_coarse_af <- match.data(match_out_af_coarse)

# ---- Treatment - Exact ----
match_out_treatment_exact <- matchit(
  as.formula(paste0("treatment ~ ", paste0(exact_cols, collapse = " + "))),
  data = quarterly_cem, method = "exact"
)
message("Matched data for exact_cols")
# Store matched data for exact_cols
cem_data_exact_treatment <- match.data(match_out_treatment_exact)

# ---- PDB Submissions - Less Strict Matching ----

# quarterly_cem with no na values for the pdb_cols
quarterly_cem_pdb <- quarterly_cem %>%
  filter(complete.cases(across(pdb_cols)))

# Use quartiles instead of binary threshold for PDB matching
match_out_pdb <- matchit(
  as.formula(
    paste0(
      "high_pdb_pre2021 ~ ",
      paste0(pdb_cols, collapse = " + ")
    )
  ),
  data = quarterly_cem_pdb, method = "cem",
  k2k = FALSE
)
message("Matched data for pdb_cols")
cem_data_pdb <- match.data(match_out_pdb)

# Combine the matched results
combined_cem_data_treatment <- intersect(
  cem_data_coarse_treatment$author,
  cem_data_exact_treatment$author
)

# union with af
combined_cem_data <- union(
  combined_cem_data_treatment,
  cem_data_coarse_af$author
)

# union with af and pdb
combined_cem_data <- union(
  combined_cem_data,
  cem_data_pdb$author
)

matched_data <- nonecr_data %>%
  filter(author %in% combined_cem_data)


# Print number of unique authors for each tool type
message(paste0(
  "Number of AlphaFold authors: ",
  matched_data %>% filter(af == 1) %>% pull(author) %>% n_distinct()
))
message(paste0(
  "Number of AI tool authors: ",
  matched_data %>% filter(ct_ai == 1) %>% pull(author) %>% n_distinct()
))
message(paste0(
  "Number of protein prediction tool authors: ",
  matched_data %>% filter(ct_pp == 1) %>% pull(author) %>% n_distinct()
))
message(paste0(
  "Number of structure biology tool authors: ",
  matched_data %>% filter(ct_sb == 1) %>% pull(author) %>% n_distinct()
))
message(paste0(
  "Number of other authors: ",
  matched_data %>% filter(other == 1) %>% pull(author) %>% n_distinct()
))


# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------
# Define sub_samples as a list of samples

matched_data <- matched_data %>%
  select(
    # Sample-defining variables
    "is_applied",
    "high_pdb_pre2021",

    # Basic identifiers and time
    "quarter",
    "year",
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
    "cited_by_count",
    "patent_count",
    "patent_citation",
    "ca_count",
    "ln1p_fwci",
    "ln1p_mesh_C",

    # Covid
    "covid_share_2020",

    # Structure quality metrics
    "ln1p_resolution",
    "ln1p_R_free",

    # Basic usage
    "af",
    "ct_ai",
    "ct_pp",
    "ct_sb",
    "other",

    # Intent usage
    "af_with_intent",
    "ct_ai_with_intent",
    "ct_pp_with_intent",
    "ct_sb_with_intent",

    # Intent variables
    "af_intent_strong",
    "af_intent_weak",
    "af_intent_mixed",
    "ct_ai_intent_strong",
    "ct_ai_intent_weak",
    "ct_ai_intent_mixed",
    "ct_pp_intent_strong",
    "ct_pp_intent_weak",
    "ct_pp_intent_mixed",
    "ct_sb_intent_strong",
    "ct_sb_intent_weak",
    "ct_sb_intent_mixed",

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
    "ln1p_organism_rarity_max",
    "ln1p_organism_rarity_mean",
    "ln1p_max_score",
    "ln1p_max_tmscore",
    "ln1p_max_fident",
    "normalised_max_score",
    "normalised_max_tmscore",
    "normalised_max_fident",

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
sub_groups <- c("All PDB", "High PDB")
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
      if (scope_lvl == "All") {
        # drop strong
        sub_sample <- sub_sample %>%
          select(-matches("strong[01]$|_strong"))
      }

      # Apply field filter
      if (field != "All Fields") {
        sub_sample <- subset(sub_sample, primary_field == field)
      }

      # Apply sub_group filter
      if (grepl("High PDB", sub_group)) {
        sub_sample <- subset(sub_sample, high_pdb_pre2021 == TRUE)
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
