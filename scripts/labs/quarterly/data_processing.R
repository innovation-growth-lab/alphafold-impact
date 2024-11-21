# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 500)

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
pathdir <- "data/05_model_output/labs/quarterly/data/"

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

credentials <- yaml.load_file("conf/base/credentials.yml")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key, # nolint
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret, # nolint
  "AWS_DEFAULT_REGION" = "eu-west-2" # nolint
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
foundational_path <- "oct/04_output/analysis/foundational_labs/staggered/outputs_quarterly.parquet" # nolint
applied_path <- "oct/04_output/analysis/applied_labs/staggered/outputs_quarterly.parquet" # nolint

# Fetch the data from the S3 bucket
foundational_labs_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = foundational_path,
  bucket = bucket
)

applied_labs_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = applied_path,
  bucket = bucket
)

# ------------------------------------------------------------------------------
# Function
# ------------------------------------------------------------------------------

process_high_pdb_labs_data <- function(data, quantile_threshold = 0.75) { # nolint
  filtered_data <- data %>% # nolint
    filter(time %in% -9:-1) # nolint

  pdb_count_data <- filtered_data %>% # nolint
    group_by(pi_id) %>% # nolint
    summarise(pdb_count = sum(!is.na(R_free))) %>% # nolint
    ungroup() # nolint

  pdb_count_data <- pdb_count_data %>% # nolint
    mutate( # nolint
      high_pdb = ifelse( # nolint
        pdb_count > quantile(pdb_count, quantile_threshold, na.rm = TRUE), # nolint
        1, 0
      )
    )

  processed_data <- data %>% # nolint
    left_join(pdb_count_data, by = "pi_id") # nolint

  return(processed_data)
}

# ------------------------------------------------------------------------------
# Data Prep - Foundational
# ------------------------------------------------------------------------------

foundational_labs_data <- process_high_pdb_labs_data(foundational_labs_data)

# create a new factor variable
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    strong_af = ifelse(is.na(strong_cumul_af), 0, strong_cumul_af),
    strong_ct_ai = ifelse(is.na(strong_cumul_ct_ai), 0, strong_cumul_ct_ai),
    strong_ct_noai = ifelse(is.na(strong_cumul_ct_noai), 0, strong_cumul_ct_noai), # nolint
    af_ind = ifelse(cum_af > 0, 1, 0),
    ct_ai_ind = ifelse(cum_ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(cum_ct_noai > 0, 1, 0),
    "af:ct_ai_ind" = ifelse(cum_af > 0 & cum_ct_ai > 0, 1, 0),
    "af:ct_noai_ind" = ifelse(cum_af > 0 & cum_ct_noai > 0, 1, 0),
    strong_af_ind = ifelse(strong_af > 0, 1, 0),
    strong_ct_ai_ind = ifelse(strong_ct_ai > 0, 1, 0),
    strong_ct_noai_ind = ifelse(strong_ct_noai > 0, 1, 0),
    "strong_af:strong_ct_ai_ind" = ifelse(strong_af > 0 & strong_ct_ai > 0, 1, 0), # nolint
    "strong_af:strong_ct_noai_ind" = ifelse(strong_af > 0 & strong_ct_noai > 0, 1, 0), # nolint
    af = cum_af,
    ct_ai = cum_ct_ai,
    ct_noai = cum_ct_noai,
    depth = "foundational"
  )

# relabel value alphafold in seed to af
foundational_labs_data$seed <- str_replace(
  foundational_labs_data$seed, "alphafold", "af"
)

# creating quantiles for institution controls
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    institution_2yr_mean_citedness = factor(
      ntile(institution_2yr_mean_citedness, 4)
    ),
    institution_h_index = factor(ntile(institution_h_index, 4)),
    institution_i10_index = factor(ntile(institution_i10_index, 4)),
    institution_cited_by_count = factor(
      ntile(institution_cited_by_count, 4)
    )
  )

# fill with nan
foundational_labs_data <- foundational_labs_data %>%
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
    institution = ifelse(is.na(institution_works_count), 0, institution_works_count), # nolint
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb)),
    covid_share_2020 = ifelse(is.na(covid_share_2020), 0, covid_share_2020)
  )

# create factors, log transforms
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    pi_id = as.factor(pi_id),
    tech_group = as.factor(seed),
    institution = as.factor(institution),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    logit_cit_norm_perc = log(
      citation_normalized_percentile_value /
        (1 - citation_normalized_percentile_value)
    ),
    num_pdb_submissions = pdb_submission,
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    ln1p_resolution = log1p(as.numeric(resolution)),
    ln1p_R_free = log1p(as.numeric(R_free)),
    quarter_year = as.factor(time)
  )

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
foundational_labs_data$primary_field <- recode(
  foundational_labs_data$primary_field, !!!field_mapping
)

# ------------------------------------------------------------------------------
# Data Prep - Applied
# ------------------------------------------------------------------------------

applied_labs_data <- process_high_pdb_labs_data(applied_labs_data)

# create a new factor variable
applied_labs_data <- applied_labs_data %>%
  mutate(
    strong_af = ifelse(is.na(strong_cumul_af), 0, strong_cumul_af),
    strong_ct_ai = ifelse(is.na(strong_cumul_ct_ai), 0, strong_cumul_ct_ai),
    strong_ct_noai = ifelse(is.na(strong_cumul_ct_noai), 0, strong_cumul_ct_noai), # nolint
    af_ind = ifelse(cum_af > 0, 1, 0),
    ct_ai_ind = ifelse(cum_ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(cum_ct_noai > 0, 1, 0),
    "af:ct_ai_ind" = ifelse(cum_af > 0 & cum_ct_ai > 0, 1, 0),
    "af:ct_noai_ind" = ifelse(cum_af > 0 & cum_ct_noai > 0, 1, 0),
    strong_af_ind = ifelse(strong_af > 0, 1, 0),
    strong_ct_ai_ind = ifelse(strong_ct_ai > 0, 1, 0),
    strong_ct_noai_ind = ifelse(strong_ct_noai > 0, 1, 0),
    "strong_af:strong_ct_ai_ind" = ifelse(strong_af > 0 & strong_ct_ai > 0, 1, 0), # nolint
    "strong_af:strong_ct_noai_ind" = ifelse(strong_af > 0 & strong_ct_noai > 0, 1, 0), # nolint
    af = cum_af,
    ct_ai = cum_ct_ai,
    ct_noai = cum_ct_noai,
    depth = "applied"
  )

# relabel value alphafold in seed to af
applied_labs_data$seed <- str_replace(
  applied_labs_data$seed, "alphafold", "af"
)

# creating quantiles for institution controls
applied_labs_data <- applied_labs_data %>%
  mutate(
    institution_2yr_mean_citedness = factor(
      ntile(institution_2yr_mean_citedness, 4)
    ),
    institution_h_index = factor(ntile(institution_h_index, 4)),
    institution_i10_index = factor(ntile(institution_i10_index, 4)),
    institution_cited_by_count = factor(
      ntile(institution_cited_by_count, 4)
    )
  )

# fill with nan
applied_labs_data <- applied_labs_data %>%
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
    institution = ifelse(is.na(institution_works_count), 0, institution_works_count), # nolint
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb)),
    covid_share_2020 = ifelse(is.na(covid_share_2020), 0, covid_share_2020)
  )

# create factors, log transforms
applied_labs_data <- applied_labs_data %>%
  mutate(
    pi_id = as.factor(pi_id),
    tech_group = as.factor(seed),
    institution = as.factor(institution),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    logit_cit_norm_perc = log(
      citation_normalized_percentile_value /
        (1 - citation_normalized_percentile_value)
    ),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    num_publications = num_publications,
    num_pdb_submissions = pdb_submission,
    primary_field = as.factor(primary_field),
    ln1p_resolution = log1p(as.numeric(resolution)),
    ln1p_R_free = log1p(as.numeric(R_free)),
    quarter_year = as.factor(time)
  )

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
applied_labs_data$primary_field <- recode(
  applied_labs_data$primary_field, !!!field_mapping
)

# ------------------------------------------------------------------------------
# Combine Foundational and Applied Data
# ------------------------------------------------------------------------------

quarterly_lab_data <- bind_rows(foundational_labs_data, applied_labs_data)

# last changes
quarterly_lab_data <- quarterly_lab_data %>%
  group_by(pi_id) %>%
  fill(
    institution_country_code, covid_share_2020,
    starts_with("institution_"),
    .direction = "downup"
  )

# drop rows with NA in the fields_, mesh_ columns
selected_columns <- quarterly_lab_data %>%
  select(
    quarter_year,
    primary_field,
    starts_with("field_"),
    starts_with("mesh_")
  )

# Filter rows with complete cases for the selected columns
quarterly_lab_data <- quarterly_lab_data %>%
  filter(
    complete.cases(
      across(c(
        quarter_year, primary_field, starts_with("field_"),
        starts_with("mesh_")
      ))
    )
  )

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching)
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
coarse_cols <- c(
  "ln1p_cited_by_count", "num_publications", "covid_share_2020"
)

exact_cols <- c("institution_country_code")

mode_function <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Filter and prepare data for collapsing
quarterly_cem <- quarterly_lab_data %>%
  group_by(pi_id) %>%
  mutate(af_ind = max(af_ind)) %>%
  ungroup() %>%
  filter(complete.cases(across(coarse_cols))) %>%
  filter(quarter_year %in% -9:-1) %>%
  select(af_ind, pi_id, all_of(coarse_cols), all_of(exact_cols)) %>%
  group_by(pi_id) %>%
  summarise(
    af_ind = max(af_ind),
    across(coarse_cols, \(x) mean(x, na.rm = TRUE)),
    across(exact_cols, mode_function)
  )

# ------------------------------------------------------------------------------
# Select Columns
# ------------------------------------------------------------------------------
# drop field_ and mesh_ columns
quarterly_lab_data <- quarterly_lab_data %>%
  select(
    "quarter_year",
    "pi_id",
    "institution",
    "institution_type",
    "institution_country_code",
    "institution_cited_by_count",
    "institution_2yr_mean_citedness",
    "institution_h_index",
    "institution_i10_index",
    "num_publications",
    "ln1p_cited_by_count",
    "ln1p_cit_0",
    "ln1p_cit_1",
    "ln1p_fwci",
    "logit_cit_norm_perc",
    "patent_count",
    "patent_citation",
    "ca_count",
    "ln1p_resolution",
    "ln1p_R_free",
    "num_pdb_submissions",
    "af_ind",
    "ct_ai_ind",
    "ct_noai_ind",
    "af",
    "ct_ai",
    "ct_noai",
    "strong_af_ind",
    "strong_ct_ai_ind",
    "strong_ct_noai_ind",
    "strong_af",
    "strong_ct_ai",
    "strong_ct_noai",
    "depth",
    "primary_field",
    "high_pdb",
    "covid_share_2020",
    grep("^field_", names(quarterly_lab_data), value = TRUE)
  )

colnames(quarterly_lab_data) <- gsub(",", "", colnames(quarterly_lab_data))

# free up memory, drop the individual dataframes
rm(foundational_labs_data, applied_labs_data)
gc()
# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------

# Define sub_samples as a list of samples
sub_samples <- list()
sub_groups <- c("All PDB - CEM", "High PDB - CEM")
unique_depths <- c("All Groups", "Foundational", "Applied")
unique_fields <- c(
  "All Fields",
  "Molecular Biology",
  "Medicine"
)

# Create subsets for all combinations of depth, field, and sub_group
for (depth_lvl in unique_depths) { # nolint
  for (field in unique_fields) {
    for (sub_group in sub_groups) {
      sample_name <- paste0(
        "depth_", depth_lvl, "__field_",
        field, "__subgroup_", sub_group
      )
      message("Creating sample: ", sample_name)

      # Start with the full dataset
      sub_sample <- quarterly_lab_data

      # Apply depth filter
      if (depth_lvl == "Foundational") {
        sub_sample <- subset(sub_sample, depth == "foundational")
      } else if (depth_lvl == "Applied") {
        sub_sample <- subset(sub_sample, depth == "applied")
      }

      # Apply field filter
      if (field != "All Fields") {
        sub_sample <- subset(sub_sample, primary_field == field)
      }

      # Apply sub_group filter
      if (grepl("High PDB", sub_group)) {
        sub_sample <- subset(sub_sample, high_pdb == 1)
      }

      if (grepl("CEM", sub_group)) {
        # CEM matching on the collapsed data using coarse_cols
        match_out_af_coarse <- matchit(
          as.formula(paste0(
            "af_ind ~ ",
            paste0(coarse_cols, collapse = " + ")
          )),
          data = quarterly_cem, method = "cem", k2k = FALSE
        )

        # Store matched data for coarse_cols
        cem_data_coarse <- match.data(match_out_af_coarse)

        # Exact matching on the collapsed data using exact_cols
        match_out_af_exact <- matchit(
          as.formula(paste0("af_ind ~ ", paste0(exact_cols, collapse = " + "))),
          data = quarterly_cem, method = "exact"
        )

        # Store matched data for exact_cols
        cem_data_exact <- match.data(match_out_af_exact)

        # Combine the matched results
        combined_cem_data <- intersect(
          cem_data_coarse$pi_id,
          cem_data_exact$pi_id
        )

        # Sample based on the combined matched group
        sub_sample <- sub_sample %>% filter(pi_id %in% combined_cem_data)
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
