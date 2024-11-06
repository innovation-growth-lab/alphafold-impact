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
pathdir <- "data/05_model_output/labs/articles/data/"

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
foundational_path <- "oct/04_output/analysis/foundational_labs/staggered/outputs_primary.parquet" # nolint
applied_path <- "oct/04_output/analysis/applied_labs/staggered/outputs_primary.parquet" # nolint

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

process_high_pdb_labs_data <- function(data, year_threshold = 2018, quantile_threshold = 0.75) { # nolint
  filtered_data <- data %>% # nolint
    filter(year < year_threshold) # nolint

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

foundational_labs_data <- foundational_labs_data %>%
  mutate(
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = as.factor(publication_date)
  )

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

# fill with nan
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    institution_type = ifelse(
      is.na(institution_type), "unknown", institution_type
    ),
    institution_country_code = ifelse(
      is.na(institution_country_code), "unknown", institution_country_code
    ),
    institution = ifelse(is.na(institution), "unknown", institution),
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb))
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
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )

# Compute the per-pi_id number of publications at each given quarter_year
pubs_per_quarter <- foundational_labs_data %>%
  group_by(pi_id, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original foundational_labs_data DataFrame
foundational_labs_data <- foundational_labs_data %>%
  left_join(pubs_per_quarter, by = c("pi_id", "quarter_year"))

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

applied_labs_data <- applied_labs_data %>%
  mutate(
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = as.factor(publication_date)
  )

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

# fill with nan
applied_labs_data <- applied_labs_data %>%
  mutate(
    institution_type = ifelse(
      is.na(institution_type), "unknown", institution_type
    ),
    institution_country_code = ifelse(
      is.na(institution_country_code), "unknown", institution_country_code
    ),
    institution = ifelse(is.na(institution), "unknown", institution),
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb))
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
    primary_field = as.factor(primary_field),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )

# Compute the per-pi_id number of publications at each given quarter_year
pubs_per_quarter <- applied_labs_data %>%
  group_by(pi_id, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original applied_labs_data DataFrame
applied_labs_data <- applied_labs_data %>%
  left_join(pubs_per_quarter, by = c("pi_id", "quarter_year"))

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

papers_lab_data <- bind_rows(foundational_labs_data, applied_labs_data)

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching)
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
cols <- c(
  "cited_by_count", "fwci", "citation_normalized_percentile_value",
  "patent_count"
)

# Filter and prepare data for collapsing
papers_cem <- papers_lab_data %>%
  group_by(pi_id) %>%
  mutate(af_ind = max(af_ind)) %>%
  ungroup() %>%
  filter(complete.cases(fwci, citation_normalized_percentile_value)) %>%
  filter(quarter_year %in% c(
    "2017Q1", "2017Q2", "2017Q3", "2017Q4",
    "2018Q1", "2018Q2", "2018Q3", "2018Q4"
  )) %>%
  select(af_ind, pi_id, all_of(cols)) %>%
  group_by(pi_id) %>%
  summarise(
    af_ind = max(af_ind),
    across(cols, \(x) mean(x, na.rm = TRUE))
  )

# ------------------------------------------------------------------------------
# Select Columns
# ------------------------------------------------------------------------------
# drop field_ and mesh_ columns
papers_lab_data <- papers_lab_data %>%
  select(
    "num_publications",
    "quarter_year",
    "pi_id",
    "institution",
    "institution_type",
    "institution_country_code",
    "ln1p_cited_by_count",
    "ln1p_cit_0",
    "ln1p_cit_1",
    "ln1p_fwci",
    "logit_cit_norm_perc",
    "ln1p_patent_count",
    "ln1p_patent_citation",
    "ln1p_ca_count",
    "resolution",
    "R_free",
    "pdb_submission",
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
    grep("^field_", names(papers_lab_data), value = TRUE),
    grep("^mesh_", names(papers_lab_data), value = TRUE)
  )

# free up memory, drop the individual dataframes
rm(foundational_labs_data, applied_labs_data)
gc()
# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------

# Define sub_samples as a list of samples
sub_samples <- list()
sub_groups <- c("All PDB", "High PDB", "CEM")
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
      sub_sample <- papers_lab_data

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
      if (sub_group == "High PDB") {
        sub_sample <- subset(sub_sample, high_pdb == 1)
      }

      if (sub_group == "CEM") {
        # CEM matching on the collapsed data
        match_out_af <- matchit(
          as.formula(paste0("af_ind ~ ", paste0(cols, collapse = " + "))),
          data = papers_cem, method = "cem", k2k = TRUE
        )

        # Store matched data
        cem_data <- match.data(match_out_af)

        # Sample based on the matched group
        qtly_cem <- papers_lab_data %>% semi_join(cem_data, by = "pi_id")
        sub_sample <- qtly_cem

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
