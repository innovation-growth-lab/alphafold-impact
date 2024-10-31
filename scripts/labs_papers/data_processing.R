# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 100)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "aws.s3", "yaml", "zoo"
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
pathdir <- "data/05_model_output/labs_papers/"

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
# Data Prep - Foundational
# ------------------------------------------------------------------------------

# create a new factor variable
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    ct = cum_ct_ai + cum_ct_noai,
    af_ind = ifelse(cum_af > 0, 1, 0),
    ct_ai_ind = ifelse(cum_ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(cum_ct_noai > 0, 1, 0),
    af = cum_af,
    ct = cum_ct_ai + cum_ct_noai,
    depth = "foundational"
  )

# relabel value alphafold in seed to af
foundational_labs_data$seed <- str_replace(foundational_labs_data$seed, "alphafold", "af")

# fill with nan
foundational_labs_data <- foundational_labs_data %>%
  mutate(
    ca_count = ifelse(is.na(ca_count), 0, ca_count)
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
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = as.factor(publication_date),
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
  "Biochemistry, Genetics and Molecular Biology" = "biochem_genetics_molecular_biology", # nolint
  "Medicine" = "medicine",
  "Chemistry" = "chemistry",
  "Immunology and Microbiology" = "immunology_microbiology"
)

# Apply the mapping to the primary_field column
foundational_labs_data$primary_field <- recode(
  foundational_labs_data$primary_field, !!!field_mapping
)

# Filter the DataFrame to include only observations with quarter_year < 2018
filtered_data <- foundational_labs_data %>%
  filter(year < 2018)

# Group by pi_id and count the non-NA R_free values for each pi_id
pdb_count_data <- filtered_data %>%
  group_by(pi_id) %>%
  summarise(pdb_count = sum(!is.na(R_free))) %>%
  ungroup()

pdb_count_data <- pdb_count_data %>%
  mutate(
    high_pdb = ifelse(
      pdb_count > quantile(pdb_count, 0.75, na.rm = TRUE),
      1, 0
    )
  )

# Merge the resulting count back into the original DataFrame
foundational_labs_data <- foundational_labs_data %>%
  left_join(pdb_count_data, by = "pi_id")

# ------------------------------------------------------------------------------
# Data Prep - Applied
# ------------------------------------------------------------------------------

# create a new factor variable
applied_labs_data <- applied_labs_data %>%
  mutate(
    ct = cum_ct_ai + cum_ct_noai,
    af_ind = ifelse(cum_af > 0, 1, 0),
    ct_ai_ind = ifelse(cum_ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(cum_ct_noai > 0, 1, 0),
    depth = "applied",
    af = cum_af,
    ct = cum_ct_ai + cum_ct_noai
  )

# relabel value alphafold in seed to af
applied_labs_data$seed <- str_replace(applied_labs_data$seed, "alphafold", "af")

applied_labs_data <- applied_labs_data %>%
  mutate(
    ca_count = ifelse(is.na(ca_count), 0, ca_count)
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
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = publication_date,
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
  "Biochemistry, Genetics and Molecular Biology" = "biochem_genetics_molecular_biology", # nolint
  "Medicine" = "medicine",
  "Chemistry" = "chemistry",
  "Immunology and Microbiology" = "immunology_microbiology"
)

# Apply the mapping to the primary_field column
applied_labs_data$primary_field <- recode(
  applied_labs_data$primary_field, !!!field_mapping
)

# work to create high pdb normal pdb
filtered_data <- applied_labs_data %>%
  filter(year < 2018)

# Group by pi_id and count the non-NA R_free values for each pi_id
pdb_count_data <- filtered_data %>%
  group_by(pi_id) %>%
  summarise(pdb_count = sum(!is.na(R_free))) %>%
  ungroup()

pdb_count_data <- pdb_count_data %>%
  mutate(
    high_pdb = ifelse(
      pdb_count > quantile(pdb_count, 0.75, na.rm = TRUE),
      1, 0
    )
  )

# Merge the resulting count back into the original DataFrame
applied_labs_data <- applied_labs_data %>%
  left_join(pdb_count_data, by = "pi_id")

# ------------------------------------------------------------------------------
# Combine Foundational and Applied Data
# ------------------------------------------------------------------------------

papers_lab_data <- bind_rows(foundational_labs_data, applied_labs_data)

# drop field_ and mesh_ columns
papers_lab_data <- papers_lab_data %>%
  select(
    -starts_with("field_"),
    -starts_with("primary_field_")
  )

# free up memory, drop the individual dataframes
rm(foundational_labs_data, applied_labs_data)
# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------

# Define sub_samples as a list of samples
sub_samples <- list()
pdb_groups <- c("all", "high")
tech_groups <- c("all", "ct_ai", "ct_noai")
unique_depths <- c("all", "foundational", "applied")
unique_fields <- c(
  "biochem_genetics_molecular_biology",
  "medicine",
  "chemistry",
  "immunology_microbiology"
)

create_tech_group_subset <- function(data, tech_group) {
  if (tech_group == "all") {
    return(data)
  } else if (tech_group == "ct_ai") {
    return(data %>% filter(cum_af >= 0, cum_ct_ai >= 0, cum_ct_noai == 0)) # nolint
  } else if (tech_group == "ct_noai") {
    return(data %>% filter(cum_af >= 0, cum_ct_noai >= 0, cum_ct_ai == 0)) # nolint
  }
}

# Create subsets for all combinations of depth, field, tech_group, and pdb_group
for (depth_lvl in unique_depths) { # nolint
  for (field in c("all", unique_fields)) {
    for (tech_group in tech_groups) {
      for (pdb_group in pdb_groups) {
        sample_name <- paste0(
          "depth_", depth_lvl, "__field_",
          field, "__tech_", tech_group, "__pdb_", pdb_group
        )
        if (field == "all") {
          if (depth_lvl == "all") {
            sub_sample <- create_tech_group_subset(papers_lab_data, tech_group)
          } else {
            sub_sample <- create_tech_group_subset(
              subset(papers_lab_data, depth == depth_lvl), tech_group
            )
          }
        } else {
          if (depth_lvl == "all") {
            sub_sample <- create_tech_group_subset(
              subset(papers_lab_data, primary_field == field), tech_group
            )
          } else {
            sub_sample <- create_tech_group_subset(
              subset(papers_lab_data, depth == depth_lvl & primary_field == field), # nolint
              tech_group
            )
          }
        }
        # Filter based on pdb_group
        if (pdb_group == "high") {
          sub_sample <- subset(sub_sample, high_pdb == 1)
        }
        sub_samples[[sample_name]] <- sub_sample
      }
    }
  }
}
# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
