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

credentials <- yaml.load_file("conf/base/credentials.yml")

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
# Data Prep
# ------------------------------------------------------------------------------

# create a new factor variable
ecr_data <- ecr_data %>%
  mutate(
    ct = ct_ai + ct_noai,
    af_ind = ifelse(af > 0, 1, 0),
    ct_ind = ifelse(ct > 0, 1, 0),
    ct_ai_ind = ifelse(ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(ct_noai > 0, 1, 0),
    "af:ct_ai_ind" = ifelse(af > 0 & ct_ai > 0, 1, 0),
    "af:ct_noai_ind" = ifelse(af > 0 & ct_noai > 0, 1, 0),
  )

# fill type, country_code, institution NA as "unknown"
ecr_data <- ecr_data %>%
  mutate(
    type = ifelse(is.na(type), "unknown", type),
    country_code = ifelse(is.na(country_code), "unknown", country_code),
    institution = ifelse(is.na(institution), "unknown", institution),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb))
  )

# create factors, log transforms
ecr_data <- ecr_data %>%
  mutate(
    author = as.factor(author),
    author_position = as.factor(author_position),
    strong = as.factor(
      if_else(chain_label %in% c("strong", "partial_strong"), 1,
        if_else(chain_label %in% c("no_data", "unknown"), NA_integer_, 0)
      )
    ),
    depth = as.factor(depth),
    institution = as.factor(institution),
    institution_type = as.factor(type),
    institution_country_code = as.factor(country_code),
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
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )
# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
ecr_data$primary_field <- recode(ecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------
# Define sub_samples as a list of samples
sub_samples <- list()
pdb_groups <- c("All PDB", "High PDB")
strength_groups <- c("General Use", "Methodological Use")
unique_depths <- c("All Groups", "Foundational", "Applied")
unique_fields <- c(
  "All Fields",
  "Molecular Biology",
  "Medicine"
)

# Create subsets for all combinations of depth, field, and pdb_group
for (depth_lvl in unique_depths) { # nolint
  for (field in unique_fields) {
    for (strength_group in strength_groups) {
      for (pdb_group in pdb_groups) {
        sample_name <- paste0(
          "depth_", depth_lvl, "__field_",
          field, "__use_", strength_group, "__pdb_", pdb_group
        )
        message("Creating sample: ", sample_name)

        # Start with the full dataset
        sub_sample <- ecr_data

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

        # Apply strength_group filter
        if (strength_group == "Methodological Use") {
          sub_sample <- subset(sub_sample, strong == 1)
        }

        # Apply pdb_group filter
        if (pdb_group == "High PDB") {
          sub_sample <- subset(sub_sample, high_pdb == 1)
        }

        # Store the subset
        sub_samples[[sample_name]] <- sub_sample
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
