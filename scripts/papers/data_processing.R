# Clean out the workspace
rm(list = ls())
options(max.print = 1000)

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
pathdir <- "data/05_model_output/papers/"

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
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key,
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret,
  "AWS_DEFAULT_REGION" = "eu-west-2"
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
path <- "oct/04_output/publications/regression_inputs.parquet" # nolint


# Fetch the data from the S3 bucket
papers <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# Replace commas in column names
colnames(papers) <- gsub(",", "", colnames(papers))

# for duplicate ids, keep first. fill few nas
papers <- papers %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    type = ifelse(is.na(type), "unknown", type),
    country_code = ifelse(is.na(country_code), "unknown", country_code),
    institution = ifelse(is.na(institution), "unknown", institution),
    pdb_submission = ifelse(pdb_submission > 0, 1, 0)
  )

# Create 'strong' variable
papers <- papers %>%
  mutate(
    strong = if_else(chain_label %in% c("strong", "partial_strong"), 1, 0),
    # strong = replace_na(strong, 0), # nolint
    strong = as.factor(strong),
    depth = as.factor(if_else(level == 0, "foundational", "applied")),
    af = if_else(source == "af", 1, 0),
    ct = if_else(source %in% c("ct_ai", "ct_noai"), 1, 0),
    ct_ai = if_else(source == "ct_ai", 1, 0),
    ct_noai = if_else(source == "ct_noai", 1, 0),
    institution = as.factor(last_author_institution),
    institution_type = as.factor(type),
    institution_country_code = as.factor(country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci), # nolint
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    logit_cit_norm_perc = log(
      citation_normalized_percentile_value /
        (1 - citation_normalized_percentile_value)
    ),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field), # nolint
    publication_date = as.Date(publication_date),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free),
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = paste0(
      year(publication_date), " Q", quarter(publication_date)
    ),
  ) %>%
  rename(author = last_author, tech_group = source)


pubs_per_quarter <- papers %>%
  group_by(author, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()


# Merge the summary back into the original papers DataFrame
papers <- papers %>%
  left_join(pubs_per_quarter, by = c("author", "quarter_year"))

# Fill NA values in 'field_' and 'mesh_' prefix columns with 0
papers <- papers %>%
  mutate(across(starts_with("field_"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("mesh_"), ~ replace_na(., 0)))

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

papers$primary_field <- recode(papers$primary_field, !!!field_mapping)

# Create 'high_pdb' variable
papers <- papers %>%
  mutate(high_pdb = if_else(
    group_pdb_count >= quantile(
      group_pdb_count, 0.75,
      na.rm = TRUE
    ),
    1,
    0
  )) # nolint

# ------------------------------------------------------------------------------
# SUBSET
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
        sub_sample <- papers

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
