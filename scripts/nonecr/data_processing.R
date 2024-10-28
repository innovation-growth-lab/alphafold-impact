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
pathdir <- "data/05_model_output/nonecr/"

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
path <- "oct/03_primary/nonecr/publications.parquet" # nolint

# Fetch the data from the S3 bucket
nonecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# ------------------------------------------------------------------------------
# Data Prep
# ------------------------------------------------------------------------------

# create a new factor variable
nonecr_data <- nonecr_data %>%
  mutate(
    ct = ct_ai + ct_noai,
    af_ind = ifelse(af > 0, 1, 0),
    ct_ind = ifelse(ct > 0, 1, 0),
    ct_ai_ind = ifelse(ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(ct_noai > 0, 1, 0),
  )

# create factors, log transforms
nonecr_data <- nonecr_data %>%
  mutate(
    author = as.factor(author),
    author_position = as.factor(author_position),
    depth = as.factor(depth),
    institution = as.factor(institution),
    institution_type = as.factor(type),
    institution_country_code = as.factor(country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    ln1p_cit_2 = log1p(cit_2),
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    publication_date = as.Date(publication_date),
    quarter_year = paste0(
      year(publication_date), " Q", quarter(publication_date)
    ),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )

# Compute the per-author number of publications at each given quarter_year
pubs_per_quarter <- nonecr_data %>%
  group_by(author, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original nonecr_data DataFrame
nonecr_data <- nonecr_data %>%
  left_join(pubs_per_quarter, by = c("author", "quarter_year"))

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "biochem_genetics_molecular_biology", # nolint
  "Medicine" = "medicine",
  "Chemistry" = "chemistry",
  "Immunology and Microbiology" = "immunology_microbiology"
)

# Apply the mapping to the primary_field column
nonecr_data$primary_field <- recode(nonecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------
# Define sub_samples as a list of samples
sub_samples <- list(
  depth_all__field_all__tech_all = nonecr_data
)

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
    return(data %>% filter(af >= 0, ct_ai >= 0, ct_noai == 0)) # nolint
  } else if (tech_group == "ct_noai") {
    return(data %>% filter(af >= 0, ct_noai >= 0, ct_ai == 0)) # nolint
  }
}

# Create subsets for all combinations of depth, field, and tech_group
for (depth_lvl in unique_depths) {
  for (field in c("all", unique_fields)) {
    for (tech_group in tech_groups) {
      sample_name <- paste0(
        "depth_", depth_lvl, "__field_", field, "__tech_", tech_group
      )
      if (field == "all") {
        if (depth_lvl == "all") {
          sub_samples[[sample_name]] <- create_tech_group_subset(
            nonecr_data, tech_group
          )
        } else {
          sub_samples[[sample_name]] <- create_tech_group_subset(
            subset(nonecr_data, depth == depth_lvl), tech_group
          )
        }
      } else {
        if (depth_lvl == "all") {
          sub_samples[[sample_name]] <- create_tech_group_subset(
            subset(nonecr_data, primary_field == field), tech_group
          )
        } else {
          sub_samples[[sample_name]] <- create_tech_group_subset(
            subset(nonecr_data, depth == depth_lvl & primary_field == field),
            tech_group
          )
        }
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
