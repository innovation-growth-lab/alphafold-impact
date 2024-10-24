# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)

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
pathdir <- "data/05_model_output/ecr/"

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
path <- "oct/03_primary/ecr/publications.parquet" # nolint

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
    af_ct = ifelse(af > 0 & ct > 0, af + ct, 0),
    is_af = af > 0 & ct == 0,
    is_ct = ct > 0 & af == 0,
    is_af_ct = af_ct > 0,
    is_other = !(is_af | is_ct | is_af_ct)
  ) %>%
  mutate(
    af = ifelse(af_ct > 0, 0, af),
    ct = ifelse(af_ct > 0, 0, ct)
  )

# create factors, log transforms
ecr_data <- ecr_data %>%
  mutate(
    author = as.factor(author),
    author_position = as.factor(author_position),
    depth = as.factor(depth),
    institution = as.factor(institution),
    institution_type = as.factor(type),
    institution_country_code = as.factor(country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    ln1p_cit_2 = log1p(cit_2),
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
pubs_per_quarter <- ecr_data %>%
  group_by(author, quarter_year) %>%
  summarize(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original ecr_data DataFrame
ecr_data <- ecr_data %>%
  left_join(pubs_per_quarter, by = c("author", "quarter_year"))

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "biochem_genetics_molecular_biology", # nolint
  "Medicine" = "medicine",
  "Chemistry" = "chemistry",
  "Agricultural and Biological Sciences" = "agricultural_biological_sciences",
  "Immunology and Microbiology" = "immunology_microbiology",
  "Health Professions" = "health_professions"
)

# Apply the mapping to the primary_field column
ecr_data$primary_field <- recode(ecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------

# Define sub_samples as a list of samples
sub_samples <- list(
  depth_all__field_all = ecr_data,
  depth_foundational__field_all = subset(ecr_data, depth == "foundational"),
  depth_applied__field_all = subset(ecr_data, depth == "applied")
)

unique_depths <- c("all", "foundational", "applied")
unique_fields <- c(
  "biochem_genetics_molecular_biology", "medicine", "chemistry",
  "agricultural_biological_sciences", "immunology_microbiology",
  "health_professions"
)

for (depth_lvl in unique_depths) {
  for (field in unique_fields) {
    sample_name <- paste0("depth_", depth_lvl, "__field_", field)
    if (depth_lvl == "all") {
      sub_samples[[sample_name]] <- subset(
        ecr_data, primary_field == field
      )
    } else {
      sub_samples[[sample_name]] <- subset(
        ecr_data, depth == depth_lvl & primary_field == field
      )
    }
  }
}

# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
