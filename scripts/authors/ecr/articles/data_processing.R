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
pathdir <- "data/05_model_output/authors/ecr/articles/data/"

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
    strong_af = ifelse(is.na(strong_af), 0, strong_af),
    strong_ct_ai = ifelse(is.na(strong_ct_ai), 0, strong_ct_ai),
    strong_ct_noai = ifelse(is.na(strong_ct_noai), 0, strong_ct_noai),
    ct = ct_ai + ct_noai,
    af_ind = ifelse(af > 0, 1, 0),
    ct_ind = ifelse(ct > 0, 1, 0),
    ct_ai_ind = ifelse(ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(ct_noai > 0, 1, 0),
    "af:ct_ai_ind" = ifelse(af > 0 & ct_ai > 0, 1, 0),
    "af:ct_noai_ind" = ifelse(af > 0 & ct_noai > 0, 1, 0),
    strong_af_ind = ifelse(strong_af > 0, 1, 0),
    strong_ct_ai_ind = ifelse(strong_ct_ai > 0, 1, 0),
    strong_ct_noai_ind = ifelse(strong_ct_noai > 0, 1, 0),
    "strong_af:strong_ct_ai_ind" = ifelse(strong_af > 0 & strong_ct_ai > 0, 1, 0), # nolint
    "strong_af:strong_ct_noai_ind" = ifelse(strong_af > 0 & strong_ct_noai > 0, 1, 0) # nolint
  )

# fill type, country_code, institution NA as "unknown"
ecr_data <- ecr_data %>%
  mutate(
    type = ifelse(is.na(type), "unknown", type),
    country_code = ifelse(is.na(country_code), "unknown", country_code),
    institution = ifelse(is.na(institution), "unknown", institution),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb)),
    pdb_submission = ifelse(pdb_submission > 0, 1, 0)
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
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    ln1p_cit_2 = log1p(cit_2),
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    logit_cit_norm_perc = log(
      citation_normalized_percentile_value /
        (1 - citation_normalized_percentile_value)
    ),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    publication_date = as.Date(publication_date),
    quarter_year = paste0(
      year(publication_date), " Q", quarter(publication_date)
    ),
    year = as.integer(str_sub(publication_date, 1, 4)),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )

# Compute the per-author number of publications at each given quarter_year
pubs_per_quarter <- ecr_data %>%
  group_by(author, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original ecr_data DataFrame
ecr_data <- ecr_data %>%
  left_join(pubs_per_quarter, by = c("author", "quarter_year"))

# Define the mapping o
# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
ecr_data$primary_field <- recode(ecr_data$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching)
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
cols <- c(
  "num_publications", "cited_by_count", "cit_0", "cit_1", "patent_count"
)

# Filter and prepare data for collapsing
ecr_data_cem <- ecr_data %>%
  group_by(author) %>%
  mutate(af_ind = max(af_ind)) %>%
  ungroup() %>%
  filter(complete.cases(cit_0, cit_1)) %>%
  filter(quarter_year %in% c("2020 Q1", "2020 Q2", "2020 Q3", "2020 Q4")) %>%
  select(af_ind, author, cols) %>%
  group_by(author) %>%
  summarise(
    af_ind = max(af_ind),
    across(cols, \(x) mean(x, na.rm = TRUE))
  )

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

ecr_data <- ecr_data %>%
  select(c(
    "num_publications",
    "quarter_year",
    "author",
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
    "high_pdb"
  ))

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

      # Apply sub_group filter
      if (sub_group == "High PDB") {
        sub_sample <- subset(sub_sample, high_pdb == 1)
      }

      if (sub_group == "CEM") {
        # CEM matching on the collapsed data
        match_out_af <- matchit(
          as.formula(paste0("af_ind ~ ", paste0(cols, collapse = " + "))),
          data = ecr_data_cem, method = "cem", k2k = TRUE
        )

        # Store matched data
        cem_data <- match.data(match_out_af)

        # Sample based on the matched group
        qtly_cem <- ecr_data %>% semi_join(cem_data, by = "author")
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
