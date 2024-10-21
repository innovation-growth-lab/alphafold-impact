# Clean out the workspace
rm(list = ls())
options(max.print = 1000)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "MatchIt", "fastDummies", "aws.s3",
  "yaml", "zoo", "fixest"
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
figures <- "data/05_model_output/figures/ecr/"
tables <- "data/05_model_output/tables/ecr/"

# Create directories if they do not exist
if (!dir.exists(figures)) {
  dir.create(figures, recursive = TRUE)
}

if (!dir.exists(tables)) {
  dir.create(tables, recursive = TRUE)
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
path <- "oct/04_output/ecr/regression_inputs.parquet" # nolint

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
    category = case_when(
      af > 0 & ct == 0 ~ "af",
      ct > 0 & af == 0 ~ "ct",
      af > 0 & ct > 0 ~ "af_ct",
      TRUE ~ "other"
    )
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
    )
  )

# Define sub_samples as a list of samples
sub_samples <- list(
  all = ecr_data,
  depth_foundational = subset(ecr_data, depth == "foundational"),
  depth_applied = subset(ecr_data, depth == "applied")
)

# Add interactions of every depth with every unique primary_field
unique_depths <- unique(ecr_data$depth)
unique_fields <- unique(ecr_data$primary_field)

for (depth in unique_depths) {
  for (field in unique_fields) {
    sample_name <- paste0("depth_", depth, "_field_", field)
    sub_samples[[sample_name]] <- subset(
      ecr_data, depth == depth & primary_field == field
    )
  }
}

covs <- list()
covs[["base0"]] <- c("primary_field", "author_position")

fes <- list()
fes[["fe0"]] <- c("author", "quarter_year")
fes[["fe1"]] <- c(
  "quarter_year", "institution", "institution_type",
  "institution_country_code", "author"
)

cov_sets <- c("base0")
fe_list <- c("fe0", "fe1")
dep_vars <- c(
  "ln1p_cited_by_count", "ln1p_cit_0", "ln1p_cit_1",
  "fwci", "citation_normalized_percentile_value",
  "ln1p_patent_count", "ln1p_patent_citation", "ln1p_ca_count",
  "resolution", "R_free"
)
treat_vars <- c(
  "category",
  "category + category:af + category:ct + category:af_ct",
  "category + category:af + category:ct + category:af_ct + category:depth"
)

form_list <- list()
# Iterate over dependent variables
for (dep_var in dep_vars) { # nolint
  # Iterate over covariate sets
  for (cov_set in cov_sets) {
    # Iterate over fixed effects
    for (fe in fe_list) {
      # Iterate over treatment variables
      for (treat_var in treat_vars) {
        # Check if covs[[cov_set]] is empty
        if (length(covs[[cov_set]]) == 0) {
          # Create formula without '+' before '|'
          form_list[[
            paste0(
              dep_var, "_", cov_set, "_", fe, "_", gsub(" ", "_", treat_var)
            )
          ]] <- as.formula(
            paste0(
              dep_var, " ~ ", treat_var, " |",
              paste0(fes[[fe]], collapse = " + ")
            )
          )
        } else {
          # Create formula with '+' before '|'
          form_list[[
            paste0(
              dep_var, "_", cov_set, "_", fe, "_", gsub(" ", "_", treat_var)
            )
          ]] <- as.formula(
            paste0(
              dep_var, " ~ ", treat_var, " +",
              paste0(covs[[cov_set]], collapse = " + "),
              "|", paste0(fes[[fe]], collapse = " + ")
            )
          )
        }
      }
    }
  }
}

results <- list()
# For each subset, compute feols
for (sub in names(sub_samples)) {
  # For each formula, compute feols
  for (form in names(form_list)) {
    regression_label <- paste0(sub, "_", form)
    message("Running regression: ", regression_label)
    results[[regression_label]] <- feols(
      form_list[[form]],
      data = sub_samples[[sub]],
    )
  }
}

# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------

variable_interest <- c(
  "category", "category:af", "category:ct", "category:af_ct",
  "category:depth"
)

table_info <- list(
  "ln1p_cited_by_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cited_by_count.tex"
  ),
  "ln1p_cit_0" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cit_0.tex"
  ),
  "ln1p_cit_1" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cit_1.tex"
  ),
  "fwci" = list(
    vars_to_keep = variable_interest,
    file_name = "fwci.tex"
  ),
  "citation_normalized_percentile_value" = list(
    vars_to_keep = variable_interest,
    file_name = "citation_normalized_percentile_value.tex"
  ),
  "ln1p_patent_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_patent_count.tex"
  ),
  "ln1p_patent_citation" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_patent_citation.tex"
  ),
  "ln1p_ca_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_ca_count.tex"
  ),
  "resolution" = list(
    vars_to_keep = variable_interest,
    file_name = "resolution.tex"
  ),
  "R_free" = list(
    vars_to_keep = variable_interest,
    file_name = "R_free.tex"
  )
)

subsets <- names(sub_samples)

# Function to generate tables
generate_tables <- function(dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint
  for (dep_var in dep_vars) {
    vars_to_keep <- table_info[[dep_var]]$vars_to_keep
    file_name <- table_info[[dep_var]]$file_name
    result_names <- c()

    # Iterate over subsets, covariate sets, fixed effects, and treatment variables # nolint
    for (sub in subsets) {
      for (cov_set in cov_sets) {
        for (fe in fe_list) {
          for (treat_var in treat_vars) {
            # Build the result name
            result_name <- paste0(
              sub, "_", dep_var, "_", cov_set, "_", fe, "_", gsub(" ", "_", treat_var) # nolint
            )
            # Check if result exists
            if (result_name %in% names(results)) {
              result_names <- c(result_names, result_name)
            }
          }
        }
      }
    }

    # Use etable to output the table
    if (length(result_names) > 0) {
      fixest::etable(
        results[result_names],
        keep = vars_to_keep,
        file = paste0(tables, file_name)
      )
    }
  }
}

# Generate tables
generate_tables(
  dep_vars = dep_vars,
  table_info = table_info,
  subsets = subsets,
  cov_sets = cov_sets,
  fe_list = fe_list,
)