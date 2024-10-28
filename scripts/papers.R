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
figures <- "data/05_model_output/figures/papers/"
tables <- "data/05_model_output/tables/papers/"

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

# for duplicate ids, keep first
papers <- papers %>% distinct(id, .keep_all = TRUE)

# Create 'strong' variable
papers <- papers %>%
  mutate(
    strong = if_else(chain_label %in% c("strong", "partial_strong"), 1, 0),
    strong = replace_na(strong, 0),
    strong = as.factor(strong)
  )

# Create treatment variables
papers <- papers %>%
  mutate(
    treatment_af_dyn = if_else(source == "af", 1, 0),
    treatment_ct_dyn = if_else(source == "ct", 1, 0)
  )

# Normalize 'group_pdb_count' by 0-1 scaling
papers <- papers %>%
  mutate(group_pdb_count = group_pdb_count / max(group_pdb_count, na.rm = TRUE))

# Create 'high_pdb' variable
papers <- papers %>%
  mutate(high_pdb = if_else(group_pdb_count >= quantile(group_pdb_count, 0.75, na.rm = TRUE), 1, 0)) # nolint

# Extract quarter-year from 'publication_date'
papers <- papers %>%
  mutate(
    publication_date = as.Date(publication_date),
    time_qtly = (year(publication_date) - 2018) * 4 + quarter(publication_date),
    time_qtly = as.factor(time_qtly)
  )

# Fill NA values in 'field_' and 'mesh_' prefix columns with 0
papers <- papers %>%
  mutate(across(starts_with("field_"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("mesh_"), ~ replace_na(., 0)))

# Create log-transformed variables
papers <- papers %>%
  mutate(
    cited_by_count_ln = log1p(cited_by_count),
    patent_count_ln = log1p(patent_count),
    patent_citation_ln = log1p(patent_citation),
    ca_count_ln = log1p(ca_count)
  )

# ------------------------------------------------------------------------------
# SUBSET
# ------------------------------------------------------------------------------

papers_lvl0 <- papers %>% filter(level == "0")
papers_lvl12 <- papers %>% filter(level != "0")

percentile_75_pdb <- quantile(papers$group_pdb_count, 0.75)

### SUBSETS ###
sub_samples <- list(
  "all_lvl0" = papers_lvl0,
  # papers either in ct or af
  "all_lvl0_ct" = papers_lvl0 %>% filter(str_detect(source, "af") | str_detect(source, "ct")), # nolint
  "all_lvl0_w_high_pdb" = papers_lvl0 %>% filter(group_pdb_count >= percentile_75_pdb), # nolint
  "af_ct_lvl0_w_high_pdb" = papers_lvl0 %>% filter(group_pdb_count >= percentile_75_pdb) %>% filter(str_detect(source, "af") | str_detect(source, "ct")), # nolint
  "all_lvl12" = papers_lvl12,
  "all_lvl12_ct" = papers_lvl12 %>% filter(str_detect(source, "af") | str_detect(source, "ct")), # nolint
  "all_lvl12_w_high_pdb" = papers_lvl12 %>% filter(group_pdb_count >= percentile_75_pdb), # nolint
  "af_ct_lvl12_w_high_pdb" = papers_lvl12 %>% filter(group_pdb_count >= percentile_75_pdb) %>% filter(str_detect(source, "af") | str_detect(source, "ct")) # nolint
)

# ------------------------------------------------------------------------------
# REGRESSIONS
# ------------------------------------------------------------------------------

field_cols <- grep("^field_", names(papers), value = TRUE)
mesh_cols <- grep("^mesh_", names(papers), value = TRUE)

covs <- list()
covs[["base1"]] <- c(field_cols)
covs[["base2"]] <- c(covs$base1, mesh_cols)

fes <- list()
fes[["fe1"]] <- c("time_qtly", "group_pdb_count")

### Extensive ###

# List of dependent variables
dep_vars <- c(
  "cited_by_count_ln", "ca_count_ln", "patent_count_ln", "patent_citation_ln"
)

# List of covariate sets
cov_sets <- c("base2")

# List of fixed effects
fe_list <- c("fe1")

# List of treatment variables
treat_vars <- c(
  "treatment_af_dyn + treatment_af_dyn:strong + treatment_af_dyn:high_pdb + treatment_af_dyn:strong:high_pdb" # nolint
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

# Define variables of interest
variable_interest <- c(
  "treatment_af_dyn",
  "treatment_af_dyn:strong",
  "treatment_af_dyn:high_pdb",
  "treatment_af_dyn:strong:high_pdb"
)

# Define mapping from dependent variables to variables of interest and file names # nolint
table_info <- list(
  "cited_by_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "03_cited_by_count_ln_productivity_lvl0.tex"
  ),
  "patent_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "01_patent_count_ln_productivity_lvl0.tex"
  ),
  "patent_citation_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "02_patent_citation_ln_productivity_lvl0.tex"
  ),
  "ca_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "04_ca_count_ln_productivity_lvl0.tex"
  )
)

# Define subsets, covariate sets, fixed effects
subsets_lvl0 <- c(
  "all_lvl0", "all_lvl0_ct", "all_lvl0_w_high_pdb", "af_ct_lvl0_w_high_pdb"
)

# Define mapping from dependent variables to variables
table_info_lvl12 <- list(
  "patent_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "01_patent_count_ln_productivity_lvl12.tex"
  ),
  "cited_by_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "03_cited_by_count_ln_productivity_lvl12.tex"
  ),
  "patent_citation_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "02_patent_citation_ln_productivity_lvl12.tex"
  ),
  "ca_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "04_ca_count_ln_productivity_lvl12.tex"
  )
)

# Define subsets, covariate sets, fixed effects
subsets_lvl12 <- c(
  "all_lvl12", "all_lvl12_ct", "all_lvl12_w_high_pdb", "af_ct_lvl12_w_high_pdb"
)

cov_sets <- c("base2")
fe_list <- c("fe1")
treat_vars <- c(
  "treatment_af_dyn",
  "treatment_af_dyn + treatment_af_dyn:strong + treatment_af_dyn:high_pdb + treatment_af_dyn:strong:high_pdb" # nolint
)

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

# Generate tables for lvl0
generate_tables(
  dep_vars = c(
    "cited_by_count_ln", "patent_count_ln", "patent_citation_ln", "ca_count_ln"
  ),
  table_info = table_info,
  subsets = subsets_lvl0,
  cov_sets = cov_sets,
  fe_list = fe_list,
  treat_vars = treat_vars
)

# Generate tables for lvl12
generate_tables(
  dep_vars = c(
    "cited_by_count_ln", "patent_count_ln", "patent_citation_ln", "ca_count_ln"
  ),
  table_info = table_info_lvl12,
  subsets = subsets_lvl12,
  cov_sets = cov_sets,
  fe_list = fe_list,
  treat_vars = treat_vars
)
