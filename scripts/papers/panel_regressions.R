# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)

# Check installation & load required packages
list_of_packages <- c(
  "tidyverse", "zoo", "fixest", "ggh4x"
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
# LOAD DATA
# ------------------------------------------------------------------------------
sub_samples <- readRDS(paste0(pathdir, "sub_samples.rds"))

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

field_cols <- grep("^field_", names(sub_samples$all_lvl0), value = TRUE)
mesh_cols <- grep("^mesh_", names(sub_samples$all_lvl0), value = TRUE)

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

source("scripts/papers/utils_tables.R")

# Generate tables
generate_tables(
  results = results,
  dep_vars = dep_vars,
  table_info = table_info,
  subsets = names(sub_samples),
  treat_vars = treat_vars,
  cov_sets = cov_sets,
  fe_list = fe_list
)
