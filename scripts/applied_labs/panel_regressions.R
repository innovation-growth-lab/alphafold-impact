# %%
# Clean out the workspace
rm(list = ls())
options(max.print = 1000)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "fixest"
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

# Set working directory and paths
setwd("~/projects/alphafold-impact/")
figures <- "data/05_model_output/applied_labs/figures/"
tables <- "data/05_model_output/applied_labs/tables/"

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
# Load data
# ------------------------------------------------------------------------------
sub_samples <- readRDS("data/05_model_output/applied_sub_samples.rds")

# drop any column that starts with rel_
sub_samples <- lapply(
  sub_samples, function(x) x %>% select(-starts_with("rel_"))
)
# ------------------------------------------------------------------------------
# REGRESSIONS
# ------------------------------------------------------------------------------

field_cols <- grep("^field_", names(sub_samples$all), value = TRUE)
mesh_cols <- grep("^mesh_", names(sub_samples$all), value = TRUE)
institution_cols <- grep("^institution_", names(sub_samples$all), value = TRUE)

covs <- list()
covs[["base0"]] <- c("pdb_share_4q", "protein_share_4q", "experimental_share_4q") # nolint
covs[["base1"]] <- c(covs$base0, field_cols)
covs[["base2"]] <- c(covs$base1, mesh_cols)

fes <- list()
fes[["fe0"]] <- c("pi_id", "covid_share_2020")
fes[["fe1"]] <- c(fes$fe0, "time_qtly", "country", institution_cols)

### Extensive ###

# List of dependent variables
dep_vars <- c(
  "ca_count",
  "num_publications", "cit_0", "cit_1", "cited_by_count",
  "patent_count", "patent_citation",
  "resolution", "R_free"
)

# List of covariate sets
cov_sets <- c("base2")

# List of fixed effects
fe_list <- c("fe1")

# List of treatment variables
treat_vars <- c(
  "treatment_af_dyn",
  paste0(
    "treatment_af_dyn + treatment_af_dyn:intent_strong +",
    "treatment_af_dyn:ext_af"
  ),
  paste0(
    "treatment_af_dyn + treatment_af_dyn:intent_strong + ",
    "treatment_af_dyn:ext_af + treatment_af_dyn:intent_strong + ",
    "treatment_af_dyn:intent_strong:ext_af"
  )
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
      cluster = "pi_id"
    )
  }
}


# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------
# Define variables of interest
variable_interest <- c(
  "treatment_af_dyn",
  "treatment_af_dyn:intent_strong",
  "treatment_af_dyn:ext_af",
  "treatment_af_dyn:intent_strong:ext_af"
)

# Define mapping from dependent variables to variables of interest and file names # nolint
table_info <- list(
  "ca_count" = list(
    vars_to_keep = variable_interest,
    file_name = "01_ca_translational.tex"
  ),
  "tcc" = list(
    vars_to_keep = variable_interest,
    file_name = "02_tcc_translational.tex"
  ),
  "patent_count" = list(
    vars_to_keep = variable_interest,
    file_name = "03_patent_count_translational.tex"
  ),
  "patent_citation" = list(
    vars_to_keep = variable_interest,
    file_name = "04_patent_citation_translational.tex"
  ),
  "cit_0" = list(
    vars_to_keep = variable_interest,
    file_name = "05_cit_0_productivity.tex"
  ),
  "cit_1" = list(
    vars_to_keep = variable_interest,
    file_name = "06_cit_1_productivity.tex"
  ),
  "num_publications" = list(
    vars_to_keep = variable_interest,
    file_name = "07_num_publications_productivity.tex"
  ),
  "cited_by_count_std" = list(
    vars_to_keep = variable_interest,
    file_name = "08_cited_by_count_std_productivity.tex"
  ),
  "resolution" = list(
    vars_to_keep = variable_interest,
    file_name = "09_resolution.tex"
  ),
  "R_free" = list(
    vars_to_keep = variable_interest,
    file_name = "10_R_free.tex"
  )
)

# Define subsets, covariate sets, fixed effects, and treatment variables
subsets <- c(
  "all", "af_ct", "af_ct_ai", "af_ct_noai",
  "af_ct_w_high_pdb", "af_ct_cem"
)

# Function to generate tables
generate_tables <- function(dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint
  for (dep_var in dep_vars) {
    vars_to_keep <- table_info[[dep_var]]$vars_to_keep
    file_name <- table_info[[dep_var]]$file_name
    result_names <- c()
    col_names <- c()

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
              col_names <- c(col_names, sub)
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
        file = paste0(tables, file_name),
        headers = col_names
      )
    }
  }
}

generate_tables(
  dep_vars = c(
    "ca_count", "tcc",
    "num_publications", "cit_0", "cit_1", "cited_by_count_std",
    "patent_count", "patent_citation", "resolution", "R_free"
  ),
  table_info = table_info,
  subsets = subsets,
  cov_sets = cov_sets,
  fe_list = fe_list,
  treat_vars = treat_vars
)
