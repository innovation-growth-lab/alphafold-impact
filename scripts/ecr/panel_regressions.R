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
# LOAD DATA
# ------------------------------------------------------------------------------
sub_samples <- readRDS(paste0(pathdir, "sub_samples.rds"))
# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

covs <- list()
covs[["base0"]] <- c("num_publications")

fes <- list()
fes[["fe0"]] <- c("quarter_year")
fes[["fe1"]] <- c(
  "quarter_year", "institution", "institution_type",
  "institution_country_code"
)

cov_sets <- c("base0")
fe_list <- c("fe1")
dep_vars <- c(
  "ln1p_cited_by_count", "ln1p_cit_0", "ln1p_cit_1",
  "ln1p_fwci", "ln1p_cit_norm_perc",
  "ln1p_patent_count", "ln1p_patent_citation", "ln1p_ca_count",
  "resolution", "R_free"
)
treat_vars <- c(
  "af_ind + ct_ind", # because subgroups do not intersect (ie. ct_ai subgroups require ct_noai == 0) # nolint
  "af + ct + af^2 + ct^2 + af:ct + af^2:ct^2"
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
        if (treat_var == "af_ind + ct_ind") {
          treat_var <- paste0("af + ct")
          label_var <- "af_ind + ct_ind"
        } else {
          label_var <- treat_var
        }
        # Check if covs[[cov_set]] is empty
        if (length(covs[[cov_set]]) == 0) {
          # Create formula without '+' before '|'
          form_list[[
            paste0(
              dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", label_var)
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
              dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", label_var)
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
    regression_label <- paste0(sub, "__", form)
    message("Running regression: ", regression_label)


    # Create a local copy of the subset
    local_data <- sub_samples[[sub]]

    # If form's string includes _ind_, then rename the columns in the local copy
    if (grepl("_ind_", form)) {
      local_data <- local_data %>%
        select(-af, -ct) %>%
        rename(
          af = af_ind,
          ct = ct_ind
        )
    }

    results[[regression_label]] <- tryCatch(
      {
        feols(
          form_list[[form]],
          data = local_data,
          cluster = "author"
        )
      },
      error = function(e) {
        message("Error in regression: ", regression_label, " - ", e$message)
        return(NULL) # Return NULL if an error occurs
      }
    )
  }
}

# ------------------------------------------------------------------------------
# GENERATE TABLES
# ------------------------------------------------------------------------------

# import from utils_tables.R
source("scripts/ecr/utils_tables.R")

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


# ------------------------------------------------------------------------------
# GENERATE PLOTS
# ------------------------------------------------------------------------------

# import from utils_figures.R
source("scripts/ecr/utils_figures.R")

coef_table <- extract_coefficients(
  results = results,
  dep_vars = dep_vars,
  subsets = names(sub_samples),
  cov_sets = cov_sets,
  fe_list = fe_list,
  treat_vars = treat_vars,
  treat_var_interest = c(
    "is_afTRUE", "is_afTRUE:af", "is_afTRUE:I(af^2)", "is_ct_aiTRUE",
    "is_ct_aiTRUE:ct_ai", "is_ct_aiTRUE:I(ct_ai^2)", "is_ct_noaiTRUE",
    "is_ct_noaiTRUE:ct_noai", "is_ct_noaiTRUE:I(ct_noai^2)", "is_af_ct_aiTRUE",
    "is_af_ct_aiTRUE:is_af_ct_ai", "is_af_ct_aiTRUE:I(is_af_ct_ai^2)",
    "is_af_ct_noaiTRUE", "is_af_ct_noaiTRUE:is_af_ct_noai",
    "is_af_ct_noaiTRUE:I(is_af_ct_noai^2)"
  )
)

coef_order <- c(
  "AlphaFold + Counterfactual non-AI (ext.)^2",
  "AlphaFold + Counterfactual non-AI (ext.)",
  "AlphaFold + Counterfactual non-AI (int.)",
  "AlphaFold + Counterfactual AI (ext.)^2",
  "AlphaFold + Counterfactual AI (ext.)",
  "AlphaFold + Counterfactual AI (int.)",
  "Counterfactual non-AI (ext.)^2",
  "Counterfactual non-AI (ext.)",
  "Counterfactual non-AI (int.)",
  "Counterfactual AI (ext.)^2",
  "Counterfactual AI (ext.)",
  "Counterfactual AI (int.)",
  "AlphaFold (ext.)^2",
  "AlphaFold (ext.)",
  "AlphaFold (int.)"
)

# Change coefficient names to more readable labels
coef_labels <- c(
  "is_afTRUE" = "AlphaFold (int.)",
  "is_afTRUE:af" = "AlphaFold (ext.)",
  "is_afTRUE:I(af^2)" = "AlphaFold (ext.)^2",
  "is_ct_aiTRUE" = "Counterfactual AI (int.)",
  "is_ct_aiTRUE:ct_ai" = "Counterfactual AI (ext.)",
  "is_ct_aiTRUE:I(ct_ai^2)" = "Counterfactual AI (ext.)^2",
  "is_ct_noaiTRUE" = "Counterfactual non-AI (int.)",
  "is_ct_noaiTRUE:ct_noai" = "Counterfactual non-AI (ext.)",
  "is_ct_noaiTRUE:I(ct_noai^2)" = "Counterfactual non-AI (ext.)^2",
  "is_af_ct_aiTRUE" = "AlphaFold + Counterfactual AI (int.)",
  "is_af_ct_aiTRUE:af_ct_ai" = "AlphaFold + Counterfactual AI (ext.)",
  "is_af_ct_aiTRUE:I(af_ct_ai^2)" = "AlphaFold + Counterfactual AI (ext.)^2",
  "is_af_ct_noaiTRUE" = "AlphaFold + Counterfactual non-AI (int.)",
  "is_af_ct_noaiTRUE:af_ct_noai" = "AlphaFold + Counterfactual non-AI (ext.)", # nolint
  "is_af_ct_noaiTRUE:I(af_ct_noai^2)" = "AlphaFold + Counterfactual non-AI (ext.)^2" # nolint
)

# results with interactions
coef_table_interacted <- coef_table[
  coef_table$indep_vars != "is_af_+_is_ct_ai_+_is_ct_noai_+_is_af_ct_ai_+_is_af_ct_noai", # nolint
]

generate_coef_plots(
  coef_table_interacted,
  coef_order,
  coef_labels,
  "extensive"
)

# results without interactions
coef_table_noninteracted <- coef_table[
  coef_table$indep_vars == "is_af_+_is_ct_ai_+_is_ct_noai_+_is_af_ct_ai_+_is_af_ct_noai", # nolint
]
generate_coef_plots(
  coef_table_noninteracted,
  coef_order,
  coef_labels,
  "intensive"
)
