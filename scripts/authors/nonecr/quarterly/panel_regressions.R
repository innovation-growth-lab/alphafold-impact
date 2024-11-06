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
pathdir <- "data/05_model_output/authors/nonecr/quarterly/"

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
sub_samples <- readRDS(paste0(pathdir, "data/sub_samples.rds"))

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

covs <- list()
covs[["base0"]] <- c("num_publications")

fes <- list()
fes[["fe0"]] <- c("quarter")
fes[["fe1"]] <- c(
  "author", "quarter", "institution", "institution_type",
  "institution_country_code"
)

cov_sets <- c("base0")
fe_list <- c("fe1")
dep_vars <- c(
  "num_publications"
)

for (dep_var_out in dep_vars) { # nolint
  treat_vars <- c(
    "af_ind + ct_ai_ind + ct_noai_ind + af:ct_ai_ind + af:ct_noai_ind + strong_af_ind + strong_ct_ai_ind + strong_ct_noai_ind + strong_af:strong_ct_ai_ind + strong_af:strong_ct_noai_ind", # nolint
    "af + ct_ai + ct_noai + af:ct_ai + af:ct_noai + strong_af + strong_ct_ai + strong_ct_noai + strong_af:strong_ct_ai + strong_af:strong_ct_noai" # nolint
  )

  form_list <- list()
  # Iterate over dependent variables
  for (dep_var in dep_var_out) { # nolint
    # Iterate over covariate sets
    for (cov_set in cov_sets) {
      local_covs <- covs[[cov_set]]
      # if dep_var is num_publications, remove it from covs
      if (dep_var == "num_publications") {
        local_covs <- covs[[cov_set]][-which(covs[[cov_set]] == "num_publications")] # nolint
      } else {
        local_covs <- covs[[cov_set]]
      }
      # Iterate over fixed effects
      for (fe in fe_list) {
        # Iterate over treatment variables
        for (treat_var in treat_vars) {
          if (treat_var == "af_ind + ct_ai_ind + ct_noai_ind + af:ct_ai_ind + af:ct_noai_ind + strong_af_ind + strong_ct_ai_ind + strong_ct_noai_ind + strong_af:strong_ct_ai_ind + strong_af:strong_ct_noai_ind") { # nolint
            treat_var <- paste0("af + ct_ai + ct_noai + af:ct_ai + af:ct_noai + strong_af + strong_ct_ai + strong_ct_noai + strong_af:strong_ct_ai + strong_af:strong_ct_noai") # nolint
            label_var <- "af_ind + ct_ai_ind + ct_noai_ind + af:ct_ai_ind + af:ct_noai_ind + strong_af_ind + strong_ct_ai_ind + strong_ct_noai_ind + strong_af:strong_ct_ai_ind + strong_af:strong_ct_noai_ind" # nolint
          } else {
            label_var <- treat_var
          }
          # Check if covs[[cov_set]] is empty
          if (length(local_covs) == 0) {
            # Create formula without '+' before '|'
            form_list[[
              paste0(
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", label_var) # nolint
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
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", label_var) # nolint
              )
            ]] <- as.formula(
              paste0(
                dep_var, " ~ ", treat_var, " +",
                paste0(local_covs, collapse = " + "),
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

      # If form's string includes _ind_, then rename the columns
      if (grepl("_ind_", form)) {
        local_data <- local_data %>%
          select(-af, -ct_ai, -ct_noai) %>%
          rename(
            af = af_ind,
            ct_ai = ct_ai_ind,
            ct_noai = ct_noai_ind
          )
      }

      # consider skipping regression if saturated
      dep_var <- strsplit(form, "__")[[1]][1]

      non_na_data <- local_data[!is.na(local_data[[dep_var]]), ]

      # compute the unique number of quarter
      n_authors <- length(unique(non_na_data$author))
      n_quarters <- length(unique(non_na_data$quarter))
      n_institutions <- length(unique(non_na_data$institution))
      n_institution_types <- length(unique(non_na_data$institution_type))
      n_institution_countries <- length(
        unique(non_na_data$institution_country_code)
      )

      if (
        n_authors + n_quarters + n_institutions +
          n_institution_types + n_institution_countries
        > nrow(non_na_data)
      ) {
        message("Skipping regression: ", regression_label)
        results[[regression_label]] <- feols(
          as.formula(paste(dep_var, "~ 1")),
          data = local_data
        )
        next
      }

      # run the regression as linear, but make an exception for pdb_submission
      if (dep_var == "pdb_submission") {
        results[[regression_label]] <- tryCatch(
          {
            feols(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter"),
              family = binomial(link = "logit")
            )
          },
          error = function(e) {
            message("Error in regression: ", regression_label, " - ", e$message)
            return(feols(as.formula(paste(dep_var, "~ 1")), data = local_data))
          }
        )
      } else {
        # run the regression
        results[[regression_label]] <- tryCatch(
          {
            feols(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter")
            )
          },
          error = function(e) {
            message("Error in regression: ", regression_label, " - ", e$message)
            return(feols(as.formula(paste(dep_var, "~ 1")), data = local_data))
          }
        )
      }
    }
  }

  # ----------------------------------------------------------------------------
  # GENERATE TABLES
  # ----------------------------------------------------------------------------

  # import from utils_tables.R
  source("scripts/authors/nonecr/quarterly/utils_tables.R")
  message("Generating tables")
  tryCatch(
    {
      # Generate tables
      generate_tables(
        results = results,
        dep_vars = dep_var_out,
        table_info = table_info,
        subsets = names(sub_samples),
        treat_vars = treat_vars,
        cov_sets = cov_sets,
        fe_list = fe_list
      )
    },
    error = function(e) {
      message("Error in generating tables: ", e$message)
    }
  )

  # ----------------------------------------------------------------------------
  # GENERATE PLOTS
  # ----------------------------------------------------------------------------

  # import from utils_figures.R
  source("scripts/authors/nonecr/quarterly/utils_figures.R")
  message("Generating plots")
  tryCatch(
    {
      coef_table <- extract_coefficients(
        results = results,
        dep_vars = dep_var_out,
        subsets = names(sub_samples),
        cov_sets = cov_sets,
        fe_list = fe_list,
        treat_vars = treat_vars,
        treat_var_interest = c(
          "af", "af_ind", "ct_ai_ind", "ct_noai_ind", "ct_ai", "ct_noai",
          "af:ct_ai", "af:ct_noai",
          "strong_af", "strong_ct_ai", "strong_ct_noai",
          "strong_af:strong_ct_ai", "strong_af:strong_ct_noai",
          "strong_af_ind", "strong_ct_ai_ind", "strong_ct_noai_ind",
          "strong_af:strong_ct_ai_ind", "strong_af:strong_ct_noai_ind"
        )
      )

      generate_coef_plots(
        coef_table
      )
    },
    error = function(e) {
      message("Error in generating plots: ", e$message)
    }
  )
}
