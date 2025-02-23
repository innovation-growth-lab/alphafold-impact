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
sub_samples <- readRDS(paste0(pathdir, "data/sub_samples.rds"))

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

field_cols <- grep("^field_", names(sub_samples$all_lvl0), value = TRUE)
mesh_cols <- grep("^mesh_", names(sub_samples$all_lvl0), value = TRUE)

covs <- list()
covs[["base0"]] <- c(
  field_cols,
  "institution_type",
  "institution_cited_by_count",
  "institution_2yr_mean_citedness",
  "institution_h_index",
  "institution_i10_index",
  "institution_country_code"
)

fes <- list()
fes[["fe0"]] <- c("quarter_year")
fes[["fe1"]] <- c("author", "quarter_year")



### Extensive ###

cov_sets <- c("base0")
fe_list <- c("fe1")
dep_vars <- c(
  # "mesh_C"
  # "ln1p_cited_by_count",
  # "ln1p_fwci",
  # "ln1p_resolution",
  # "ln1p_R_free",
  # "patent_count",
  # "patent_citation",
  "num_pdb_ids"
  # "ca_count",
  # "num_uniprot_structures",
  # "num_primary_submissions",
  # "num_diseases",
  # "organism_rarity_mean",
  # "mean_tmscore",
  # "num_uniprot_structures_w_disease",
  # "num_primary_submissions_w_disease",
  # "num_uniprot_structures_w_rare_organisms",
  # "num_primary_submissions_w_rare_organisms",
  # "num_uniprot_structures_w_low_similarity",
  # "num_primary_submissions_w_low_similarity"

)
treat_vars <- paste(
  c(
    "af", "ct_ai", "ct_noai", "af:strong", "ct_ai:strong",
    "ct_noai:strong", "af:ct_ai", "af:ct_noai"
  ),
  collapse = " + "
)
for (dep_var_out in dep_vars) { # nolint
  form_list <- list()
  # Iterate over dependent variables
  for (dep_var in dep_vars) { # nolint
    # Iterate over covariate sets
    for (cov_set in cov_sets) {
      local_covs <- covs[[cov_set]]
      # Iterate over fixed effects
      for (fe in fe_list) {
        # Iterate over treatment variables
        for (treat_var in treat_vars) {
          # Check if covs[[cov_set]] is empty
          if (length(local_covs) == 0) {
            # Create formula without '+' before '|'
            form_list[[
              paste0(
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
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
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
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

      # consider skipping regression if saturated
      dep_var <- strsplit(form, "__")[[1]][1]

      non_na_data <- local_data[!is.na(local_data[[dep_var]]), ]

      # compute the unique number of quarter_year
      n_authors <- length(unique(non_na_data$author))
      n_quarters <- length(unique(non_na_data$quarter_year))

      if (
        n_authors + n_quarters
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
      if (dep_var %in% c(
        "num_publications", "num_pdb_submissions", "num_pdb_ids", 
        "ca_count", "patent_count", "patent_citation",
        "num_uniprot_structures",
        "num_primary_submissions",
        "num_diseases",
        "num_uniprot_structures_w_disease",
        "num_primary_submissions_w_disease",
        "num_uniprot_structures_w_rare_organisms",
        "num_primary_submissions_w_rare_organisms",
        "num_uniprot_structures_w_low_similarity",
        "num_primary_submissions_w_low_similarity"
      )) {
        message("Running Poisson regression")
        results[[regression_label]] <- tryCatch(
          {
            fepois(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter_year"),
              control = list(maxit = 500)
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
              cluster = c("author", "quarter_year")
            )
          },
          error = function(e) {
            message("Error in regression: ", regression_label, " - ", e$message)
            placeholder_data <- data.frame(dep_var = c(0, 1))
            colnames(placeholder_data) <- dep_var
            return(feols(
              as.formula(paste(dep_var, "~ 1")),
              data = placeholder_data
            ))
          }
        )
      }
    }
  }

  # ----------------------------------------------------------------------------
  # GENERATE TABLES
  # ----------------------------------------------------------------------------

  # import from utils_tables.R
  source("scripts/papers/utils_tables.R")
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
  source("scripts/papers/utils_figures.R")
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
          "af", "ct_ai", "ct_noai", "af:strong1",
          "ct_ai:strong1", "ct_noai:strong1"
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
