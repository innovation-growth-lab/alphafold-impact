# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)
source("scripts/utils.R")
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
pathdir <- "data/05_model_output/authors/ecr/quarterly/"

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

field_cols <- grep("^field_", names(sub_samples[[1]]), value = TRUE)
mesh_cols <- grep("^mesh_", names(sub_samples[[1]]), value = TRUE)

covs <- list()
covs[["base0"]] <- c(
  field_cols,
  mesh_cols,
  "num_publications"
)

fes <- list()
fes[["fe1"]] <- c("author", "quarter")

cov_sets <- c("base0")
fe_list <- c("fe1")
dep_vars <- c(
  "num_publications",
  "cited_by_count",
  "mesh_C",
  "ln1p_fwci",
  "ln1p_resolution",
  "ln1p_R_free",
  "patent_count",
  "patent_citation",
  "num_pdb_ids",
  "num_pdb_submissions",
  "ca_count",
  "num_uniprot_structures",
  "num_primary_submissions",
  "num_diseases",
  "ln1p_organism_rarity_mean",
  "ln1p_organism_rarity_max",
  "ln1p_max_tmscore",
  "ln1p_max_fident",
  "ln1p_max_score",
  "normalised_max_tmscore",
  "normalised_max_fident",
  "normalised_max_score",
  "num_uniprot_structures_w_disease",
  "num_primary_submissions_w_disease",
  "num_uniprot_structures_w_rare_organisms",
  "num_primary_submissions_w_rare_organisms",
  "num_uniprot_structures_w_low_similarity",
  "num_primary_submissions_w_low_similarity",
  "ln1p_maxtmscore_lt_0.405"
)

# Define base treatment vars that exist in all samples
treat_vars_base <- paste(
  c(
    "af", "ct_ai", "ct_pp", "ct_sb"
  ),
  collapse = " + "
)


treat_vars_base_w_het <- paste(
  c(
    "af", "ct_ai", "ct_pp", "ct_sb",
    "af:is_applied", "ct_ai:is_applied", "ct_pp:is_applied", "ct_sb:is_applied"
  ),
  collapse = " + "
)

# Define treatment vars with strong interactions
treat_vars_with_strong <- paste(
  c(
    "af_intent_strong",
    "af_intent_weak",
    "af_intent_mixed",
    "ct_ai_intent_strong",
    "ct_ai_intent_weak",
    "ct_ai_intent_mixed",
    "ct_pp_intent_strong",
    "ct_pp_intent_weak",
    "ct_pp_intent_mixed",
    "ct_sb_intent_strong",
    "ct_sb_intent_weak",
    "ct_sb_intent_mixed"
  ),
  collapse = " + "
)

treat_vars_with_strong_w_het <- paste(
  c(
    "af_intent_strong",
    "af_intent_weak",
    "af_intent_strong:is_applied",
    "af_intent_weak:is_applied",
    "ct_ai_intent_strong",
    "ct_ai_intent_weak",
    "ct_ai_intent_strong:is_applied",
    "ct_ai_intent_weak:is_applied",
    "ct_pp_intent_strong",
    "ct_pp_intent_weak",
    "ct_pp_intent_strong:is_applied",
    "ct_pp_intent_weak:is_applied",
    "ct_sb_intent_strong",
    "ct_sb_intent_weak",
    "ct_sb_intent_strong:is_applied",
    "ct_sb_intent_weak:is_applied",
    "af_intent_mixed",
    "ct_ai_intent_mixed",
    "ct_pp_intent_mixed",
    "ct_sb_intent_mixed"
  ),
  collapse = " + "
)

for (dep_var in dep_vars) { # nolint
  form_list <- list()
  # Iterate over dependent variables
  # for (dep_var in dep_var_out) { # nolint
  # Iterate over covariate sets
  for (cov_set in cov_sets) {
    local_covs <- covs[[cov_set]]

    # if dep_var is num_publications, remove it from covs
    if (dep_var == "num_publications") {
      local_covs <- covs[[cov_set]][
        -which(covs[[cov_set]] == "num_publications")
      ]
    } else if (dep_var %in% c("ln1p_mesh_C", "mesh_C")) {
      local_covs <- covs[[cov_set]][-which(covs[[cov_set]] == "mesh_C")]
    }
    # Iterate over fixed effects
    for (fe in fe_list) {
      # Iterate over treatment variables
      for (local_treat_vars in c(
        treat_vars_base,
        treat_vars_base_w_het,
        treat_vars_with_strong,
        treat_vars_with_strong_w_het
      )) { # Create formula name using the subset and treatment vars
        form_name <- paste0(
          dep_var, "__", cov_set, "__", fe, "__",
          gsub(" ", "_", local_treat_vars)
        )

        # Create the appropriate formula
        if (length(local_covs) == 0) {
          form_list[[form_name]] <- as.formula(
            paste0(
              dep_var, " ~ ", local_treat_vars, " |",
              paste0(fes[[fe]], collapse = " + ")
            )
          )
        } else {
          form_list[[form_name]] <- as.formula(
            paste0(
              dep_var, " ~ ", local_treat_vars, " +",
              paste0(local_covs, collapse = " + "),
              "|", paste0(fes[[fe]], collapse = " + ")
            )
          )
        }
      }
      # }
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

      # if "intent" in form name, swap the "*_with_intent" vars
      if (grepl("scope_Intent", sub)) {
        # Get columns ending with _with_intent
        intent_cols <- names(local_data)[endsWith(
          names(local_data), "_with_intent"
        )]
        # For each intent column, overwrite the base column with its value
        for (col in intent_cols) {
          base_col <- gsub("_with_intent", "", col)
          local_data[[base_col]] <- local_data[[col]]
        }
        # Drop the _with_intent columns
        local_data <- local_data %>% select(-all_of(intent_cols))
      }

      # consider skipping regression if saturated
      dep_var <- strsplit(form, "__")[[1]][1]

      if (dep_var %in% c("patent_count", "patent_citation")) {
        local_data <- local_data[
          # did not update prior to 2021
          local_data$year >= 2021 & local_data$year <= 2025,
        ]
      }

      # remove na obs on dep_var
      non_na_data <- local_data[!is.na(local_data[[dep_var]]), ]

      # compute the unique number of quarter
      n_authors <- length(unique(non_na_data$author))
      n_quarters <- length(unique(non_na_data$quarter))

      if (
        n_authors + n_quarters
        > nrow(non_na_data)
      ) {
        message("Skipping regression. Not enough data.")
        results[[regression_label]] <- feols(
          as.formula(paste(dep_var, "~ 1")),
          data = local_data
        )
        next
      }

      # skipping regression if form includes "strong" but no strong var
      if (
        grepl("strong", form) && !("af_intent_strong" %in% names(local_data))
      ) {
        message(
          "Skipping regression. No strong intent data."
        )
        next
      }

      # run the regression as linear, but make exceptions for counts and binary variables
      # so actually once you drop enough, you can get a rough 25% increase, similar to the linear reg. #nolint
      # the main thing is, using ln is odd because it assumes continuous variables and far from zero values #nolint
      if (dep_var %in% c(
        "num_publications", "num_pdb_ids", "num_pdb_submissions",
        "ca_count", "patent_count", "patent_citation",
        "num_uniprot_structures", "cited_by_count",
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
            # Apply the collinearity fix function before running the regression
            local_data <- fix_perfect_collinearity(
              local_data, fes[["fe1"]], dep_var
            )

            model <- fepois(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter"),
              fixef.iter = 250,
              nthreads = 1,
              lean = FALSE,
              mem.clean = TRUE
            )

            # Check if model converged by looking at convergence code
            if (!model$convStatus) {
              message("Model did not converge, using fallback model")
              feols(as.formula(paste(dep_var, "~ 1")), data = local_data)
            } else {
              model
            }
          },
          error = function(e) {
            message(
              "Error in regression: ", regression_label, " - ", e$message
            )
            return(
              feols(as.formula(paste(dep_var, "~ 1")), data = local_data)
            )
          }
        )
      } else if (dep_var %in% c("ln1p_maxtmscore_lt_0.405", "mesh_C")) {
        message("Running Logistic regression")
        results[[regression_label]] <- tryCatch(
          {
            # Apply the collinearity fix function before running the regression
            local_data <- fix_perfect_collinearity(
              local_data, fes[["fe1"]], dep_var
            )

            model <- feglm(
              form_list[[form]],
              data = local_data,
              family = binomial(link = "logit"),
              cluster = c("author", "quarter"),
              fixef.iter = 250,
              nthreads = 1,
              lean = FALSE,
              mem.clean = TRUE
            )

            # Check if model converged by looking at convergence code
            if (!model$convStatus) {
              message("Model did not converge, using fallback model")
              feols(as.formula(paste(dep_var, "~ 1")), data = local_data)
            } else {
              model
            }
          },
          error = function(e) {
            message(
              "Error in regression: ", regression_label, " - ", e$message
            )
            return(
              feols(as.formula(paste(dep_var, "~ 1")), data = local_data)
            )
          }
        )
      } else {
        # run the regression
        results[[regression_label]] <- tryCatch(
          {
            feols(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter"),
              lean = FALSE,
              mem.clean = TRUE
            )
          },
          error = function(e) {
            message(
              "Error in regression: ", regression_label, " - ", e$message
            )
            return(
              feols(as.formula(paste(dep_var, "~ 1")), data = local_data)
            )
          }
        )
      }
    }
  }

  # ------------------------------------------------------------------------
  # GENERATE PLOTS
  # ------------------------------------------------------------------------

  # import from utils_figures.R
  source("scripts/authors/ecr/quarterly/utils_figures.R")
  message("Generating plots")
  tryCatch(
    {
      coef_table <- extract_coefficients(
        results = results,
        dep_vars = dep_var,
        subsets = names(sub_samples),
        cov_sets = cov_sets,
        fe_list = fe_list,
        treat_vars = c(treat_vars_base, treat_vars_with_strong),
        treat_var_interest = c(
          "af", "ct_ai", "ct_pp", "ct_sb",
          "af_intent_strong", "af_intent_weak",
          "ct_ai_intent_strong", "ct_ai_intent_weak",
          "ct_pp_intent_strong", "ct_pp_intent_weak",
          "ct_sb_intent_strong", "ct_sb_intent_weak",
          "af:is_applied", "ct_ai:is_applied", "ct_pp:is_applied",
          "ct_sb:is_applied", "af_intent_strong:is_applied",
          "af_intent_weak:is_applied", "ct_ai_intent_strong:is_applied",
          "ct_ai_intent_weak:is_applied", "ct_pp_intent_strong:is_applied",
          "ct_pp_intent_weak:is_applied", "ct_sb_intent_strong:is_applied",
          "ct_sb_intent_weak:is_applied"
        )
      )

      # save coef_table
      coefplot_dir <- paste0(pathdir, "coef_tables/")

      if (!dir.exists(coefplot_dir)) {
        dir.create(coefplot_dir, recursive = TRUE)
      }

      saveRDS(
        coef_table,
        paste0(coefplot_dir, dep_var, "_coef_table.rds")
      )
    },
    error = function(e) {
      message("Error in generating tables: ", e$message)
    }
  )
  # ------------------------------------------------------------------------
  # GENERATE TABLES
  # ------------------------------------------------------------------------

  # import from utils_tables.R
  source("scripts/authors/ecr/quarterly/utils_tables.R")
  message("Generating tables")
  tryCatch(
    {
      # Generate tables
      generate_tables(
        results = results,
        dep_vars = dep_var,
        table_info = table_info,
        subsets = names(sub_samples),
        cov_sets = cov_sets,
        fe_list = fe_list,
        treat_vars = c(treat_vars_base, treat_vars_with_strong),
        intermediate_path = "base"
      )
      generate_tables(
        results = results,
        dep_vars = dep_var,
        table_info = table_info,
        subsets = names(sub_samples),
        cov_sets = cov_sets,
        fe_list = fe_list,
        treat_vars = c(treat_vars_base_w_het, treat_vars_with_strong_w_het),
        intermediate_path = "base_w_het"
      )
    },
    error = function(e) {
      message("Error in generating tables: ", e$message)
    }
  )
}
