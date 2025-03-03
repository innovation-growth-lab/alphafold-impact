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

field_cols <- grep("^field_", names(sub_samples[[1]]), value = TRUE)

covs <- list()
covs[["base0"]] <- c(
  field_cols
)

fes <- list()
fes[["fe1"]] <- c("author", "quarter")

cov_sets <- c("base0")
fe_list <- c("fe1")
dep_vars <- c(
  "mesh_C",
  "num_publications",
  "ln1p_cited_by_count",
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
  "organism_rarity_mean",
  "mean_tmscore",
  "num_uniprot_structures_w_disease",
  "num_primary_submissions_w_disease",
  "num_uniprot_structures_w_rare_organisms",
  "num_primary_submissions_w_rare_organisms",
  "num_uniprot_structures_w_low_similarity",
  "num_primary_submissions_w_low_similarity"
)

# Define base treatment vars that exist in all samples
treat_vars_base <- paste(
  c(
    "af", "ct_ai", "ct_noai",
    "af:ct_ai", "af:ct_noai", "ct_ai:ct_noai",
    "af:ct_ai:ct_noai"
  ),
  collapse = " + "
)

# Define treatment vars with strong interactions
treat_vars_with_strong <- paste(
  c(
    "af_strong0",
    "af_strong1",
    "ct_ai_strong0",
    "ct_ai_strong1",
    "ct_noai_strong0",
    "ct_noai_strong1",
    "af_ct_ai_strong0",
    "af_ct_ai_strong1",
    "af_ct_noai_strong0",
    "af_ct_noai_strong1",
    "ct_ai_ct_noai_strong0",
    "ct_ai_ct_noai_strong1"
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
      local_covs <- covs[[cov_set]][-which(covs[[cov_set]] == "num_publications")] # nolint
    } else if (dep_var == "mesh_C") {
      local_covs <- covs[[cov_set]][-which(covs[[cov_set]] == "mesh_C")] # nolint
    } else {
      local_covs <- covs[[cov_set]]
    }
    # Iterate over fixed effects
    for (fe in fe_list) {
      # Iterate over treatment variables
      for (local_treat_vars in c(treat_vars_base, treat_vars_with_strong)) {
        # Create formula name using the subset and treatment vars
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

      # consider skipping regression if saturated
      dep_var <- strsplit(form, "__")[[1]][1]

      non_na_data <- local_data[!is.na(local_data[[dep_var]]), ]

      # compute the unique number of quarter
      n_authors <- length(unique(non_na_data$author))
      n_quarters <- length(unique(non_na_data$quarter))

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

      # skipping regression if form includes "strong" but no strong var
      if (grepl("strong", form) && !("af_strong0" %in% names(local_data))) {
        message("Skipping regression: ", regression_label)
        next
      }

      # run the regression as linear, but make exceptions for counts
      # so actually once you drop enough, you can get a rough 25% increase, similar to the linear reg. #nolint
      # the main thing is, using ln is odd because it assumes continuous variables and far from zero values #nolint
      if (dep_var %in% c(
        "num_publications", "num_pdb_ids", "num_pdb_submissions",
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
        message("Running Negative Binomial regression")
        results[[regression_label]] <- tryCatch(
          {
            # values for whome the fixed-effects-grouped data is not all 0
            local_data <- local_data %>%
              group_by(author) %>%
              filter(sum(!!sym(dep_var)) > 0) %>%
              ungroup() %>%
              group_by(quarter) %>%
              filter(sum(!!sym(dep_var)) > 0) %>%
              ungroup()

            # Apply the collinearity fix function before running the regression
            local_data <- fix_perfect_collinearity(
              local_data, fes[["fe1"]], dep_var
            )

            model <- fenegbin(
              form_list[[form]],
              data = local_data,
              cluster = c("author", "quarter_year"),
              fixef.iter = 100,
              nthreads = 1,
              lean = FALSE,
              mem.clean = TRUE
            )


            # Check if model converged by looking at convergence code
            if (!is.numeric(model$se[1])) {
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
  # GENERATE TABLES
  # ------------------------------------------------------------------------

  # import from utils_tables.R
  source("scripts/authors/nonecr/quarterly/utils_tables.R")
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
        treat_vars = c(treat_vars_base, treat_vars_with_strong)
      )
    },
    error = function(e) {
      message("Error in generating tables: ", e$message)
    }
  )

  # ------------------------------------------------------------------------
  # GENERATE PLOTS
  # ------------------------------------------------------------------------

  # import from utils_figures.R
  source("scripts/authors/nonecr/quarterly/utils_figures.R")
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
          "af", "ct_ai", "ct_noai",
          "af_strong0", "af_strong1",
          "ct_ai_strong0", "ct_ai_strong1",
          "ct_noai_strong0", "ct_noai_strong1"
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
}
