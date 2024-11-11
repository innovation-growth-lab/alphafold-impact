# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)

# Check installation & load required packages
list_of_packages <- c(
  "tidyverse", "zoo", "fixest", "ggh4x", "stringr",
  "broom", "ggplot2", "tibble", "purrr"
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
pathdir <- "data/05_model_output/authors/nonecr/quarterly/es/"

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
sub_samples <- readRDS(
  "data/05_model_output/authors/nonecr/quarterly/data/es_sub_samples.rds"
)

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

lagsleads <- c(-9:-3, -1:9)

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
  "num_publications",
  "ln1p_cited_by_count", "ln1p_cit_0", "ln1p_cit_1",
  "ln1p_fwci", "logit_cit_norm_perc",
  "ln1p_patent_count", "ln1p_patent_citation", "ln1p_ca_count",
  "resolution", "R_free", "pdb_submission"
)

treat_vars <- list()
for (group in c(
  "rel_treat_af"
  # "rel_treat_strong_af", "rel_treat_ct_ai", "rel_treat_ct_noai",
  # "rel_treat_strong_ct_ai", "rel_treat_strong_ct_noai"
)) {
  treat_vars[[group]] <- c(paste0("`", group, "_", lagsleads, "`"))
}


results <- list()
for (dep_var_out in dep_vars) { # nolint

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
        for (treat_var in names(treat_vars)) {
          # Check if covs[[cov_set]] is empty
          # Create formula with '+' before '|'
          if (length(local_covs) == 0) {
            # Create formula without '+' before '|'
            form_list[[
              paste0(
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
              )
            ]] <- as.formula(
              paste0(
                dep_var, " ~ ", paste(
                  treat_vars[[treat_var]],
                  collapse = " + "
                ),
                "|", paste0(fes[[fe]], collapse = " + ")
              )
            )
          } else {
            form_list[[
              paste0(
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
              )
            ]] <- as.formula(
              paste0(
                dep_var, " ~ ", paste(
                  treat_vars[[treat_var]],
                  collapse = " + "
                ), " +",
                paste0(local_covs, collapse = " + "),
                "|", paste0(fes[[fe]], collapse = " + ")
              )
            )
          }
        }
      }
    }
  }

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
            placeholder_data <- data.frame(dep_var = c(0, 1))
            colnames(placeholder_data) <- dep_var
            return(feols(
              as.formula(paste(dep_var, "~ 1")),
              data = placeholder_data
            ))
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
}


plot_es <- function(model, yvar) {
  # Tidy the model and filter the terms
  tidy_model <- tidy(model, conf.int = TRUE) %>% # nolint
    filter(grepl("rel_treat", term)) %>% # nolint
    filter(!(grepl("rel_treat_9|`rel_treat_-9`", term))) %>%
    dplyr::mutate(t = c(-6:-3, -1:8)) %>% # nolint
    bind_rows(tibble(t = -2, estimate = 0, conf.low = 0, conf.high = 0)) %>% # nolint
    mutate(group = as.factor(case_when( # nolint
      t < 0 ~ 1,
      t >= 0 ~ 2
    )))

  # Calculate dynamic ylims based on confidence intervals
  ylims <- range(c(tidy_model$conf.low, tidy_model$conf.high), na.rm = TRUE)

  # Create the plot
  return(tidy_model %>% # nolint
    ggplot(aes(x = t, y = estimate)) + # nolint
    geom_hline(yintercept = 0, linetype = "longdash", color = "gray") + # nolint
    geom_vline(xintercept = 0, linetype = "longdash", color = "gray") + # nolint
    geom_rect( # nolint
      aes(
        xmin = -1.8, xmax = -1.2, ymin = -Inf, ymax = Inf
      ),
      fill = "gray", alpha = .1
    ) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), # nolint
      linetype = "solid", show.legend = FALSE,
      color = "darkgray", width = .1, linewidth = 1
    ) +
    geom_point(fill = "black", shape = 21, size = 2) + # nolint
    geom_line(linewidth = 1) + # nolint
    ggtitle("") + # nolint
    ylim(ylims) + # nolint
    labs(y = yvar, x = "Quarters Relative to First AF Citation") + # nolint
    scale_x_continuous(breaks = seq(-8, 8, by = 2)) + # nolint
    theme_classic() + # nolint
    theme( # nolint
      plot.title = element_text(hjust = 0.5), # nolint
      axis.title = element_text(size = 8)
    ))
}

figures <- "data/05_model_output/authors/nonecr/quarterly/figures/es/"

if (!dir.exists(figures)) {
  dir.create(figures, recursive = TRUE)
}

# # Iterate over all es_results and save a figure for each
for (result in names(results)) {
  tryCatch(
    {
      # split string into "__" separated parts
      parts <- strsplit(result, "__")[[1]]

      # dependent variable is fourth last
      dep_var <- parts[length(parts) - 3]

      model <- results[[result]]

      # Create the plot
      plot <- plot_es(model, dep_var)

      # Save the plot with a unique filename
      ggsave(
        paste0(figures, result, ".png"),
        plot,
        width = 12, height = 6
      )
    },
    error = function(e) {
      message(paste("Error in processing", result, ":", e$message))
    }
  )
}

# save the "results" object
saveRDS(
  results,
  "data/05_model_output/authors/nonecr/quarterly/data/es_results.rds"
)
