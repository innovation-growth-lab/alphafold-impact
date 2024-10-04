# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
source("scripts/utils.R")

# Check installation & LOAD PACKAGES
list_of_packages <- c(
  "arrow", "tidyverse", "ggplot2", "data.table", "bacondecomp",
  "fixest", "broom", "stargazer", "kableExtra", "patchwork",
  "extrafont", "RColorBrewer", "plotrix", "MatchIt", "cem", "zoo"
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

setwd("~/projects/alphafold-impact/")

figures <- "data/05_model_output/figures/"
tables <- "data/05_model_output/tables/"
figs <- list()

getwd()

select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize

################################################################################
################################################################################
############################### DATA PREP ######################################
################################################################################
################################################################################

papers <- read_parquet(
  "data/04_outputs/descriptive/ecr_data.parquet"
)

# create factors from af, ct, other, is_ecr
papers <- papers %>%
  mutate(
    af = as.factor(af),
    ct = as.factor(ct),
    other = as.factor(other),
    is_ecr = as.factor(is_ecr)
  )

# create ln1p count
papers <- papers %>%
  mutate(
    ln1p_count = log1p(count)
  )

# create factors for institution, country_code, and year
papers <- papers %>%
  mutate(
    institution = as.factor(institution),
    country_code = as.factor(country_code),
    type = as.factor(type),
    year = as.factor(year)
  )

# create treatment var which is 1 if af and 2021 onwards and 0 otherwise
papers$treatment_af_dyn <- papers$af
sub_samples <- list(
  "all_ecr" = papers %>% filter(is_ecr == 1),
  "af_ct_ecr" = papers %>% filter(is_ecr == 1 & (af == 1 | ct == 1))
)

covs <- list()
covs[["base0"]] <- c("type") # "institution",

fes <- list()
fes[["fe0"]] <- c("year", "country_code")

cov_sets <- c("base0") # c("base0", "base1", "base2")
fe_list <- c("fe0")
dep_vars <- c("ln1p_count")
treat_vars <- c("treatment_af_dyn")


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
    results[[paste0(sub, "_", form)]] <- feols(
      form_list[[form]],
      data = sub_samples[[sub]]
    )
  }
}

variable_interest <- c(
  "treatment_af_dyn", "typeeducation"
)

etable(results[c("all_ecr_ln1p_count_base0_fe0_treatment_af_dyn", "af_ct_ecr_ln1p_count_base0_fe0_treatment_af_dyn")], keep = variable_interest, file = "data/05_model_output/tables/feols_results.tex")
