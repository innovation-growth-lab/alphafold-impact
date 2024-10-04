# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
source("scripts/utils.R")

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "fixest", "stringr",
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

# Set working directory and paths
setwd("~/projects/alphafold-impact/")
figures <- "data/05_model_output/figures/"
tables <- "data/05_model_output/tables/"
figs <- list()


# Assign commonly used dplyr functions
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
sub_samples <- readRDS("data/05_model_output/applied_sub_samples.rds")

# ------------------------------------------------------------------------------
# Event Study Regressions
# ------------------------------------------------------------------------------

field_cols <- grep("^field_", names(sub_samples$all), value = TRUE)
mesh_cols <- grep("^mesh_", names(sub_samples$all), value = TRUE)
institution_cols <- grep("^institution_", names(sub_samples$all), value = TRUE)

# List of dependent variables
dep_vars <- c(
  "ca_count", # "tcc",
  "num_publications", "ct0", "ct1", "cited_by_count",
  "patent_count", "patent_citation"
)

# Define the lags and leads for the event study
lagsleads <- c(-9:-3, -1:9)
regs <- list()
for (group in c(
  "rel_treat_af", "rel_treat_af_strong"
)) {
  regs[[group]] <- c(paste0("`", group, "_", lagsleads, "`"))
}

# Baseline covariates
es_covs <- c("covid_share_2020", "high_pdb", field_cols, mesh_cols, institution_cols)

# # Add field and mesh columns
# interacted_covs <- expand.grid(
#   time = "time_qtly", cols = c(field_cols, mesh_cols)
# )
# interacted_covs <- apply(
#   interacted_covs, 1, function(x) paste(x, collapse = ":")
# )

# Add additional interactions
interacted_covs <- c(
  # interacted_covs,
  paste0("time_qtly:", c("protein_share", "experimental_share"))
)

# Final set of covariates
es_covs <- c(es_covs, interacted_covs)

# Add fixed effects
es_fes <- c("pi_id", "time_qtly", "time_qtly^country")

# Create formulas for each regs group
es_form_list <- list()
for (group in names(regs)) {
  for (dep_var in dep_vars) {
    es_form_list[[paste0(dep_var, "__", group)]] <- as.formula(
      paste0(
        dep_var, " ~ ",
        paste(regs[[group]], collapse = " + "), " + ",
        paste(es_covs, collapse = " + "), " | ",
        paste(es_fes, collapse = " + ")
      )
    )
  }
}

# Run the event study regressions systematically
es_results <- list()
for (sub in names(sub_samples)) {
  # skip if group has ai or noai
  if (str_detect(sub, "ai") || str_detect(sub, "noai")) {
    next
  }
  for (group in names(es_form_list)) {
    print(paste0(sub, "__", group))
    es_results[[paste0(sub, "__", group)]] <- feols(
      es_form_list[[group]],
      data = sub_samples[[sub]],
      cluster = "pi_id"
    )
  }
}

# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

plot_es <- function(model, yvar) {
  # Tidy the model and filter the terms
  tidy_model <- tidy(model, conf.int = TRUE) %>% # nolint
    filter(grepl("rel_treat", term)) %>% # nolint
    filter(!(grepl("rel_treat_9|`rel_treat_-9`", term))) %>%
    dplyr::mutate(t = c(-9:-3, -1:9)) %>% # nolint
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

# Iterate over all es_results and save a figure for each
for (result in names(es_results)) {
  tryCatch(
    {
      model <- es_results[[result]]

      # Extract the dependent variable
      dep_var <- str_extract(result, paste0(dep_vars, collapse = "|"))

      # Create the plot
      plot <- plot_es(model, dep_var)

      # Save the plot with a unique filename
      ggsave(
        paste0(figures, "exploration/applied/", result, ".png"),
        plot,
        width = 12, height = 6
      )
    },
    error = function(e) {
      message(paste("Error in processing", result, ":", e$message))
    }
  )
}
