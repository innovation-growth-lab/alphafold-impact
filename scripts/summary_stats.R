# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)

# Check installation & LOAD PACKAGES
list_of_packages <- c(
  "arrow", "tidyverse", "ggplot2", "data.table", "bacondecomp",
  "fixest", "broom", "stargazer", "kableExtra", "patchwork",
  "extrafont", "RColorBrewer", "plotrix", "MatchIt", "cem", "scales"
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

nonecr <- readRDS(
  "data/05_model_output/authors/nonecr/quarterly/data/sub_samples.rds"
)[[1]]
ecr <- readRDS(
  "data/05_model_output/authors/ecr/quarterly/data/sub_samples.rds"
)[[1]]
labs <- readRDS(
  "data/05_model_output/labs/quarterly/data/sub_samples.rds"
)[[1]]

# break labs groups
labs <- labs %>%
  ungroup() %>%
  mutate(
    quarter = quarter_year,
    author = pi_id,
    field_biochemist = field_biochemistry_genetics_and_molecular_biology
  )
# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

# Define the quarters of interest for each dataset
quarters_nonecr <- 192:195
quarters_ecr <- 200:204
quarters_labs <- -9:-6

# Define the variables to consider
variables <- c(
  "num_publications",
  "ln1p_cited_by_count",
  "ln1p_fwci",
  "ca_count",
  "patent_count",
  "num_pdb_submissions",
  # "field_biochemist",
  # "field_medicine",
  "institution_h_index",
  "institution_2yr_mean_citedness"
)
# Function to get the last observation of each variable for each author
get_last_observation <- function(data) {
  data %>% # nolint
    group_by(author) %>% # nolint
    filter(row_number() == n()) %>% # nolint
    select(author, af, ct_ai, ct_noai) %>% # nolint
    ungroup() # nolint
}

# Get the last observation of af, ct_ai, and ct_noai for each author
last_obs_nonecr <- get_last_observation(nonecr)
last_obs_ecr <- get_last_observation(ecr)
last_obs_labs <- get_last_observation(labs)

# Function to assign treatment groups
assign_treatment_group <- function(data, last_obs) {
  data <- data %>% # nolint
    left_join(last_obs, by = "author", suffix = c("", "_last")) %>% # nolint
    mutate( # nolint
      treatment_group = case_when( # nolint
        af_last > 0 ~ "AF",
        ct_ai_last > 0 ~ "CT AI",
        ct_noai_last > 0 ~ "CT NOAI",
        TRUE ~ "Baseline Control"
      )
    )
  return(data)
}

# Assign treatment groups for each dataset
nonecr <- assign_treatment_group(nonecr, last_obs_nonecr)
ecr <- assign_treatment_group(ecr, last_obs_ecr)
labs <- assign_treatment_group(labs, last_obs_labs)

# Function to summarize data for specified quarters and variables
summarise_data <- function(data, quarters, variables) {
  data %>% # nolint
    filter(quarter %in% quarters) %>% # nolint
    mutate(across(all_of(variables), as.numeric, .names = "numeric_{col}")) %>% # nolint
    select(starts_with("numeric_"), treatment_group) %>% # nolint
    group_by(treatment_group) %>% # nolint
    summarise( # nolint
      across(
        starts_with("numeric_"), # nolint
        list(
          mean = \(x) mean(x, na.rm = TRUE),
          sd = \(x) sd(x, na.rm = TRUE)
        )
      )
    )
}

summary_nonecr <- summarise_data(nonecr, quarters_nonecr, variables)

summary_ecr <- summarise_data(ecr, quarters_ecr, variables)

summary_labs <- summarise_data(labs, quarters_labs, variables)

# combine the summaries
summary_combined <- bind_rows(
  summary_nonecr %>% mutate(dataset = "Established"),
  summary_ecr %>% mutate(dataset = "Early Career"),
  summary_labs %>% mutate(dataset = "Laboratories")
)

summary_pivoted <- summary_combined %>%
  pivot_longer(
    cols = starts_with("numeric_"),
    names_to = c("variable", ".value"),
    names_pattern = "numeric_(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = treatment_group,
    values_from = c(mean, sd),
    names_glue = "{.value}_{treatment_group}"
  )

summary_pivoted <- summary_pivoted %>%
  mutate(
    diff_AF_CT_AI = mean_AF - `mean_CT AI`,
    diff_AF_CT_NOAI = mean_AF - `mean_CT NOAI`,
    diff_AF_Baseline_Control = mean_AF - `mean_Baseline Control`
  )

# Create a single column that averages the differences
summary_pivoted <- summary_pivoted %>%
  mutate(
    avg_diff_AF = round(rowMeans(select(., diff_AF_CT_AI, diff_AF_CT_NOAI, diff_AF_Baseline_Control), na.rm = TRUE), 2)
  ) %>%
  select(-diff_AF_CT_AI, -diff_AF_CT_NOAI, -diff_AF_Baseline_Control)


# order columns to be mean, sd for each treatment group
summary_pivoted <- summary_pivoted %>%
  select(
    dataset,
    variable,
    mean_AF,
    sd_AF,
    "mean_CT AI",
    "sd_CT AI",
    "mean_CT NOAI",
    "sd_CT NOAI",
    "mean_Baseline Control",
    "sd_Baseline Control",
    avg_diff_AF
  )

# ------------------------------------------------------------------------------
# Addibng the number of observations
# ------------------------------------------------------------------------------

count_unique_authors <- function(data) {
  data %>% # nolint
    group_by(treatment_group, depth) %>% # nolint
    summarise(n_obs = n_distinct(author)) # nolint
}

n_obs_nonecr <- count_unique_authors(nonecr)
n_obs_ecr <- count_unique_authors(ecr)
n_obs_labs <- count_unique_authors(labs)

n_obs_combined <- bind_rows(
  n_obs_nonecr %>% mutate(dataset = "Established"),
  n_obs_ecr %>% mutate(dataset = "Early Career"),
  n_obs_labs %>% mutate(dataset = "Laboratories")
)

# ------------------------------------------------------------------------------
#  Create the final table
# ------------------------------------------------------------------------------

# name the Variable strings with more reasonable names
summary_pivoted <- summary_pivoted %>%
  mutate(
    variable = case_when(
      variable == "num_publications" ~ "Number of Publications",
      variable == "ln1p_cited_by_count" ~ "Log Citations",
      variable == "ln1p_fwci" ~ "Log Field-Weighted Citation Index",
      variable == "ca_count" ~ "Number of Clinical Citations",
      variable == "patent_count" ~ "Number of Patents",
      variable == "num_pdb_submissions" ~ "Number of PDB Submissions",
      variable == "institution_h_index" ~ "Institution h-index",
      variable == "institution_2yr_mean_citedness" ~ "Institution 2-Year Mean Citedness", # nolint
    ),
  ) %>%
  mutate_at(vars(mean_AF:"avg_diff_AF"), ~ round(., 2))

# Remove the dataset column
summary_pivotedend <- summary_pivoted %>%
  select(-dataset)

# Change column names to have two rows
colnames(summary_pivotedend) <- c(
  "Variable",
  "AF_Mean", "AF_SD",
  "CT AI_Mean", "CT AI_SD",
  "CT NOAI_Mean", "CT NOAI_SD",
  "Baseline Control_Mean", "Baseline Control_SD",
  "Avg Diff AF"
)

column_width <- "1.8cm"

# Create a LaTeX table using kable
table <- kable(
  summary_pivotedend, "latex", 
align = "lccccccccc", 
booktabs = TRUE, 
linesep = ""
) %>%
  kable_styling(latex_options = c("striped")) %>%
  add_header_above(
    c(" ", "Mean", "SE", "Mean", "SE", "Mean", "SE", "Mean", "SE", "")
  ) %>%
  add_header_above(
    c(
      "\\\\textit{Variable}" = 1,
      "Alphafold" = 2,
      "Counterfactual (AI)" = 2,
      "Counterfactual (no AI)" = 2,
      "Baseline Structural Biology" = 2,
      "Mean Diff" = 1
    ),
    escape = FALSE
  ) %>%
  column_spec(2, width = column_width, bold = TRUE) %>%
  column_spec(3, width = column_width) %>%
  column_spec(4, width = column_width, bold = TRUE) %>%
  column_spec(5, width = column_width) %>%
  column_spec(6, width = column_width, bold = TRUE) %>%
  column_spec(7, width = column_width) %>%
  column_spec(8, width = column_width, bold = TRUE) %>%
  column_spec(9, width = column_width) %>%
  column_spec(10, width = column_width)


# Save the table to a LaTeX file
save_kable(table, file = "data/05_model_output/summary/stats.tex")
