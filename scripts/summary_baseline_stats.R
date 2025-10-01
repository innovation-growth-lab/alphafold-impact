# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 250)

# Check installation & LOAD PACKAGES
list_of_packages <- c(
  "arrow", "tidyverse", "ggplot2", "data.table", "bacondecomp",
  "fixest", "broom", "stargazer", "kableExtra", "patchwork",
  "extrafont", "RColorBrewer", "plotrix", "MatchIt", "cem", "scales",
  "effsize"
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

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

baseline_years <- 2015:2020

# Define the variables to consider
variables_base <- c(
  "num_publications",
  "cited_by_count",
  "ln1p_fwci",
  "ca_count",
  "patent_count",
  "num_pdb_submissions",
  "ln1p_max_tmscore",
  "institution_h_index",
  "institution_2yr_mean_citedness"
)

# papers data lack num_publications; we'll tailor per dataset later

# Function to get the last observation of each variable for each author
get_last_observation <- function(data) {
  data %>% # nolint
    group_by(author) %>% # nolint
    filter(row_number() == n()) %>% # nolint
    select(author, af, ct_ai, ct_pp, ct_sb) %>% # nolint
    ungroup() # nolint
}

# Get the last observation of af, ct_ai, ct_pp, ct_sb for each author
last_obs_nonecr <- get_last_observation(nonecr)
last_obs_labs <- get_last_observation(labs)
last_obs_ecr <- get_last_observation(ecr)

# Function to assign treatment groups
assign_treatment_group <- function(data, last_obs) {
  data <- data %>% # nolint
    left_join(last_obs, by = "author", suffix = c("", "_last")) %>% # nolint
    mutate( # nolint
      treatment_group = case_when( # nolint
        af_last > 0 ~ "AF",
        ct_ai_last > 0 ~ "AI Fr.",
        ct_pp_last > 0 ~ "PP Fr.",
        ct_sb_last > 0 ~ "SB Fr.",
        TRUE ~ "Baseline Control"
      )
    )
  return(data)
}

# Assign treatment groups for each dataset
nonecr <- assign_treatment_group(nonecr, last_obs_nonecr)
labs <- assign_treatment_group(labs, last_obs_labs)
ecr <- assign_treatment_group(ecr, last_obs_ecr)

# Function to summarize data for specified quarters and variables
summarise_data <- function(data, years, variables) {
  var_present <- intersect(variables, names(data))
  data %>% # nolint
    filter(year %in% years) %>% # nolint
    mutate(
      across(
        all_of(var_present),
        ~ {
          if (is.numeric(.)) {
            .
          } else {
            suppressWarnings(as.numeric(as.character(.)))
          }
        },
        .names = "numeric_{col}"
      )
    ) %>%
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

# Build variable list by dataset (keep only columns that exist in each data frame)
variables_ecr <- intersect(variables_base, names(ecr))
variables_nonecr <- intersect(variables_base, names(nonecr))
variables_labs <- intersect(variables_base, names(labs))

summary_nonecr <- summarise_data(nonecr, baseline_years, variables_nonecr)
summary_ecr <- summarise_data(ecr, baseline_years, variables_ecr)
summary_labs <- summarise_data(labs, baseline_years, variables_labs)

# Re-order: Publications first, then Established, then Laboratories
summary_combined <- bind_rows(
  summary_ecr %>% mutate(dataset = "Early Career"),
  summary_nonecr %>% mutate(dataset = "Established"),
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
    diff_AF_AI = mean_AF - `mean_AI Fr.`,
    diff_AF_PP = mean_AF - `mean_PP Fr.`,
    diff_AF_SB = mean_AF - `mean_SB Fr.`,
    diff_AF_Baseline_Control = mean_AF - `mean_Baseline Control`
  )

# Create a single column that averages the differences
summary_pivoted <- summary_pivoted %>%
  mutate(
    avg_diff_AF = round(rowMeans(select(., diff_AF_AI, diff_AF_PP, diff_AF_SB, diff_AF_Baseline_Control), na.rm = TRUE), 2)
  ) %>%
  select(-diff_AF_AI, -diff_AF_PP, -diff_AF_SB, -diff_AF_Baseline_Control)


# order columns to be mean, sd for each treatment group
summary_pivoted <- summary_pivoted %>%
  select(
    dataset,
    variable,
    mean_AF,
    sd_AF,
    "mean_AI Fr.",
    "sd_AI Fr.",
    "mean_PP Fr.",
    "sd_PP Fr.",
    "mean_SB Fr.",
    "sd_SB Fr.",
    "mean_Baseline Control",
    "sd_Baseline Control",
    avg_diff_AF
  )

# ------------------------------------------------------------------------------
# Addibng the number of observations
# ------------------------------------------------------------------------------

count_unique_authors <- function(data) {
  data %>% # nolint
    group_by(treatment_group) %>% # nolint
    summarise(n_obs = n_distinct(author)) # nolint
}

n_obs_nonecr <- count_unique_authors(nonecr)
n_obs_labs <- count_unique_authors(labs)
n_obs_ecr <- count_unique_authors(ecr)

# Same ordering for counts
n_obs_combined <- bind_rows(
  n_obs_ecr %>% mutate(dataset = "Early Career"),
  n_obs_nonecr %>% mutate(dataset = "Established"),
  n_obs_labs %>% mutate(dataset = "Laboratories")
)


# ---------------------------------------------------------------------------
# 1.  Reshape each dataset to long format separately
# ---------------------------------------------------------------------------
analysis_data_nonecr <- nonecr %>%
  mutate(dataset = "Established") %>%
  mutate(across(all_of(variables_nonecr), as.numeric)) %>%
  filter(year %in% baseline_years) %>% # keep baseline yrs
  select(dataset, author, treatment_group, all_of(variables_nonecr)) %>%
  pivot_longer(
    cols = all_of(variables_nonecr),
    names_to = "variable",
    values_to = "value"
  )

analysis_data_ecr <- ecr %>%
  mutate(dataset = "Early Career") %>%
  mutate(across(all_of(variables_ecr), as.numeric)) %>%
  filter(year %in% baseline_years) %>% # keep baseline yrs
  select(dataset, author, treatment_group, all_of(variables_ecr)) %>%
  pivot_longer(
    cols = all_of(variables_ecr),
    names_to = "variable",
    values_to = "value"  )

analysis_data_labs <- labs %>%
  mutate(dataset = "Laboratories") %>%
  mutate(across(all_of(variables_labs), as.numeric)) %>%
  filter(year %in% baseline_years) %>% # keep baseline yrs
  select(dataset, author, treatment_group, all_of(variables_labs)) %>%
  pivot_longer(
    cols = all_of(variables_labs),
    names_to = "variable",
    values_to = "value"
  )


# ---------------------------------------------------------------------------
# 2.  Function that returns p-values (Wilcoxon rank-sum test) + effect sizes
# ---------------------------------------------------------------------------
get_stats <- function(df) {
  af <- df$value[df$treatment_group == "AF"]
  baseline <- df$value[df$treatment_group == "Baseline Control"]

  # Remove missing values
  af <- af[!is.na(af)]
  baseline <- baseline[!is.na(baseline)]

  # Sample sizes
  n_af <- length(af)
  n_baseline <- length(baseline)

  # Skip if too few observations
  if (n_af < 10 || n_baseline < 10) {
    return(tibble(
      p_value = NA,
      n_af = n_af,
      n_baseline = n_baseline,
      cohens_d = NA,
      effect_size_interpretation = "Insufficient data"
    ))
  }

  # Wilcoxon test
  test <- wilcox.test(af, baseline)

  # Cohen's d (effect size)
  pooled_sd <- sqrt(((n_af - 1) * var(af) + (n_baseline - 1) * var(baseline)) /
    (n_af + n_baseline - 2))
  cohens_d <- (mean(af) - mean(baseline)) / pooled_sd

  # Effect size interpretation
  effect_interpretation <- case_when(
    abs(cohens_d) < 0.2 ~ "Negligible",
    abs(cohens_d) < 0.5 ~ "Small",
    abs(cohens_d) < 0.8 ~ "Medium",
    TRUE ~ "Large"
  )

  tibble(
    p_value = test$p.value,
    n_af = n_af,
    n_baseline = n_baseline,
    cohens_d = round(cohens_d, 3),
    effect_size_interpretation = effect_interpretation
  )
}

# Calculate stats for each dataset separately
stats_tbl_nonecr <- analysis_data_nonecr %>%
  group_by(variable) %>%
  group_modify(~ get_stats(.x)) %>%
  ungroup() %>%
  mutate(dataset = "Established")

stats_tbl_labs <- analysis_data_labs %>%
  group_by(variable) %>%
  group_modify(~ get_stats(.x)) %>%
  ungroup() %>%
  mutate(dataset = "Laboratories")

stats_tbl_ecr <- analysis_data_ecr %>%
  group_by(variable) %>%
  group_modify(~ get_stats(.x)) %>%
  ungroup() %>%
  mutate(dataset = "Early Career")

# Combine the stats tables
stats_tbl <- bind_rows(stats_tbl_ecr, stats_tbl_nonecr, stats_tbl_labs)

# ---------------------------------------------------------------------------
# 3.  Combine with the summary table you already built
# ---------------------------------------------------------------------------
summary_pivoted <- summary_pivoted %>%
  left_join(stats_tbl, by = c("variable", "dataset")) %>%
  mutate(
    # Create clearer significance indicators
    statistical_sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    practical_sig = case_when(
      abs(cohens_d) >= 0.8 ~ "Large",
      abs(cohens_d) >= 0.5 ~ "Medium",
      abs(cohens_d) >= 0.2 ~ "Small",
      TRUE ~ "Negligible"
    ),
    # Combined significance (both statistical AND practical)
    combined_sig = case_when(
      p_value < 0.001 & abs(cohens_d) >= 0.2 ~ "***",
      p_value < 0.01 & abs(cohens_d) >= 0.2 ~ "**",
      p_value < 0.05 & abs(cohens_d) >= 0.2 ~ "*",
      TRUE ~ ""
    ),
    # Format p-values more clearly
    p_value_display = case_when(
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01 ~ paste0("<", format(round(p_value, 3), nsmall = 3)),
      p_value < 0.05 ~ format(round(p_value, 3), nsmall = 3),
      p_value < 0.1 ~ format(round(p_value, 3), nsmall = 3),
      TRUE ~ format(round(p_value, 2), nsmall = 2)
    )
  )


# ------------------------------------------------------------------------------
#  Create the final table
# ------------------------------------------------------------------------------

# name the Variable strings with more reasonable names
summary_pivoted <- summary_pivoted %>%
  mutate(
    variable = case_when(
      variable == "num_publications" ~ "Number of Publications",
      variable == "cited_by_count" ~ "Citations per Publication",
      variable == "ln1p_fwci" ~ "Log Field-Weighted Citation Index",
      variable == "ca_count" ~ "Number of Clinical Citations",
      variable == "patent_count" ~ "Number of Patents",
      variable == "num_pdb_submissions" ~ "Number of PDB Submissions",
      variable == "institution_h_index" ~ "Institution h-index",
      variable == "institution_2yr_mean_citedness" ~ "Institution 2-Year Mean Citedness", # nolint
      variable == "num_uniprot_structures" ~ "Number of UniProt Structures",
      variable == "num_primary_submissions" ~ "Number of Primary Submissions",
      variable == "ln1p_max_tmscore" ~ "Log Max TM-Score",
      variable == "field_biochemist" ~ "Biochemistry",
      variable == "field_medicine" ~ "Medicine"
    ),
  ) %>%
  mutate_at(vars(mean_AF:"avg_diff_AF"), ~ round(., 2))

# Remove the dataset column and add just Cohen's d and practical significance
summary_pivotedend <- summary_pivoted %>%
  select(
    Variable = variable,
    AF_Mean = mean_AF, AF_SD = sd_AF,
    mean_AI = `mean_AI Fr.`, sd_AI = `sd_AI Fr.`,
    mean_PP = `mean_PP Fr.`, sd_PP = `sd_PP Fr.`,
    mean_SB = `mean_SB Fr.`, sd_SB = `sd_SB Fr.`,
    mean_Base = `mean_Baseline Control`, sd_Base = `sd_Baseline Control`,
    Avg_Diff_AF = avg_diff_AF,
    cohens_d,
    effect_size = practical_sig,
    p_value = p_value_display,
    combined_sig
  )

# Change column names
colnames(summary_pivotedend) <- c(
  "Variable",
  "AF_Mean", "AF_SD",
  "AI Fr._Mean", "AI Fr._SD",
  "PP Fr._Mean", "PP Fr._SD",
  "SB Fr._Mean", "SB Fr._SD",
  "Baseline_Mean", "Baseline_SD",
  "Mean Diff",
  "Cohen's d", "Effect Size",
  "p-value", "Practical Sig"
)

column_width <- "1.8cm"

# Create a LaTeX table using kable - back to simpler format
table <- kable(
  summary_pivotedend, "latex",
  align = "lcccccccccccccccc",
  booktabs = TRUE,
  linesep = ""
) %>%
  kable_styling(latex_options = c("striped")) %>%
  add_header_above(
    c(" ", "Mean", "SD", "Mean", "SD", "Mean", "SD", "Mean", "SD", "Mean", "SD", " ", " ", " ", " ", " ")
  ) %>%
  add_header_above(
    c(
      "\\\\textit{Variable}" = 1,
      "Alphafold" = 2,
      "AI Frontiers" = 2,
      "PP Frontiers" = 2,
      "SB Frontiers" = 2,
      "Baseline" = 2,
      "Effect Analysis" = 4
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
  column_spec(10, width = column_width, bold = TRUE) %>%
  column_spec(11, width = column_width) %>%
  column_spec(12, width = column_width) %>%
  column_spec(13, width = column_width, bold = TRUE) %>%
  column_spec(14, width = column_width) %>%
  column_spec(15, width = column_width, bold = TRUE) %>%
  column_spec(16, width = column_width, bold = TRUE)


# Save the table to a LaTeX file
save_kable(table, file = "data/05_model_output/summary/stats.tex")
