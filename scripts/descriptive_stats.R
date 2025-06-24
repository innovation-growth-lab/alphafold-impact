# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 300)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "aws.s3", "yaml", "zoo", "MatchIt"
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


# Assign commonly used dplyr functions
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows


# --------------------------- Paths & Inputs ---------------------------- #
pathdir <- "data/05_model_output/descriptive_stats/"

# Create directories if they do not exist
if (!dir.exists(pathdir)) {
  dir.create(pathdir, recursive = TRUE)
}

credentials <- yaml.load_file("conf/base/credentials.yml")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key, # nolint
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret, # nolint
  "AWS_DEFAULT_REGION" = "eu-west-2" # nolint
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
foundational_path <- "2025Q1/04_output/analysis/foundational_labs/individual/outputs_quarterly.parquet" # nolint
applied_path <- "2025Q1/04_output/analysis/applied_labs/individual/outputs_quarterly.parquet" # nolint
nonecr_path <- "2025Q1/03_primary/nonecr/publications_quarterly.parquet" # nolint # mistake saved to oct
ecr_path <- "2025Q1/03_primary/ecr/publications_quarterly.parquet" # nolint
papers_path <- "2025Q1/04_output/publications/regression_inputs.parquet" # nolint

# Fetch the data from the S3 bucket
foundational_labs_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = foundational_path,
  bucket = bucket
)

applied_labs_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = applied_path,
  bucket = bucket
)

nonecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = nonecr_path,
  bucket = bucket
)

ecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = ecr_path,
  bucket = bucket
)

papers_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = papers_path,
  bucket = bucket
)

# --------------------------- Transforms - Pre ---------------------------- #
# create column is_applied
foundational_labs_data <- foundational_labs_data %>%
  mutate(is_applied = 0)
applied_labs_data <- applied_labs_data %>%
  mutate(is_applied = 1)

# merge the data
labs_data <- bind_rows(foundational_labs_data, applied_labs_data)

# drop duplicates by author, quarter
labs_data <- labs_data %>%
  distinct(author, quarter, .keep_all = TRUE)

# drop obs after 2025Q1 (ie. 2025Q2)
labs_data <- labs_data %>%
  filter(quarter <= "2025Q1")

# --------------------------- Seek to merge --------------------- #
papers_data <- papers_data %>% mutate(
  author = row_number(),
  af = ifelse(source == "af", 1, 0),
  ct_ai = ifelse(source == "ct_ai", 1, 0),
  ct_pp = ifelse(source == "ct_pp", 1, 0),
  ct_sb = ifelse(source == "ct_sb", 1, 0),
  other = ifelse(source == "other", 1, 0),
  dataset = "papers",
  num_publications = 1,
  depth = as.factor(if_else(level == 0, "foundational", "applied")),
  is_applied = if_else(depth == "applied", 1, 0),
  num_pdb_submissions = ifelse(pdb_submission > 0, 1, 0),
  quarter = paste0(year(publication_date), "Q", quarter(publication_date))
)

papers_data <- papers_data %>%
  arrange(af, ct_ai, ct_pp, ct_sb, other) %>%
  distinct(id, .keep_all = TRUE)

nonecr_data <- nonecr_data %>%
  mutate(
    af = af_strong + af_weak + af_mixed + af_unknown,
    ct_ai = ct_ai_strong + ct_ai_weak + ct_ai_mixed + ct_ai_unknown,
    ct_pp = ct_pp_strong + ct_pp_weak + ct_pp_mixed + ct_pp_unknown,
    ct_sb = ct_sb_strong + ct_sb_weak + ct_sb_mixed + ct_sb_unknown,
    other = other_strong + other_weak + other_mixed + other_unknown,
    dataset = "nonecr",
    is_applied = if_else(depth == "applied", 1, 0),
    num_pdb_submissions = ifelse(is.na(pdb_submission), 0, pdb_submission)
  )

ecr_data <- ecr_data %>%
  mutate(
    af = af_strong + af_weak + af_mixed + af_unknown,
    ct_ai = ct_ai_strong + ct_ai_weak + ct_ai_mixed + ct_ai_unknown,
    ct_pp = ct_pp_strong + ct_pp_weak + ct_pp_mixed + ct_pp_unknown,
    ct_sb = ct_sb_strong + ct_sb_weak + ct_sb_mixed + ct_sb_unknown,
    other = other_strong + other_weak + other_mixed + other_unknown,
    dataset = "ecr",
    is_applied = if_else(depth == "applied", 1, 0),
    num_pdb_submissions = ifelse(is.na(pdb_submission), 0, pdb_submission)
  )

labs_data <- labs_data %>%
  mutate(
    af = af_strong + af_weak + af_mixed + af_unknown,
    ct_ai = ct_ai_strong + ct_ai_weak + ct_ai_mixed + ct_ai_unknown,
    ct_pp = ct_pp_strong + ct_pp_weak + ct_pp_mixed + ct_pp_unknown,
    ct_sb = ct_sb_strong + ct_sb_weak + ct_sb_mixed + ct_sb_unknown,
    other = other_strong + other_weak + other_mixed + other_unknown,
    dataset = "labs",
    num_pdb_submissions = ifelse(is.na(pdb_submission), 0, pdb_submission)
  )

# Select relevant columns for analysis
relevant_cols <- c(
  "author", "dataset", "is_applied", "quarter",
  "af", "ct_ai", "ct_pp", "ct_sb", "other",
  "num_pdb_submissions", "num_publications", "ca_count", "patent_count"
)

# filter data
papers_data <- papers_data %>% select(all_of(relevant_cols))
ecr_data <- ecr_data %>% select(all_of(relevant_cols))
nonecr_data <- nonecr_data %>% select(all_of(relevant_cols))
labs_data <- labs_data %>% select(all_of(relevant_cols))

# Append dataframes and select only relevant columns
raw_data <- papers_data %>%
  rbind(ecr_data) %>%
  rbind(nonecr_data) %>%
  rbind(labs_data)

# Add time period classification
raw_data <- raw_data %>%
  mutate(
    time_period = ifelse(quarter <= "2021Q2", "Before/Including 2021Q2", "After 2021Q2")
  )

# Determine ever-treated status per author
raw_data <- raw_data %>%
  group_by(author, dataset) %>%
  mutate(
    group = case_when(
      any(af > 0, na.rm = TRUE) ~ "AlphaFold Papers",
      any(ct_ai > 0, na.rm = TRUE) ~ "AI-intensive Frontiers",
      any(ct_pp > 0, na.rm = TRUE) ~ "Protein Prediction Frontiers",
      any(ct_sb > 0, na.rm = TRUE) ~ "Other Struct. Biol. Frontiers",
      TRUE ~ "Other Struct. Biol. Research"
    )
  ) %>%
  ungroup()

# ------------------------ Create Summary Tables ---------------------- #

# Load additional packages for table creation
if (!require(knitr)) install.packages("knitr")
if (!require(kableExtra)) install.packages("kableExtra")
library(knitr)
library(kableExtra)

# Function to create summary table for a given dataset
create_summary_table <- function(data, dataset_name) {
  # Filter for the specific dataset
  df <- data %>% filter(dataset == dataset_name)

  # Define technology groups and their labels
  tech_groups <- c(
    "AlphaFold Papers", "AI-intensive Frontiers", "Protein Prediction Frontiers",
    "Other Struct. Biol. Frontiers", "Other Struct. Biol. Research"
  )
  tech_labels <- c("AF", "CT-AI", "CT-PP", "CT-SB", "Other-SB")

  # Define metrics and time periods
  metrics <- c("num_publications", "num_pdb_submissions", "ca_count", "patent_count")
  metric_labels <- c("Publications", "PDB Submissions", "Clinical Citations", "Patent Citations")
  time_periods <- c("Before/Including 2021Q2", "After 2021Q2")

  # Create summary matrix with 8 rows (4 metrics x 2 time periods)
  results <- matrix(0, nrow = 8, ncol = 11) # Added one column for Period
  colnames(results) <- c("2021Q2", paste(rep(tech_labels, each = 2), rep(c("Adjacent", "Downstream"), 5), sep = "_"))

  # Create row names without time periods and add Period column
  row_names <- character(8)
  period_col <- character(8)
  for (i in seq_along(metric_labels)) {
    row_names[2 * i - 1] <- metric_labels[i]
    row_names[2 * i] <- metric_labels[i]
    period_col[2 * i - 1] <- "Pre" # <= 2021Q2 (6 chars max)
    period_col[2 * i] <- "Post" # > 2021Q2 (6 chars max)
  }
  rownames(results) <- row_names
  results[, "2021Q2"] <- period_col

  # Fill the matrix
  for (i in seq_along(tech_groups)) {
    tech_group <- tech_groups[i]
    tech_label <- tech_labels[i]

    # Get data for this technology group
    group_data <- df %>% filter(group == tech_group)

    if (nrow(group_data) > 0) {
      for (j in seq_along(metrics)) {
        metric <- metrics[j]

        # For each time period
        for (k in seq_along(time_periods)) {
          time_period <- time_periods[k]
          row_idx <- 2 * j - 2 + k # Calculate row index

          # Adjacent (foundational, is_applied = 0)
          adjacent_data <- group_data %>%
            filter(is_applied == 0, time_period == !!time_period)
          if (nrow(adjacent_data) > 0) {
            results[row_idx, paste0(tech_label, "_Adjacent")] <-
              format(sum(adjacent_data[[metric]], na.rm = TRUE), big.mark = ",")
          }

          # Downstream (applied, is_applied = 1)
          downstream_data <- group_data %>%
            filter(is_applied == 1, time_period == !!time_period)
          if (nrow(downstream_data) > 0) {
            results[row_idx, paste0(tech_label, "_Downstream")] <-
              format(sum(downstream_data[[metric]], na.rm = TRUE), big.mark = ",")
          }
        }
      }
    }
  }

  return(results)
}

# Function to write LaTeX table
write_latex_table <- function(mat, filename, caption) {
  # Create column headers including Period column
  col_headers <- c("2021Q2", rep(c("Adjacent", "Downstream"), 5))

  # Create the table
  latex_table <- kable(mat,
    format = "latex",
    booktabs = TRUE,
    align = paste0("cl", paste(rep("c", 10), collapse = "")), # center align Period, left align row names, center for data
    col.names = col_headers,
    caption = caption
  ) %>%
    add_header_above(c(" " = 1, "Pre" = 1, "AF" = 2, "CT-AI" = 2, "CT-PP" = 2, "CT-SB" = 2, "Other-SB" = 2)) %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    # Add some visual separation between time periods for each metric
    row_spec(c(2, 4, 6), extra_latex_after = "\\midrule")

  # Write to file
  writeLines(latex_table, filename)
  message("Table written to: ", filename)
}

# --------------------------- Generate Tables ---------------------------- #

message("Creating summary tables...")

# Get unique datasets
datasets <- unique(raw_data$dataset)

# Create tables for each dataset
for (dataset_name in datasets) {
  message("Processing ", dataset_name, "...")

  # Create summary table
  summary_mat <- create_summary_table(raw_data, dataset_name)

  # Write LaTeX table
  filename <- file.path(pathdir, paste0("table_", dataset_name, "_quarterly_time_split.tex"))
  caption <- paste("Quarterly counts by technology, depth, and time period -", str_to_title(dataset_name))
  write_latex_table(summary_mat, filename, caption)

  # Print summary statistics
  dataset_data <- raw_data %>% filter(dataset == dataset_name)
  message("Summary for ", dataset_name, ":")
  message("  Total observations: ", nrow(dataset_data))
  message("  Adjacent (foundational): ", sum(dataset_data$is_applied == 0, na.rm = TRUE))
  message("  Downstream (applied): ", sum(dataset_data$is_applied == 1, na.rm = TRUE))

  # Print time period breakdown
  time_summary <- dataset_data %>%
    count(time_period) %>%
    arrange(time_period)
  message("  Time period breakdown:")
  for (i in 1:nrow(time_summary)) {
    message("    ", time_summary$time_period[i], ": ", time_summary$n[i])
  }

  # Print technology group breakdown by time period
  group_summary <- dataset_data %>%
    count(group, is_applied, time_period) %>%
    pivot_wider(names_from = c(is_applied, time_period), values_from = n, values_fill = 0)
  print(group_summary)
  message("")
}

# --------------------------- Overall Summary ---------------------------- #

message("Overall summary across all datasets:")
overall_summary <- raw_data %>%
  group_by(dataset, group, is_applied, time_period) %>%
  summarise(
    total_pubs = sum(num_publications, na.rm = TRUE),
    total_pdb = sum(num_pdb_submissions, na.rm = TRUE),
    total_clinical = sum(ca_count, na.rm = TRUE),
    total_patents = sum(patent_count, na.rm = TRUE),
    .groups = "drop"
  )

print(overall_summary)

message("Analysis complete! Time-split LaTeX tables saved to: ", pathdir)
