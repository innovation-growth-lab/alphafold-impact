rm(list = ls())
options(max.print = 1000, width = 250)

packages <- c("tidyverse", "ggplot2", "scales", "showtext", "zoo")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

# Google font used in ES plots
font_add_google("Mulish", "muli")
showtext_auto()

# --------------------------- Paths & Inputs ---------------------------- #
setwd("~/projects/alphafold-impact/")


args <- commandArgs(trailingOnly = TRUE)
data_type <- ifelse(length(args) > 0, args[1], "labs")
# optional second arg: grouping mode (ever | current). default is ever
group_mode <- ifelse(length(args) > 1, args[2], "ever")
# validate
if (!group_mode %in% c("ever", "current")) {
  stop("group_mode must be 'ever' or 'current'")
}

if (data_type == "labs") {
  data_path <- "data/05_model_output/labs/quarterly/"
} else if (data_type == "authors") {
  data_path <- "data/05_model_output/authors/nonecr/quarterly/"
} else {
  stop("data_type must be 'labs' or 'authors'")
}

output_dir <- paste0(data_path, "quarterly_trends/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Loading", data_type, "data...\n")
sub_samples <- readRDS(paste0(data_path, "data/sub_samples.rds"))

# Use the first dataset (All Fields, All PDB) – same as ES script
raw_data <- sub_samples[[1]]

cat("Loaded data dimensions:\n")
print(dim(raw_data))


# --------------------------- Add Treatment Groups --------------------- #

if (group_mode == "ever") {
  # Determine ever-treated status per author
  raw_data <- raw_data %>%
    group_by(author) %>%
    mutate(
      ever_af = any(af > 0, na.rm = TRUE),
      ever_ct_ai = any(ct_ai > 0, na.rm = TRUE),
      ever_ct_pp = any(ct_pp > 0, na.rm = TRUE),
      ever_ct_sb = any(ct_sb > 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      group = case_when(
        ever_af ~ "AlphaFold Papers",
        ever_ct_ai ~ "AI-intensive Frontiers",
        ever_ct_pp ~ "Protein Prediction Frontiers",
        ever_ct_sb ~ "Other Struct. Biol. Frontiers",
        TRUE ~ "Other Struct. Biol. Research"
      )
    )
} else {
  # current quarter status
  raw_data <- raw_data %>%
    mutate(
      group = case_when(
        af > 0 ~ "AlphaFold Papers",
        ct_ai > 0 ~ "AI-intensive Frontiers",
        ct_pp > 0 ~ "Protein Prediction Frontiers",
        ct_sb > 0 ~ "Other Struct. Biol. Frontiers",
        TRUE ~ "Other Struct. Biol. Research"
      )
    )
}

# set factor levels
raw_data <- raw_data %>%
  mutate(group = factor(group, levels = c(
    "AlphaFold Papers", "AI-intensive Frontiers", "Protein Prediction Frontiers",
    "Other Struct. Biol. Frontiers", "Other Struct. Biol. Research"
  )))

# --------------------------- Variables & Labels ------------------------ #

outcome_vars <- c(
  "ln1p_fwci", "num_publications",
  "ln1p_max_tmscore", "ln1p_max_fident",
  "num_pdb_submissions", "num_diseases",
  "ca_count", "patent_count"
)

outcome_labels <- list(
  "ln1p_fwci" = "Field-Weighted Citation Impact (ln)",
  "num_publications" = "Number of Publications",
  "ln1p_max_tmscore" = "Max TM-Score (ln)",
  "ln1p_max_fident" = "Max Fold Identity (ln)",
  "num_pdb_submissions" = "PDB Submissions",
  "num_diseases" = "Disease-relevant Structures",
  "ca_count" = "Clinical Article Citations",
  "patent_count" = "Patent-Paper Citations"
)

data_type_label <- ifelse(
  data_type == "labs", "Laboratories", "Established Authors"
)

# --------------------------- Helper Functions -------------------------- #

# Convert the quarter string (e.g. "2020Q1") into a numeric time value
# (year + (quarter-1)/4), compatible with the ES script's tname variable
quarter_to_num <- function(q) {
  yr <- as.numeric(substr(q, 1, 4))
  qtr <- as.numeric(substr(q, 6, 6))
  yr + (qtr - 1) / 4
}

# Clean & aggregate data for each variable – annual version
prep_trend_data <- function(df, var) {
  df %>% # nolint
    select(author, quarter, group, !!sym(var)) %>% # nolint
    filter(!is.na(.data[[var]])) %>% # nolint
    mutate(year = as.integer(substr(quarter, 1, 4))) %>% # nolint
    group_by(year, group) %>% # nolint
    summarise( # nolint
      mean_val = mean(.data[[var]], na.rm = TRUE),
      se_val = sd(.data[[var]], na.rm = TRUE) / sqrt(n()), # nolint
      .groups = "drop"
    ) %>%
    mutate(
      lower = mean_val - 1.96 * se_val, # nolint
      upper = mean_val + 1.96 * se_val, # nolint
      year_date = as.Date(paste0(year, "-01-01"))
    ) %>%
    arrange(group, year) # nolint
}

# Clean & aggregate data for each variable –4qa
prep_quarterly_trend_data <- function(df, var) {
  df %>% # nolint
    select(author, quarter, group, !!sym(var)) %>% # nolint
    filter(!is.na(.data[[var]])) %>% # nolint
    mutate( # nolint
      year = as.integer(substr(quarter, 1, 4)),
      qtr = as.integer(substr(quarter, 6, 6)),
      quarter_date = as.Date(
        paste0(year, "-", sprintf("%02d", (qtr - 1) * 3 + 1), "-01") # nolint
      )
    ) %>%
    group_by(quarter, quarter_date, group) %>% # nolint
    summarise( # nolint
      mean_val = mean(.data[[var]], na.rm = TRUE), # nolint
      se_val = sd(.data[[var]], na.rm = TRUE) / sqrt(n()), # nolint
      .groups = "drop"
    ) %>%
    arrange(group, quarter_date) %>% # nolint
    group_by(group) %>%
    mutate(
      # 4-quarter rolling average
      mean_val_smooth = rollmean(mean_val, k = 4, fill = NA, align = "right"), # nolint
      se_val_smooth = rollmean(se_val, k = 4, fill = NA, align = "right"), # nolint
      lower = mean_val_smooth - 1.96 * se_val_smooth, # nolint
      upper = mean_val_smooth + 1.96 * se_val_smooth # nolint
    ) %>%
    ungroup() %>% # nolint
    filter(!is.na(mean_val_smooth)) %>%
    arrange(group, quarter_date)
}

# --------------------------- Plotting Loop ----------------------------- #

for (var in outcome_vars) {
  cat("Creating trend plot for:", var, "...\n")

  trend_df <- prep_trend_data(raw_data, var) %>%
    filter(year >= 2016, year < 2025) # start from 2016

  # Determine y-axis label
  y_lab <- ifelse(var %in% names(outcome_labels),
    outcome_labels[[var]], var
  )

  palette <- c(
    "AlphaFold Papers" = "#D55E00",
    "AI-intensive Frontiers" = "#0072B2",
    "Protein Prediction Frontiers" = "#009E73",
    "Other Struct. Biol. Frontiers" = "#F0E442",
    "Other Struct. Biol. Research" = "#7B4FA3"
  )

  dodge <- position_dodge(width = 90) # 3 months separation for clarity

  plt <- ggplot(
    trend_df, aes(x = year_date, y = mean_val, colour = group)
  ) +
    geom_rect(
      aes(
        xmin = as.Date(-Inf), xmax = as.Date("2020-07-01"),
        ymin = -Inf, ymax = Inf
      ),
      fill = "#F0F0F0", alpha = 0.1, inherit.aes = FALSE
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 30,
      linewidth = 0.9,
      alpha = 0.8
    ) +
    geom_point(size = 3.2, position = dodge) +
    geom_line(linewidth = 1, alpha = 0.8, position = dodge) +
    scale_colour_manual(values = palette, name = "Group") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste("Annual Trend:", y_lab, "|", data_type_label),
      x = "Year",
      y = y_lab
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "muli", size = 46),
      axis.text = element_text(size = 46),
      axis.title = element_text(size = 50),
      legend.text = element_text(size = 42),
      legend.title = element_text(size = 46, face = "bold"),
      plot.title = element_text(size = 56, face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

  ggsave(
    filename = paste0(output_dir, "trend_annual_", var, ".png"),
    plot = plt,
    width = 8, height = 6, dpi = 300
  )
}

# ---------------------- Quarterly Plotting Loop ------------------------- #

for (var in outcome_vars) {
  cat("Creating quarterly trend plot for:", var, "...\n")

  quarterly_df <- prep_quarterly_trend_data(raw_data, var) %>%
    filter(
      quarter_date >= as.Date("2016-01-01"),
      quarter_date < as.Date("2025-01-01")
    )

  # Skip if no data after filtering
  if (nrow(quarterly_df) == 0) {
    cat("No data available for", var, "after filtering. Skipping...\n")
    next
  }

  # Determine y-axis label
  y_lab <- ifelse(var %in% names(outcome_labels),
    outcome_labels[[var]], var
  )

  palette <- c(
    "AlphaFold Papers" = "#D55E00",
    "AI-intensive Frontiers" = "#0072B2",
    "Protein Prediction Frontiers" = "#009E73",
    "Other Struct. Biol. Frontiers" = "#F0E442",
    "Other Struct. Biol. Research" = "#7B4FA3"
  )

  dodge <- position_dodge(width = 30) # 1 month separation for clarity

  plt_quarterly <- ggplot(
    quarterly_df, aes(x = quarter_date, y = mean_val_smooth, colour = group)
  ) +
    geom_rect(
      aes(
        xmin = as.Date(-Inf), xmax = as.Date("2020-07-01"),
        ymin = -Inf, ymax = Inf
      ),
      fill = "#F0F0F0", alpha = 0.1, inherit.aes = FALSE
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 15,
      linewidth = 0.9,
      alpha = 0.8
    ) +
    geom_point(size = 3.2, position = dodge) +
    geom_line(linewidth = 1, alpha = 0.8, position = dodge) +
    scale_colour_manual(values = palette, name = "Group") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste(
        "Quarterly Trend (4Q Rolling Avg):", y_lab, "|", data_type_label
      ),
      x = "Year",
      y = y_lab
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "muli", size = 46),
      axis.text = element_text(size = 46),
      axis.title = element_text(size = 50),
      legend.text = element_text(size = 42),
      legend.title = element_text(size = 46, face = "bold"),
      plot.title = element_text(size = 56, face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )

  ggsave(
    filename = paste0(output_dir, "trend_quarterly_", var, ".png"),
    plot = plt_quarterly,
    width = 14, height = 8, dpi = 300
  )
}

# ---------------------- Faceted Plot: All Variables -------------------- #

cat("Creating faceted plot for all variables...\n")

# Prepare combined data for all variables (using annual data)
facet_data_list <- list()

for (var in outcome_vars) {
  annual_df <- prep_trend_data(raw_data, var) %>%
    filter(year >= 2016, year < 2025)

  if (nrow(annual_df) > 0) {
    annual_df$variable <- outcome_labels[[var]]
    facet_data_list[[var]] <- annual_df
  }
}

# Combine all data
combined_facet_data <- bind_rows(facet_data_list)

# Set factor levels for variables to match outcome_vars order
variable_levels <- sapply(outcome_vars, function(var) outcome_labels[[var]])
combined_facet_data$variable <- factor(combined_facet_data$variable, levels = variable_levels)

# Create the faceted plot
if (nrow(combined_facet_data) > 0) {
  palette <- c(
    "AlphaFold Papers" = "#D55E00",
    "AI-intensive Frontiers" = "#0072B2",
    "Protein Prediction Frontiers" = "#009E73",
    "Other Struct. Biol. Frontiers" = "#F0E442",
    "Other Struct. Biol. Research" = "#7B4FA3"
  )

  dodge <- position_dodge(width = 90)

    plt_faceted <- ggplot(
    combined_facet_data, 
    aes(x = year_date, y = mean_val, colour = group)
  ) +
    geom_rect(
      aes(
        xmin = as.Date(-Inf), xmax = as.Date("2020-07-01"),
        ymin = -Inf, ymax = Inf
      ),
      fill = "#F0F0F0", alpha = 0.1, inherit.aes = FALSE
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 30,
      linewidth = 0.7,
      alpha = 0.7
    ) +
    geom_point(size = 2, position = dodge, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8, position = dodge) +
    scale_colour_manual(values = palette, name = "Group") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    labs(
      title = paste("Annual Trends |", data_type_label),
      x = "Year",
      y = "Value"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "muli", size = 32),
      axis.text = element_text(size = 28),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 36),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 32, face = "bold"),
      plot.title = element_text(size = 40, face = "bold"),
      strip.text = element_text(size = 30, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(colour = guide_legend(nrow = 3))

  ggsave(
    filename = paste0(output_dir, "trend_annual_faceted_all.png"),
    plot = plt_faceted,
    width = 8, height = 12, dpi = 300
  )

  cat("Faceted plot saved successfully.\n")
} else {
  cat("No data available for faceted plot.\n")
}

cat("\n", rep("=", 60), "\n")
cat("ANNUAL AND QUARTERLY TREND PLOTS COMPLETE\n")
cat(rep("=", 60), "\n")
cat("Plots saved to:", output_dir, "\n")
