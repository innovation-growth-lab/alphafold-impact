# Event Study Analysis for Multiple Overlapping Treatments
# Using Callaway and Sant'Anna (2021) difference-in-differences methodology

# Clean workspace and set options
rm(list = ls())
options(max.print = 1000, width = 250)

# Load required packages
packages <- c("tidyverse", "did", "ggplot2", "scales", "showtext")
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

# Set up fonts
font_add_google("Mulish", "muli")
showtext_auto()

# Set working directory and paths
setwd("~/projects/alphafold-impact/")

# Allow script to work for both labs and authors data
args <- commandArgs(trailingOnly = TRUE)
data_type <- ifelse(length(args) > 0, args[1], "labs") # Default to "labs"

if (data_type == "labs") {
  data_path <- "data/05_model_output/labs/quarterly/"
} else if (data_type == "authors") {
  data_path <- "data/05_model_output/authors/nonecr/quarterly/"
} else {
  stop("data_type must be 'labs' or 'authors'")
}

output_dir <- paste0(data_path, "did_plots/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Loading", data_type, "data...\n")
sub_samples <- readRDS(paste0(data_path, "data/sub_samples.rds"))

# Use the first dataset (All Fields, All PDB)
data <- sub_samples[[1]]

cat("Loaded data dimensions:\n")
print(dim(data))

# Treatment and outcome variables
treatments <- c("af")# , "ct_ai", "ct_pp", "ct_sb")

# Full list of dependent variables as defined in panel_regressions.R
panel_dep_vars <- c(
  "ln1p_fwci",
  "cited_by_count",
  "ln1p_mesh_C",
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
  "num_primary_submissions_w_low_similarity"
)

# Keep only those variables present in the loaded dataset to avoid errors
outcome_vars <- intersect(panel_dep_vars, names(data))

cat("Number of outcome variables to analyse:", length(outcome_vars), "\n")

# Treatment and outcome labels for plots
treatment_labels <- list(
  "af" = "AlphaFold Papers",
  "ct_ai" = "AI-intensive Frontiers",
  "ct_pp" = "Protein Prediction Frontiers",
  "ct_sb" = "Other Struct. Biol. Frontiers",
  "other" = "Other Struct. Biol. Research"
)

outcome_labels <- list(
  "cited_by_count" = "Citations",
  "ln1p_mesh_C" = "MeSH Disease Terms (ln)",
  "ln1p_fwci" = "Field-Weighted Citation Impact (ln)",
  "ln1p_resolution" = "Resolution (ln)",
  "ln1p_R_free" = "R-free Value (ln)",
  "patent_count" = "Patent-Paper Citations",
  "patent_citation" = "Patent-Paper Citations",
  "num_pdb_ids" = "PDB Structures",
  "num_pdb_submissions" = "PDB Submissions",
  "ca_count" = "Clinical Article Citations",
  "num_uniprot_structures" = "UniProt Structures",
  "num_primary_submissions" = "Primary Structure Submissions",
  "num_diseases" = "Disease Annotations",
  "ln1p_organism_rarity_mean" = "Mean Organism Rarity (ln)",
  "ln1p_organism_rarity_max" = "Max Organism Rarity (ln)",
  "ln1p_max_tmscore" = "Max TM-Score (ln)",
  "ln1p_max_fident" = "Max Fold Identity (ln)",
  "ln1p_max_score" = "Max Score (ln)",
  "normalised_max_tmscore" = "Normalized Max TM-Score",
  "normalised_max_fident" = "Normalized Max Fold Identity",
  "normalised_max_score" = "Normalized Max Score",
  "num_uniprot_structures_w_disease" = "UniProt Structures / Disease",
  "num_primary_submissions_w_disease" = "Primary Submissions / Disease",
  "num_uniprot_structures_w_rare_organisms" = "UniProt Structures / Rare Organisms",  # nolint
  "num_primary_submissions_w_rare_organisms" = "Primary Submissions / Rare Organisms", # nolint
  "num_uniprot_structures_w_low_similarity" = "UniProt Structures / Low Similarity", # nolint
  "num_primary_submissions_w_low_similarity" = "Primary Submissions / Low Similarity" # nolint
)
# Create binary treatment indicators and timing variables
for (treatment in treatments) {
  data[[paste0("ever_", treatment)]] <- data[[treatment]] > 0

  # Find first adoption quarter for each author
  first_adoption <- data %>%
    filter(!!sym(treatment) > 0) %>%
    group_by(author) %>%
    summarise(first_quarter = min(quarter), .groups = "drop") %>%
    rename(!!paste0("first_", treatment, "_quarter") := first_quarter)

  data <- data %>%
    left_join(first_adoption, by = "author")
}

cat("Available quarters:\n")
print(unique(data$quarter))

# Function to prepare data for DiD analysis
prepare_did_data <- function(data, treatment_type) {
  cat("\nPreparing data for treatment:", toupper(treatment_type), "\n")

  # Create treatment timing variable
  first_col <- paste0("first_", treatment_type, "_quarter")

  treatment_data <- data %>%  # nolint
    mutate(  # nolint
      id = as.numeric(as.factor(author)),  # nolint
      tname = (
        as.numeric(
          substr(quarter, 1, 4)  # nolint
        ) + (as.numeric(substr(quarter, 6, 6)) - 1) * 0.25
      ),
      gname = ifelse(is.na(!!sym(first_col)), 0,  # nolint
        as.numeric(substr(!!sym(first_col), 1, 4)) +
          (as.numeric(substr(!!sym(first_col), 6, 6)) - 1) * 0.25
      )
    ) %>%
    arrange(id, tname)  # nolint

  # Remove missing treatment timing
  treatment_data <- treatment_data[!is.na(treatment_data$gname), ]

  # Balance sample to avoid computational issues
  n_treated <- sum(treatment_data$gname > 0)  # nolint
  n_control <- sum(treatment_data$gname == 0)  # nolint

  return(treatment_data)
}

# Prepare data for each treatment
did_data_list <- list()
for (treatment in treatments) {
  did_data_list[[treatment]] <- prepare_did_data(data, treatment)
}

# Function to run DiD analysis
run_did_analysis <- function(data, outcome_var, treatment_type) {
  cat("\n--- DiD: ", toupper(treatment_type), " -> ", outcome_var, " ---\n")

  # Get field columns
  covars_cols <- c("field_biochemist", "field_medicine", "covid_share_2020")

  if (length(covars_cols) > 0) {
    xformla <- as.formula(paste("~", paste(covars_cols, collapse = " + ")))
  } else {
    xformla <- ~1
  }

  cat("Using", length(covars_cols), " controls\n")

  # Run att_gt
  att_result <- tryCatch(
    {
      clean_data <- data[!is.na(data[[outcome_var]]), ]

      # Remove rows with missing covariates
      if (!identical(xformla, ~1)) {
        covar_names <- all.vars(xformla)
        for (covar in covar_names) {
          if (covar %in% colnames(clean_data)) {
            clean_data <- clean_data[!is.na(clean_data[[covar]]), ]
          }
        }
      }

      cat("Clean data:", nrow(clean_data), "observations\n")

      att_gt(  # nolint
        yname = outcome_var,
        tname = "tname",
        idname = "id",
        gname = "gname",
        xformla = ~1,
        data = clean_data,
        control_group = "notyettreated",
        anticipation = 0,
        clustervars = "id",
        panel = TRUE,
        allow_unbalanced_panel = TRUE,
        est_method = "dr"
      )
    },
    error = function(e) {
      cat("Error:", e$message, "\n")
      return(NULL)
    }
  )

  if (is.null(att_result)) {
    cat("Analysis failed\n")
    return(NULL)
  }

  # Create dynamic effects plot

  dyn_result <- aggte(att_result, type = "dynamic", na.rm = TRUE)  # nolint

  # Calculate y-axis limits including confidence intervals
  egt_values <- dyn_result$egt
  att_values <- dyn_result$att.egt
  se_values <- dyn_result$se.egt
  crit_val <- dyn_result$crit.val.egt

  visible_mask <- egt_values >= -1.5 & egt_values <= 2
  visible_att <- att_values[visible_mask]
  visible_se <- se_values[visible_mask]

  visible_ci_lower <- visible_att - crit_val * visible_se
  visible_ci_upper <- visible_att + crit_val * visible_se

  if (length(visible_att) > 0 && any(!is.na(visible_att))) {
    y_min <- min(c(visible_att, visible_ci_lower), na.rm = TRUE)
    y_max <- max(c(visible_att, visible_ci_upper), na.rm = TRUE)
    y_range <- y_max - y_min
    y_padding <- max(y_range * 0.1, 0.02)
    y_limits <- c(y_min - y_padding, y_max + y_padding)
  } else {
    y_limits <- NULL
  }

  # Create plot
  treatment_label <- treatment_labels[[treatment_type]]
  outcome_label <- ifelse(outcome_var %in% names(outcome_labels),
    outcome_labels[[outcome_var]], outcome_var
  )

  dyn_plot <- ggdid(dyn_result) +  # nolint
    labs(  # nolint
      title = paste("Event Study:", treatment_label, "→", outcome_label),
      x = "Event Time (Years Relative to Treatment)",
      y = "Treatment Effect (ATT)"
    ) +
    coord_cartesian(  # nolint
      xlim = c(-1.55, 2.05),
      ylim = if (!is.null(y_limits)) y_limits else NULL
    ) +
    theme_classic() +  # nolint
    theme(  # nolint
      text = element_text(family = "muli", size = 56),  # nolint
      axis.text = element_text(size = 56),
      axis.title = element_text(size = 60),
      plot.title = element_text(size = 66, face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"),  # nolint
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey"),  # nolint
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # nolint
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")  # nolint
    )

  # check folder output_dir / treatment_type exists
  if (!dir.exists(paste0(output_dir, treatment_type))) {
    dir.create(paste0(output_dir, treatment_type), recursive = TRUE)
  }

  ggsave(  # nolint
    filename = paste0(
      output_dir, treatment_type, "/event_study_", outcome_var, ".png"
    ),
    plot = dyn_plot,
    width = 14, height = 8, dpi = 300
  )

  return(list(att_gt = att_result, dynamic = dyn_result, formula = xformla))
}

# Run DiD analysis for each treatment and outcome
did_results <- list()

for (treatment in treatments) {
  cat("\n", rep("=", 60), "\n")
  cat("ANALYZING TREATMENT:", toupper(treatment), "\n")
  cat(rep("=", 60), "\n")

  did_results[[treatment]] <- list()

  for (var in outcome_vars) {
    did_results[[treatment]][[var]] <- run_did_analysis(
      did_data_list[[treatment]], var, treatment
    )
  }
}

# Save results
saveRDS(did_results, paste0(output_dir, "did_results_multiple_treatments.rds"))
saveRDS(did_data_list, paste0(output_dir, "did_data_multiple_treatments.rds"))

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")
cat("Results saved to:", output_dir, "\n")
