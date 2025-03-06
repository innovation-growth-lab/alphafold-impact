# Load necessary library
library(dplyr)
library(ggplot2)
library(showtext)

font_add_google("Mulish", "muli")
showtext_auto()

options(max.print = 1000)
options(width = 250)

# Define the paths containing the RDS files and their origins
base_paths <- list(
  nonecr = "data/05_model_output/authors/nonecr/quarterly/coef_tables/",
  ecr = "data/05_model_output/authors/ecr/quarterly/coef_tables/",
  labs = "data/05_model_output/labs/quarterly/coef_tables/",
  papers = "data/05_model_output/papers/coef_tables/" # nolint
)

# Function to collect and combine RDS files
combine_rds_files_with_source <- function(paths) {
  combined_table <- list() # Initialize list to hold data

  for (source_origin in names(paths)) {
    path <- paths[[source_origin]]

    # Get all RDS files in the directory
    rds_files <- list.files(
      path,
      pattern = "\\.rds$", full.names = TRUE, recursive = TRUE
    )

    # Load each RDS file and store in a list
    for (file in rds_files) {
      message("Processing file: ", file)
      table <- readRDS(file) # Load RDS file
      table$source_origin <- source_origin # Add origin source column
      combined_table[[file]] <- table
    }
  }

  # Combine all tables into one large table
  combined_table <- bind_rows(combined_table, .id = "source_file")

  return(combined_table)
}

# Run the function
all_tables_combined <- combine_rds_files_with_source(base_paths)

# if source_origin is papers, change pdb_submission to num_pdb_submissions
all_tables_combined <- all_tables_combined %>%
  mutate( # nolint
    treat_var = ifelse(
      source_origin == "papers" & treat_var == "pdb_submission",
      "num_pdb_submissions",
      treat_var
    )
  )

# count values of each source origin
all_tables_combined %>%
  count(source_origin) %>%
  print()

# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #
# ------------------------- #


dep_var_labels <- c(
  "ln1p_cited_by_count" = "ln (1 + Cited by count)",
  "ln1p_fwci" = "ln (1 + Field-Weighted Citation Impact)",
  "logit_cit_norm_perc" = "logit (Citation percentile)",
  "patent_count" = "Patent citation count",
  "patent_citation" = "Patent-to-patent citation count",
  "ca_count" = "Clinical article citation count",
  "ln1p_resolution" = "ln (1 + Resolution)",
  "ln1p_R_free" = "ln (1 + R-free)",
  "num_pdb_ids" = "PDB ID submissions",
  "num_pdb_submissions" = "PDB submissions",
  "num_publications" = "Publications",
  "num_uniprot_structures" = "Uniprot structures",
  "num_primary_submissions" = "Primary submissions",
  "num_diseases" = "Number of Diseases",
  "organism_rarity_mean" = "Organism rarity mean",
  "mean_tmscore" = "Mean TM-score",
  "num_uniprot_structures_w_disease" = "Uniprot structures with disease relevance", # nolint
  "num_primary_submissions_w_disease" = "Primary submissions with disease relevance", # nolint
  "num_uniprot_structures_w_rare_organisms" = "Uniprot structures with rare organisms", # nolint
  "num_primary_submissions_w_rare_organisms" = "Primary submissions with rare organisms", # nolint
  "num_uniprot_structures_w_low_similarity" = "Uniprot structures with low similarity", # nolint
  "num_primary_submissions_w_low_similarity" = "Primary submissions with low similarity" # nolint
)

# Updated coefficient labels with new naming
coef_labels <- c(
  # Core coefficients
  "af" = "AlphaFold",
  "ct_ai" = "AI Frontier",
  "ct_noai" = "No AI Frontier",

  # Strong coefficients
  "af_strong0" = "AlphaFold - Bkg.",
  "af_strong1" = "AlphaFold - Method",
  "ct_ai_strong0" = "AI Frontier - Bkg.",
  "ct_ai_strong1" = "AI Frontier - Method",
  "ct_noai_strong0" = "No AI Frontier - Bkg.",
  "ct_noai_strong1" = "No AI Frontier - Method"
)

strip_colors <- c(
  "Laboratories" = "#00B2A2",
  "Established Researchers" = "#FF5836",
  # "Early Career Researchers" = "#FAB61B", # nolint
  "Citation Papers" = "#1F5DAD"
)

field_colors <- c(
  "All Fields" = "black",
  "Medicine" = "#FF5836",
  "Molecular Biology" = "#1F5DAD"
)

# --- Variable definitions ---
# Set desired orders for variables and names
subgroup_order <- c("All PDB", "High PDB")
scope_order <- c("scope_All", "scope_Intent")
field_order <- c(
  "field_Medicine",
  "field_Molecular Biology",
  "field_All Fields"
)

# Updated observation order - removing ecr
core_observation_order <- c("papers", "nonecr", "labs")

observation_labels <- c(
  "papers" = "Cit Chains",
  "labs" = "Laboratories",
  "nonecr" = "Established",
  "ecr" = "Early Career Researchers"
)

# Define core coefficients and strong coefficients
core_coef_vars <- c("af", "ct_ai", "ct_noai")
strong_coef_vars <- c(
  "af_strong0", "af_strong1",
  "ct_ai_strong0", "ct_ai_strong1",
  "ct_noai_strong0", "ct_noai_strong1"
)

# --- Function to generate core coefficient plots ---
generate_core_coef_plots <- function(coef_table) {
  # Define variables locally to avoid "no visible global binding" warnings
  local_coef_labels <- coef_labels
  local_core_coef_vars <- core_coef_vars
  local_core_observation_order <- core_observation_order
  local_observation_labels <- observation_labels
  local_scope_order <- scope_order
  local_field_order <- field_order
  local_subgroup_order <- subgroup_order
  local_strip_colors <- strip_colors
  local_field_colors <- field_colors
  local_dep_var_labels <- dep_var_labels

  # rename vars
  coef_table <- coef_table %>%
    mutate(
      dep_var = recode(dep_var, !!!local_dep_var_labels) # nolint
    )

  unique_dep_vars <- unique(coef_table$dep_var)

  for (single_dep_var in unique_dep_vars) {
    # Filter for core coefficients only
    coef_plot_data <- coef_table %>%
      filter(
        treat_var %in% local_core_coef_vars, # nolint
        dep_var == single_dep_var, # nolint
        source_origin %in% local_core_observation_order # nolint
      ) %>%
      mutate(
        scope = factor(gsub("scope_", "", scope), # nolint
          levels = gsub("scope_", "", local_scope_order)
        ), # nolint
        field = factor(gsub("field_", "", field), # nolint
          levels = gsub("field_", "", local_field_order)
        ), # nolint
        subgroup = factor(gsub("subgroup_", "", subgroup), # nolint
          levels = gsub("subgroup_", "", local_subgroup_order)
        ), # nolint
        scope_subgroup = factor(paste(scope, subgroup, sep = " - "),
          levels = unique(paste(scope, subgroup, sep = " - "))
        ), # nolint
        source_origin = factor(source_origin, levels = local_core_observation_order), # nolint
        source_origin = recode(source_origin, !!!local_observation_labels), # nolint
        treat_var = factor(treat_var, levels = local_core_coef_vars), # nolint
        treat_var = recode(treat_var, !!!local_coef_labels) # nolint
      )

    coef_plot_data <- coef_plot_data %>%
      filter(!is.na(field) & !is.na(scope) & !is.na(subgroup)) # nolint

    for (scopesubgroup in unique(coef_plot_data$scope_subgroup)) { # nolint
      subgroup_coef_plot_data <- coef_plot_data %>%
        filter(scope_subgroup == scopesubgroup) # nolint

      # drop levels in column scope_subgroup if they have zero values
      vertical_levels <- levels(
        subgroup_coef_plot_data$scope_subgroup # nolint
      )[levels(subgroup_coef_plot_data$scope_subgroup) %in% # nolint
        unique(subgroup_coef_plot_data$scope_subgroup[
          subgroup_coef_plot_data$scope_subgroup != ""
        ])] # nolint
      vertical_strip_colors <- local_strip_colors[vertical_levels] # nolint

      horizontal_levels <- levels(
        subgroup_coef_plot_data$source_origin # nolint
      )[levels(subgroup_coef_plot_data$source_origin) %in% # nolint
        unique(subgroup_coef_plot_data$source_origin[
          subgroup_coef_plot_data$source_origin != ""
        ])] # nolint
      horizontal_strip_colors <- local_strip_colors[horizontal_levels] # nolint

      # Check if subgroup_coef_plot_data is empty
      if (nrow(subgroup_coef_plot_data) == 0) {
        message("No data for dep_var: ", single_dep_var)
        next
      }

      subgroup_coef_plot_data <- subgroup_coef_plot_data %>%
        filter(n_obs >= 200) # nolint

      # Create the plot
      coeffplot <- ggplot(
        subgroup_coef_plot_data,
        aes(x = estimate, y = treat_var, color = field) # nolint
      ) +
        geom_point(
          size = 4,
          position = position_dodge(width = 0.8)
        ) +
        geom_errorbarh(
          aes(
            xmin = estimate - 1.645 * std_error, # nolint
            xmax = estimate + 1.645 * std_error
          ),
          height = 0, linewidth = 1, # thicker for 10% significance
          position = position_dodge(width = 0.8)
        ) +
        geom_errorbarh(
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.2,
          linewidth = 0.5,
          position = position_dodge(width = 0.8) # thinner for 5% significance
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
        ggh4x::facet_grid2(
          . ~ source_origin,
          scales = "free",
          independent = "x",
          space = "fixed",
          strip = ggh4x::strip_themed(
            background_y = ggh4x::elem_list_rect(
              fill = vertical_strip_colors
            ),
            background_x = ggh4x::elem_list_rect(
              fill = horizontal_strip_colors
            )
          )
        ) +
        scale_color_manual(values = local_field_colors) +
        labs(
          title = paste(
            "Dependent Variable:", single_dep_var, " | ", scopesubgroup # # nolint
          ),
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() +
        theme(
          axis.text.y = element_text(size = 40),
          axis.title.x = element_text(size = 36),
          axis.title.y = element_text(size = 36),
          strip.text = element_text(size = 30),
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
          panel.spacing = unit(2, "lines"),
          plot.margin = margin(1, 1, 1, 1, "cm"),
          legend.position = "bottom",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30)
        )

      all_fields_data <- subgroup_coef_plot_data %>%
        filter(field == "All Fields") # nolint

      # add counts of obs (subgroup all and ext.)
      coeffplot <- coeffplot + geom_text(
        data = all_fields_data,
        aes(label = paste0("n = ", n_obs)), # nolint
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 40,
        size = 12, color = "black"
      ) + theme(
        text = element_text(family = "muli", size = 40),
        axis.text = element_text(size = 36),
        axis.title = element_text(size = 40),
        strip.text = element_text(size = 40, face = "bold", colour = "white"),
        plot.title = element_text(size = 48, face = "bold"),
        plot.subtitle = element_text(size = 44)
      )

      # Create the directory if it doesn't exist
      pathdir <- paste0("data/08_reporting/november/core_coefplots/", single_dep_var, "/") # nolint
      outfile <- paste0(pathdir, scopesubgroup, ".png") # nolint
      message("Saving plot to: ", outfile)
      if (!dir.exists(pathdir)) {
        message("Creating directory: ", pathdir)
        dir.create(pathdir, recursive = TRUE)
      }

      n_y_facet_rows <- length(unique(subgroup_coef_plot_data$scope_subgroup)) # nolint
      plot_height <- n_y_facet_rows * 6

      ggsave(
        outfile,
        coeffplot,
        width = 15,
        height = plot_height,
        dpi = 300
      )
    }
  }
}

# --- Function to generate strong coefficient plots (6 coefficients) ---
generate_strong_coef_plots <- function(coef_table) {
  # Define variables locally to avoid "no visible global binding" warnings
  local_coef_labels <- coef_labels
  local_strong_coef_vars <- strong_coef_vars
  local_observation_labels <- observation_labels
  local_scope_order <- scope_order
  local_field_order <- field_order
  local_subgroup_order <- subgroup_order
  local_strip_colors <- strip_colors
  local_field_colors <- field_colors
  local_dep_var_labels <- dep_var_labels

  # rename vars
  coef_table <- coef_table %>%
    mutate(
      dep_var = recode(dep_var, !!!local_dep_var_labels) # nolint
    )

  unique_dep_vars <- unique(coef_table$dep_var)

  for (single_dep_var in unique_dep_vars) {
    # Filter for strong coefficients only AND scope_Intent
    coef_plot_data <- coef_table %>%
      filter(
        treat_var %in% local_strong_coef_vars, # nolint
        dep_var == single_dep_var, # nolint
        scope == "scope_Intent" # Only include scope_Intent data for strong coefficients # nolint
      ) %>%
      mutate(
        scope = factor(gsub("scope_", "", scope),
          levels = gsub("scope_", "", local_scope_order)
        ), # nolint
        field = factor(gsub("field_", "", field), # nolint
          levels = gsub("field_", "", local_field_order)
        ), # nolint
        subgroup = factor(gsub("subgroup_", "", subgroup), # nolint
          levels = gsub("subgroup_", "", local_subgroup_order)
        ), # nolint
        scope_subgroup = factor(paste(scope, subgroup, sep = " - "),
          levels = unique(paste(scope, subgroup, sep = " - "))
        ), # nolint
        source_origin = factor(source_origin, levels = names(local_observation_labels)), # nolint
        source_origin = recode(source_origin, !!!local_observation_labels), # nolint
        treat_var = factor(treat_var, levels = local_strong_coef_vars), # nolint
        treat_var = recode(treat_var, !!!local_coef_labels) # nolint
      )

    coef_plot_data <- coef_plot_data %>%
      filter(!is.na(field) & !is.na(scope) & !is.na(subgroup)) # nolint

    for (scopesubgroup in unique(coef_plot_data$scope_subgroup)) {
      subgroup_coef_plot_data <- coef_plot_data %>%
        filter(scope_subgroup == scopesubgroup) # nolint

      # drop levels in column scope_subgroup if they have zero values
      vertical_levels <- levels(
        subgroup_coef_plot_data$scope_subgroup # nolint
      )[levels(subgroup_coef_plot_data$scope_subgroup) %in% # nolint
        unique(subgroup_coef_plot_data$scope_subgroup[
          subgroup_coef_plot_data$scope_subgroup != ""
        ])] # nolint
      vertical_strip_colors <- local_strip_colors[vertical_levels] # nolint

      horizontal_levels <- levels(
        subgroup_coef_plot_data$source_origin # nolint
      )[levels(subgroup_coef_plot_data$source_origin) %in% # nolint
        unique(subgroup_coef_plot_data$source_origin[
          subgroup_coef_plot_data$source_origin != ""
        ])] # nolint
      horizontal_strip_colors <- local_strip_colors[horizontal_levels] # nolint

      # Check if subgroup_coef_plot_data is empty
      if (nrow(subgroup_coef_plot_data) == 0) {
        message("No data for dep_var: ", single_dep_var)
        next
      }

      subgroup_coef_plot_data <- subgroup_coef_plot_data %>%
        filter(n_obs >= 200) # nolint

      # Create the plot
      coeffplot <- ggplot(
        subgroup_coef_plot_data,
        aes(x = estimate, y = treat_var, color = field) # nolint
      ) +
        geom_point(
          size = 4,
          position = position_dodge(width = 0.8)
        ) +
        geom_errorbarh(
          aes(xmin = estimate - 1.645 * std_error, xmax = estimate + 1.645 * std_error), # nolint
          height = 0, linewidth = 1, # thicker for 10% significance
          position = position_dodge(width = 0.8)
        ) +
        geom_errorbarh(
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.2, linewidth = 0.5, position = position_dodge(width = 0.8) # nolint
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
        ggh4x::facet_grid2(
          . ~ source_origin,
          scales = "free",
          independent = "x",
          space = "fixed",
          strip = ggh4x::strip_themed(
            background_y = ggh4x::elem_list_rect(
              fill = vertical_strip_colors
            ),
            background_x = ggh4x::elem_list_rect(
              fill = horizontal_strip_colors
            )
          )
        ) +
        scale_color_manual(values = local_field_colors) +
        labs(
          title = paste("Dependent Variable:", single_dep_var, " | ", scopesubgroup),  # nolint
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() +
        theme(
          axis.text.y = element_text(size = 40),
          axis.title.x = element_text(size = 36),
          axis.title.y = element_text(size = 36),
          strip.text = element_text(size = 30),
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
          panel.spacing = unit(2, "lines"),
          plot.margin = margin(1, 1, 1, 1, "cm"),
          legend.position = "bottom",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30)
        )

      all_fields_data <- subgroup_coef_plot_data %>%
        filter(field == "All Fields") # nolint

      # add counts of obs (subgroup all and ext.)
      coeffplot <- coeffplot + geom_text(
        data = all_fields_data,
        aes(label = paste0("n = ", n_obs)), # nolint
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 40,
        size = 12, color = "black"
      ) + theme(
        text = element_text(family = "muli", size = 40),
        axis.text = element_text(size = 36),
        axis.title = element_text(size = 40),
        strip.text = element_text(size = 40, face = "bold", colour = "white"),
        plot.title = element_text(size = 48, face = "bold"),
        plot.subtitle = element_text(size = 44)
      )

      # Create the directory if it doesn't exist
      pathdir <- paste0("data/08_reporting/november/strong_coefplots/", single_dep_var, "/") # nolint
      outfile <- paste0(pathdir, scopesubgroup, ".png")  # nolint
      message("Saving plot to: ", outfile)
      if (!dir.exists(pathdir)) {
        message("Creating directory: ", pathdir)
        dir.create(pathdir, recursive = TRUE)
      }

      n_y_facet_rows <- length(unique(subgroup_coef_plot_data$scope_subgroup))  # nolint
      plot_height <- n_y_facet_rows * 10 # Taller for more coefficients

      ggsave(
        outfile,
        coeffplot,
        width = 15,
        height = plot_height,
        dpi = 300
      )
    }
  }
}

# Run both functions
generate_core_coef_plots(all_tables_combined)
generate_strong_coef_plots(all_tables_combined)
