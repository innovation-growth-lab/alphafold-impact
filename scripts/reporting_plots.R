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
# ------------------------- #


dep_var_labels <- c(
  # Bibliometric & Impact Metrics
  "cited_by_count" = "Number of Citations Received",
  "ln1p_fwci" = "Field-Weighted Citation Impact (log-scale)",
  "patent_count" = "Number of Patent-Paper Citations",
  "patent_citation" = "Number of Patent-Patent Citations",
  "ca_count" = "Number of Clinical Article Citations",
  "num_publications" = "Number of Publications",

  # Disease & Clinical Relevance Metrics
  "ln1p_mesh_C" = "Disease Relevance (MeSH C-index, log-scale)",
  "num_diseases" = "Number of Structure-Associated Diseases",
  "num_uniprot_structures_w_disease" = "Number of Disease-Relevant Structures",
  "num_primary_submissions_w_disease" = "Number of Novel Disease-Relevant Structures", # nolint

  # Experimental Output & Quality Metrics
  "ln1p_resolution" = "Crystallographic Resolution (log-scale)",
  "ln1p_R_free" = "R-free Value (log-scale)",
  "num_pdb_ids" = "Number of PDB Structure IDs",
  "num_pdb_submissions" = "Number of PDB Submissions",
  "num_uniprot_structures" = "Number of UniProt Structures",

  # Novelty & Rarity Metrics
  "num_primary_submissions" = "Number of Novel UniProt Structures",
  "ln1p_organism_rarity_max" = "Max. Organism Rarity (log-scale)",
  "ln1p_max_tmscore" = "Str. Similarity (TM-score, log-scale)",
  "ln1p_max_fident" = "Max. Sequence Identity (log-scale)",
  "ln1p_max_score" = "Global Shape Similarity (RCSB, log-scale)",
  "num_uniprot_structures_w_rare_organisms" = "UniProt Structures / Rare Organisms", # nolint
  "num_primary_submissions_w_rare_organisms" = "Primary PDB Submissions / Rare Organisms", # nolint
  "num_uniprot_structures_w_low_similarity" = "UniProt Structures / Low Similarity", # nolint
  "num_primary_submissions_w_low_similarity" = "Primary PDB Submissions / Low Similarity", # nolint
  "ln1p_maxtmscore_lt_0.405" = "Max TM-score < 0.5"
)

# Updated coefficient labels with new naming
coef_labels <- c(
  # Core coefficients
  "af" = "AlphaFold",
  "ct_ai" = "AI Frontier",
  "ct_pp" = "PP Frontier",
  "ct_sb" = "SB Frontier",

  # Strong coefficients
  "af_intent_strong" = "AlphaFold - Method",
  "af_intent_weak" = "AlphaFold - Bkg.",
  "ct_ai_intent_strong" = "AI Frontier - Method",
  "ct_ai_intent_weak" = "AI Frontier - Bkg.",
  "ct_pp_intent_strong" = "PP Frontier - Method",
  "ct_pp_intent_weak" = "PP Frontier - Bkg.",
  "ct_sb_intent_strong" = "SB Frontier - Method",
  "ct_sb_intent_weak" = "SB Frontier - Bkg."
)

strip_colors <- c(
  "Laboratories" = "#00B2A2",
  "Researchers" = "#FF5836",
  # "Early Career Researchers" = "#FAB61B", # nolint
  "Citing Papers" = "#1F5DAD"
)

grouped_strip_colors <- c(
  "Laboratories" = "#00B2A2",
  "Researchers" = "#FF5836",
  "Early Career Researchers" = "#FAB61B",
  "Citing Papers" = "#1F5DAD"
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
  # "field_Medicine",
  # "field_Biochemistry",
  "field_All Fields"
)

# Updated observation order - removing ecr
observation_order <- c("papers", "nonecr", "ecr", "labs")

grouped_observation_order <- c("papers", "nonecr", "labs", "ecr")

observation_labels <- c(
  "papers" = "Citing Papers",
  "labs" = "Laboratories",
  "nonecr" = "Researchers",
  "ecr" = "Early Career Researchers"
)

grouped_observation_labels <- c(
  "papers" = "Citing Papers",
  "labs" = "Laboratories",
  "nonecr" = "Researchers",
  "ecr" = "Early Career Researchers"
)

# Define core coefficients and strong coefficients
core_coef_vars <- c("af", "ct_ai", "ct_pp", "ct_sb")
strong_coef_vars <- c(
  "af_intent_strong", "af_intent_weak",
  "ct_ai_intent_strong", "ct_ai_intent_weak",
  "ct_pp_intent_strong", "ct_pp_intent_weak",
  "ct_sb_intent_strong", "ct_sb_intent_weak"
)


# --- Function to generate coefficient plots ---
generate_core_coef_plots <- function(coef_table) { # nolint

  # rename vars
  coef_table <- coef_table %>%
    filter(dep_var %in% names(dep_var_labels)) %>% # nolint
    mutate(
      dep_var = recode(dep_var, !!!dep_var_labels) # nolint
    )

  unique_dep_vars <- unique(coef_table$dep_var)

  # rename "Molecular Biology" to "Biochemistry"
  coef_table <- coef_table %>% # nolint
    mutate( # nolint
      field = recode(field, "field_Molecular Biology" = "field_Biochemistry")
    )

  for (single_dep_var in unique_dep_vars) {
    coef_plot_data <- coef_table %>% # nolint
      filter(
        treat_var %in% core_coef_vars, # nolint
        source_origin %in% names(observation_labels), # nolint
        dep_var == single_dep_var # nolint
      ) %>%
      mutate( # nolint
        scope = factor(gsub("scope_", "", scope), # nolint
          levels = gsub("scope_", "", scope_order) # nolint
        ),
        field = factor(
          gsub("field_", "", field),
          levels = gsub("field_", "", field_order)
        ),
        subgroup = factor(
          gsub("subgroup_", "", subgroup), # nolint
          levels = gsub("subgroup_", "", subgroup_order)
        ),
        scope_subgroup = factor(paste(scope, subgroup, sep = " - "),
          levels = unique(paste(scope, subgroup, sep = " - "))
        ), # nolint
        source_origin = factor(source_origin, levels = observation_order), # nolint
        source_origin = recode(source_origin, !!!observation_labels), # nolint
        treat_var = factor(
          treat_var,
          levels = rev(core_coef_vars) # nolint
        ),
        treat_var = recode(treat_var, !!!coef_labels) # nolint
      )

    coef_plot_data <- coef_plot_data %>%
      filter(!is.na(field) & !is.na(scope) & !is.na(subgroup)) # nolint

    for (scopesubgroup in unique(coef_plot_data$scope_subgroup)) { # nolint
      subgroup_coef_plot_data <- coef_plot_data %>%
        filter(scope_subgroup == scopesubgroup) %>% # nolint
        filter(n_obs >= 1000, estimate >= -10, estimate <= 10) # nolint

      for (field in unique(subgroup_coef_plot_data$field)) { # nolint
        field_subgroup_coef_plot_data <- subgroup_coef_plot_data %>%
          filter(field == !!field)

        horizontal_levels <- levels(
          field_subgroup_coef_plot_data$source_origin
        )[levels(field_subgroup_coef_plot_data$source_origin) %in%
          unique(field_subgroup_coef_plot_data$source_origin[field_subgroup_coef_plot_data$source_origin != ""])] # nolint
        horizontal_strip_colors <- strip_colors[horizontal_levels]

        # Check if subgroup_coef_plot_data is empty
        if (nrow(field_subgroup_coef_plot_data) == 0) {
          message("No data for dep_var: ", single_dep_var) # nolint
          next
        }

        # Create the plot
        coeffplot <- ggplot(
          field_subgroup_coef_plot_data,
          aes(x = estimate, y = treat_var) # nolint
        ) +
          geom_point( # nolint
            size = 4
          ) +
          geom_errorbarh( # nolint
            aes(xmin = estimate - 1.645 * std_error, xmax = estimate + 1.645 * std_error), # nolint
            height = 0, linewidth = 1
          ) +
          geom_errorbarh( # nolint
            aes(xmin = conf_low, xmax = conf_high), # nolint
            height = 0.2, linewidth = 0.5
          ) +
          # scale_color_manual(values = depth_colors) +
          geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
          ggh4x::facet_grid2(
            . ~ source_origin,
            scales = "free",
            independent = "x",
            space = "fixed",
            labeller = label_wrap_gen(20),
            strip = ggh4x::strip_themed(
              background_x = ggh4x::elem_list_rect(
                fill = horizontal_strip_colors
              )
            )
          ) + # nolint
          scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
          labs( # nolint
            title = paste("Dependent Variable:", single_dep_var), # nolint
            x = "Estimate (with 95% CI)",
            y = "Coefficient Variable"
          ) +
          theme_classic() + # nolint
          theme( # nolint
            axis.text.y = element_text(size = 56), # nolint
            axis.title.x = element_text(size = 56), # nolint
            axis.title.y = element_text(size = 56), # nolint
            strip.text = element_text(size = 45), # nolint
            panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"), # nolint
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
            panel.spacing = unit(1.5, "lines"), # nolint
            legend.position = "none",
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm") # nolint
          )

        # add counts of obs (subgroup all and ext.)
        coeffplot <- coeffplot + geom_text( # nolint
          data = field_subgroup_coef_plot_data, # nolint
          aes(label = paste0("n = ", n_obs)), # nolint
          x = Inf, y = Inf,
          hjust = 1.1, vjust = 21.75,
          size = 16, color = "black"
        ) + theme(
          text = element_text(family = "muli", size = 56), # Base font size
          axis.text = element_text(size = 56), # Axis tick labels
          axis.title = element_text(size = 60), # Axis titles
          strip.text = element_text(size = 60, face = "bold", colour = "white"), # nolint
          plot.title = element_text(size = 66, face = "bold"), # Title
          plot.subtitle = element_text(size = 44)
        )

        # Create the directory if it doesn't exist
        pathdir <- paste0(
          "data/08_reporting/2025Q1/coefplots/core/",
          single_dep_var, "/", field, "/"
        )
        outfile <- paste0(pathdir, scopesubgroup, ".png")
        message("Saving plot to: ", outfile)
        if (!dir.exists(pathdir)) {
          message("Creating directory: ", pathdir)
          dir.create(pathdir, recursive = TRUE)
        }

        n_y_facet_rows <- length(
          unique(field_subgroup_coef_plot_data$scope_subgroup)
        )
        plot_height <- n_y_facet_rows * 4.5 # nolint

        # Calculate width based on number of source_origin facets
        n_source_facets <- length(unique(field_subgroup_coef_plot_data$source_origin))
        plot_width <- n_source_facets * 4.5 # 4.5 units per facet

        ggsave( # nolint
          outfile, # nolint
          coeffplot,
          width = plot_width,
          height = plot_height,
          dpi = 300
        )
      }
    }
  }
}


# --- Function to generate strong coefficient plots ---
generate_strong_coef_plots <- function(coef_table) { # nolint

  # rename vars
  coef_table <- coef_table %>%
    filter(dep_var %in% names(dep_var_labels)) %>% # nolint
    mutate(
      dep_var = recode(dep_var, !!!dep_var_labels) # nolint
    )

  unique_dep_vars <- unique(coef_table$dep_var)

  # rename "Molecular Biology" to "Biochemistry"
  coef_table <- coef_table %>% # nolint
    mutate( # nolint
      field = recode(field, "field_Molecular Biology" = "field_Biochemistry")
    )

  for (single_dep_var in unique_dep_vars) {
    # For specific variables, filter source origins
    if (single_dep_var == "Number of Publications") {
      allowed_sources <- c("nonecr", "labs")
    } else if (single_dep_var %in% c(
      "Number of Citations Received",
      "Field-Weighted Citation Impact (log-scale)",
      "Number of Clinical Article Citations",
      "Number of Patent-Paper Citations"
    )) {
      allowed_sources <- c("papers", "nonecr", "labs")
    } else {
      allowed_sources <- names(observation_labels)
    }

    coef_plot_data <- coef_table %>% # nolint
      filter(
        treat_var %in% strong_coef_vars, # nolint
        source_origin %in% allowed_sources, # nolint
        dep_var == single_dep_var # nolint
      ) %>%
      mutate( # nolint
        scope = factor(gsub("scope_", "", scope), # nolint
          levels = gsub("scope_", "", scope_order) # nolint
        ),
        field = factor(
          gsub("field_", "", field),
          levels = gsub("field_", "", field_order)
        ),
        subgroup = factor(
          gsub("subgroup_", "", subgroup), # nolint
          levels = gsub("subgroup_", "", subgroup_order)
        ),
        scope_subgroup = factor(paste(scope, subgroup, sep = " - "),
          levels = unique(paste(scope, subgroup, sep = " - "))
        ), # nolint
        source_origin = factor(source_origin, levels = observation_order), # nolint
        source_origin = recode(source_origin, !!!observation_labels), # nolint
        treat_var = factor(
          treat_var,
          levels = rev(strong_coef_vars) # nolint
        ),
        treat_var = recode(treat_var, !!!coef_labels) # nolint
      )

    coef_plot_data <- coef_plot_data %>%
      filter(!is.na(field) & !is.na(scope) & !is.na(subgroup)) # nolint

    for (scopesubgroup in unique(coef_plot_data$scope_subgroup)) { # nolint
      subgroup_coef_plot_data <- coef_plot_data %>%
        filter(scope_subgroup == scopesubgroup) %>% # nolint
        filter(n_obs >= 1000, estimate >= -10, estimate <= 10) # nolint

      for (field in unique(subgroup_coef_plot_data$field)) { # nolint
        field_subgroup_coef_plot_data <- subgroup_coef_plot_data %>%
          filter(field == !!field)

        horizontal_levels <- levels(
          field_subgroup_coef_plot_data$source_origin
        )[levels(field_subgroup_coef_plot_data$source_origin) %in%
          unique(field_subgroup_coef_plot_data$source_origin[field_subgroup_coef_plot_data$source_origin != ""])] # nolint
        horizontal_strip_colors <- strip_colors[horizontal_levels]

        # Check if subgroup_coef_plot_data is empty
        if (nrow(field_subgroup_coef_plot_data) == 0) {
          message("No data for dep_var: ", single_dep_var) # nolint
          next
        }

        # Create the plot
        coeffplot <- ggplot(
          field_subgroup_coef_plot_data,
          aes(x = estimate, y = treat_var) # nolint
        ) +
          geom_point( # nolint
            size = 4
          ) +
          geom_errorbarh( # nolint
            aes(xmin = estimate - 1.645 * std_error, xmax = estimate + 1.645 * std_error), # nolint
            height = 0, linewidth = 1
          ) +
          geom_errorbarh( # nolint
            aes(xmin = conf_low, xmax = conf_high), # nolint
            height = 0.2, linewidth = 0.5
          ) +
          # scale_color_manual(values = depth_colors) +
          geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
          ggh4x::facet_grid2(
            . ~ source_origin,
            scales = "free",
            independent = "x",
            space = "fixed",
            labeller = label_wrap_gen(20),
            strip = ggh4x::strip_themed(
              background_x = ggh4x::elem_list_rect(
                fill = horizontal_strip_colors
              )
            )
          ) + # nolint
          scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
          labs( # nolint
            title = paste("Dependent Variable:", single_dep_var), # nolint
            x = "Estimate (with 95% CI)",
            y = "Coefficient Variable"
          ) +
          theme_classic() + # nolint
          theme( # nolint
            axis.text.y = element_text(size = 56), # nolint
            axis.title.x = element_text(size = 56), # nolint
            axis.title.y = element_text(size = 56), # nolint
            strip.text = element_text(size = 45), # nolint
            panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"), # nolint
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
            panel.spacing = unit(1.5, "lines"), # nolint
            legend.position = "none",
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm") # nolint
          )

        # add counts of obs (subgroup all and ext.)
        coeffplot <- coeffplot + geom_text( # nolint
          data = field_subgroup_coef_plot_data, # nolint
          aes(label = paste0("n = ", n_obs)), # nolint
          x = Inf, y = Inf,
          hjust = 1.1, vjust = 21.75,
          size = 16, color = "black"
        ) + theme(
          text = element_text(family = "muli", size = 56), # Base font size
          axis.text = element_text(size = 56), # Axis tick labels
          axis.title = element_text(size = 60), # Axis titles
          strip.text = element_text(size = 60, face = "bold", colour = "white"), # nolint
          plot.title = element_text(size = 66, face = "bold"), # Title
          plot.subtitle = element_text(size = 44)
        )

        # Create the directory if it doesn't exist
        pathdir <- paste0(
          "data/08_reporting/2025Q1/coefplots/strong/",
          single_dep_var, "/", field, "/"
        )
        outfile <- paste0(pathdir, scopesubgroup, ".png")
        message("Saving plot to: ", outfile)
        if (!dir.exists(pathdir)) {
          message("Creating directory: ", pathdir)
          dir.create(pathdir, recursive = TRUE)
        }

        n_y_facet_rows <- length(
          unique(field_subgroup_coef_plot_data$scope_subgroup)
        )
        plot_height <- n_y_facet_rows * 4.5 # nolint

        # Calculate width based on number of source_origin facets
        n_source_facets <- length(unique(field_subgroup_coef_plot_data$source_origin))
        plot_width <- n_source_facets * 4.5 # 4.5 units per facet

        ggsave( # nolint
          outfile, # nolint
          coeffplot,
          width = plot_width,
          height = plot_height,
          dpi = 300
        )
      }
    }
  }
}

# --- Function to generate grouped coefficient plots ---
generate_grouped_coef_plots <- function(coef_table) {
  # Prepare the data
  plot_data <- coef_table %>%
    filter(
      treat_var %in% core_coef_vars, # nolint
      source_origin %in% names(grouped_observation_labels), # nolint
      dep_var %in% names(dep_var_labels), # nolint
      scope == "scope_All", # nolint
      subgroup == "subgroup_All PDB" # nolint
    ) %>%
    mutate(
      dep_var = recode(dep_var, !!!dep_var_labels),
      dep_var = factor(dep_var, levels = unname(dep_var_labels)),
      source_origin = factor(source_origin, levels = grouped_observation_order),
      source_origin = recode(source_origin, !!!grouped_observation_labels),
      treat_var = factor(treat_var, levels = rev(core_coef_vars)),
      treat_var = recode(treat_var, !!!coef_labels)
    ) %>%
    filter(
      n_obs >= 1000, # nolint
      estimate >= -10, # nolint
      estimate <= 10 # nolint
    )

  # Generate plot for each source origin
  for (source in unique(plot_data$source_origin)) {
    for (ufield in unique(plot_data$field)) {
      if (ufield == "field_All Fields") {
        source_data <- plot_data %>%
          filter(source_origin == source, field == ufield) # nolint

        # Create faceted plot
        grouped_plot <- ggplot(
          source_data,
          aes(x = estimate, y = treat_var) # nolint
        ) +
          geom_point(size = 12) +
          geom_errorbarh(
            aes(
              xmin = estimate - 1.645 * std_error, # nolint
              xmax = estimate + 1.645 * std_error # nolint
            ),
            height = 0,
            linewidth = 2
          ) +
          geom_errorbarh(
            aes(xmin = conf_low, xmax = conf_high), # nolint
            height = 0.4,
            linewidth = 1.2
          ) +
          geom_vline(xintercept = 0, color = "black", linewidth = 1) +
          facet_wrap(
            ~dep_var,
            scales = "free_x",
            ncol = 3
          ) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
          labs(
            title = paste("Coefficient Estimates for", source, "|", ufield),
            x = "Estimate (with 95% CI)",
            y = "Coefficient Variable"
          ) +
          theme_classic() +
          theme(
            text = element_text(family = "muli"),
            axis.text.y = element_text(size = 84),
            axis.text.x = element_text(size = 90),
            axis.title = element_text(size = 90),
            strip.text = element_text(
              size = 72, # Increased from 94
              face = "bold",
              color = "white",
              margin = margin(t = 10, b = 10) # Add some vertical padding
            ),
            strip.background = element_rect(
              fill = grouped_strip_colors[source],
              color = "black", # Add border
              linewidth = 1
            ),
            panel.grid.major.x = element_line(linewidth = 1.2, color = "grey"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # nolint
            panel.spacing = unit(2, "lines"),
            plot.title = element_text(size = 102, face = "bold"),
            plot.margin = margin(1, 1, 1, 1, "cm")
          )

        # Add observation counts
        grouped_plot <- grouped_plot + geom_text(
          aes(label = paste0("n = ", n_obs), x = Inf, y = -Inf), # nolint
          hjust = 1.1,
          vjust = -0.5,
          size = 24
        )

        # Create directory and save plot
        pathdir <- paste0(
          "data/08_reporting/2025Q1/coefplots/grouped/",
          ufield, "/"
        )
        if (!dir.exists(pathdir)) {
          dir.create(pathdir, recursive = TRUE)
        }

        outfile <- paste0(pathdir, source, "_grouped.png")
        message("Saving plot to: ", outfile)

        # Calculate dimensions based on number of dependent variables
        n_deps <- length(unique(source_data$dep_var))
        n_cols <- 3
        n_rows <- ceiling(n_deps / n_cols)
        plot_height <- max(12, n_rows * 4) # Increased base height
        plot_width <- min(24, n_cols * 8) # Increased base width

        ggsave(
          outfile,
          grouped_plot,
          width = plot_width,
          height = plot_height,
          dpi = 300
        )
      }
    }
  }
}


# --- Function to generate grouped coefficient plots ---
generate_grouped_highpdb_coef_plots <- function(coef_table) {
  # Prepare the data
  plot_data <- coef_table %>%
    filter(
      treat_var %in% core_coef_vars, # nolint
      source_origin %in% names(grouped_observation_labels), # nolint
      dep_var %in% names(dep_var_labels), # nolint
      scope == "scope_All", # nolint
      subgroup == "subgroup_High PDB" # nolint
    ) %>%
    mutate(
      dep_var = recode(dep_var, !!!dep_var_labels),
      dep_var = factor(dep_var, levels = unname(dep_var_labels)),
      source_origin = factor(source_origin, levels = grouped_observation_order),
      source_origin = recode(source_origin, !!!grouped_observation_labels),
      treat_var = factor(treat_var, levels = rev(core_coef_vars)),
      treat_var = recode(treat_var, !!!coef_labels)
    ) %>%
    filter(
      n_obs >= 1000, # nolint
      estimate >= -10, # nolint
      estimate <= 10
    )

  # Generate plot for each source origin
  for (source in unique(plot_data$source_origin)) {
    for (ufield in unique(plot_data$field)) {
      source_data <- plot_data %>%
        filter(source_origin == source, field == ufield) # nolint

      # Create faceted plot
      grouped_plot <- ggplot(
        source_data,
        aes(x = estimate, y = treat_var) # nolint
      ) +
        geom_point(size = 12) +
        geom_errorbarh(
          aes(
            xmin = estimate - 1.645 * std_error, # nolint
            xmax = estimate + 1.645 * std_error
          ),
          height = 0,
          linewidth = 2
        ) +
        geom_errorbarh(
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.4,
          linewidth = 1.2
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
        facet_wrap(
          ~dep_var,
          scales = "free_x",
          ncol = 3
        ) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        labs(
          title = paste("Coefficient Estimates for", source, "| High PDB |", ufield), # nolint
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() +
        theme(
          text = element_text(family = "muli"),
          axis.text.y = element_text(size = 84),
          axis.text.x = element_text(size = 90),
          axis.title = element_text(size = 90),
          strip.text = element_text(
            size = 72, # Increased from 94
            face = "bold",
            color = "white",
            margin = margin(t = 10, b = 10) # Add some vertical padding
          ),
          strip.background = element_rect(
            fill = grouped_strip_colors[source],
            color = "black", # Add border
            linewidth = 1
          ),
          panel.grid.major.x = element_line(linewidth = 1.2, color = "grey"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # nolint
          panel.spacing = unit(2, "lines"),
          plot.title = element_text(size = 102, face = "bold"),
          plot.margin = margin(1, 1, 1, 1, "cm")
        )

      # Add observation counts
      grouped_plot <- grouped_plot + geom_text(
        aes(label = paste0("n = ", n_obs), x = Inf, y = -Inf), # nolint
        hjust = 1.1,
        vjust = -0.5,
        size = 24
      )

      # Create directory and save plot
      pathdir <- paste0(
        "data/08_reporting/2025Q1/coefplots/grouped_highpdb/",
        ufield, "/"
      )
      if (!dir.exists(pathdir)) {
        dir.create(pathdir, recursive = TRUE)
      }

      outfile <- paste0(pathdir, source, "_grouped.png")
      message("Saving plot to: ", outfile)

      # Calculate dimensions based on number of dependent variables
      n_deps <- length(unique(source_data$dep_var))
      n_cols <- 3
      n_rows <- ceiling(n_deps / n_cols)
      plot_height <- max(12, n_rows * 4) # Increased base height
      plot_width <- min(24, n_cols * 8) # Increased base width

      ggsave(
        outfile,
        grouped_plot,
        width = plot_width,
        height = plot_height,
        dpi = 300
      )
    }
  }
}



# --- Function to generate grouped coefficient plots ---
generate_grouped_strong_coef_plots <- function(coef_table) {
  # Prepare the data
  plot_data <- coef_table %>%
    filter(
      treat_var %in% strong_coef_vars, # nolint
      source_origin %in% names(grouped_observation_labels), # nolint
      dep_var %in% names(dep_var_labels), # nolint
      scope == "scope_Intent", # nolint
      subgroup == "subgroup_All PDB" # nolint
    ) %>%
    mutate(
      dep_var = recode(dep_var, !!!dep_var_labels),
      dep_var = factor(dep_var, levels = unname(dep_var_labels)),
      source_origin = factor(source_origin, levels = grouped_observation_order),
      source_origin = recode(source_origin, !!!grouped_observation_labels),
      treat_var = factor(treat_var, levels = rev(strong_coef_vars)),
      treat_var = recode(treat_var, !!!coef_labels)
    ) %>%
    filter(
      n_obs >= 1000, # nolint
      estimate >= -10, # nolint
      estimate <= 10
    )

  # Generate plot for each source origin
  for (source in unique(plot_data$source_origin)) {
    for (ufield in unique(plot_data$field)) {
      source_data <- plot_data %>%
        filter(source_origin == source, field == ufield) # nolint

      # Create faceted plot
      grouped_plot <- ggplot(
        source_data,
        aes(x = estimate, y = treat_var) # nolint
      ) +
        geom_point(size = 12) +
        geom_errorbarh(
          aes(
            xmin = estimate - 1.645 * std_error, # nolint
            xmax = estimate + 1.645 * std_error
          ),
          height = 0,
          linewidth = 2
        ) +
        geom_errorbarh(
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.4,
          linewidth = 1.2
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) +
        facet_wrap(
          ~dep_var,
          scales = "free_x",
          ncol = 3
        ) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        labs(
          title = paste("Coefficient Estimates for", source, "|", ufield),
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() +
        theme(
          text = element_text(family = "muli"),
          axis.text.y = element_text(size = 84),
          axis.text.x = element_text(size = 90),
          axis.title = element_text(size = 90),
          strip.text = element_text(
            size = 72, # Increased from 94
            face = "bold",
            color = "white",
            margin = margin(t = 10, b = 10) # Add some vertical padding
          ),
          strip.background = element_rect(
            fill = grouped_strip_colors[source],
            color = "black", # Add border
            linewidth = 1
          ),
          panel.grid.major.x = element_line(linewidth = 1.2, color = "grey"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # nolint
          panel.spacing = unit(2, "lines"),
          plot.title = element_text(size = 102, face = "bold"),
          plot.margin = margin(1, 1, 1, 1, "cm")
        )

      # Add observation counts
      grouped_plot <- grouped_plot + geom_text(
        aes(label = paste0("n = ", n_obs), x = Inf, y = -Inf), # nolint
        hjust = 1.1,
        vjust = -0.5,
        size = 24
      )

      # Create directory and save plot
      pathdir <- paste0(
        "data/08_reporting/2025Q1/coefplots/grouped_strong/",
        ufield, "/"
      )
      if (!dir.exists(pathdir)) {
        dir.create(pathdir, recursive = TRUE)
      }

      outfile <- paste0(pathdir, source, "_grouped.png")
      message("Saving plot to: ", outfile)

      # Calculate dimensions based on number of dependent variables
      n_deps <- length(unique(source_data$dep_var))
      n_cols <- 3
      n_rows <- ceiling(n_deps / n_cols)
      plot_height <- max(16, n_rows * 6) # Increased base height
      plot_width <- min(30, n_cols * 10) # Increased base width

      ggsave(
        outfile,
        grouped_plot,
        width = plot_width,
        height = plot_height,
        dpi = 300
      )
    }
  }
}

# Run the function
generate_core_coef_plots(all_tables_combined)
generate_strong_coef_plots(all_tables_combined)
generate_grouped_coef_plots(all_tables_combined)
generate_grouped_highpdb_coef_plots(all_tables_combined)
generate_grouped_strong_coef_plots(all_tables_combined)
