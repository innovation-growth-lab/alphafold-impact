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
    rds_files <- list.files(path, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

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

# if source_origin is papers, change subgroup_order "All PDB" to "All PDB - CEM"
all_tables_combined <- all_tables_combined %>%
  mutate( # nolint
    subgroup = ifelse(
      source_origin == "papers" & subgroup == "subgroup_All PDB",
      "subgroup_All PDB - CEM",
      subgroup
    )
  )

# if source origin is papers, add _ind to treat_var
all_tables_combined <- all_tables_combined %>%
  mutate( # nolint
    treat_var = ifelse(
      source_origin == "papers" & grepl("_ind", treat_var),
      paste0(treat_var, "_ind"),
      treat_var
    )
  )

# if source origin is papers, also relabel af:strong1 to strong_af, ct_ai:strong1 to strong_ct_ai, ct_noai:strong1 to strong_ct_noai
all_tables_combined <- all_tables_combined %>%
  mutate( # nolint
    treat_var = recode(treat_var, "af:strong1" = "strong_af", "ct_ai:strong1" = "strong_ct_ai", "ct_noai:strong1" = "strong_ct_noai")
  )


# count values of each source origin
all_tables_combined %>%
  count(source_origin) %>%
  print()

# change treat_var num_pdb_submissions to num_pdb_ids

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
  "ln1p_cit_0" = "ln (1 + Citations at month < 12)",
  "ln1p_cit_1" = "ln (1 + Citations at month < 24)",
  "ln1p_fwci" = "ln (1 + Field-Weighted Citation Impact)",
  "logit_cit_norm_perc" = "logit (Citation percentile)",
  "patent_count" = "Patent citation count",
  "patent_citation" = "Patent-to-patent citation count",
  "ca_count" = "Clinical article citation count",
  "ln1p_resolution" = "ln (1 + Resolution)",
  "ln1p_R_free" = "ln (1 + R-free)",
  "num_pdb_ids" = "PDB submissions",
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

coef_labels <- c(
  "af_ind" = "AlphaFold",
  "ct_ai_ind" = "AI Frontier",
  "ct_noai_ind" = "No AI Frontier",
  "strong_af_ind" = "AlphaFold - Method",
  "strong_ct_ai_ind" = "AI Frontier - Method",
  "strong_ct_noai_ind" = "No AI Frontier - Method",
  "af" = "AlphaFold",
  "ct_ai" = "AI Frontier",
  "ct_noai" = "No AI Frontier",
  "strong_af" = "AlphaFold - Method",
  "strong_ct_ai" = "AI Frontier - Method",
  "strong_ct_noai" = "No AI Frontier - Method"
)

strip_colors <- c(
  "All Groups - All PDB - CEM" = "#092640",
  "All Groups - High PDB - CEM" = "lightcoral",
  "Foundational - All PDB - CEM" = "#092640",
  "Foundational - High PDB - CEM" = "lightcoral",
  "Applied - All PDB - CEM" = "#092640",
  "Applied - High PDB - CEM" = "lightcoral",
  "Laboratories" = "#00B2A2",
  "Established Researchers" = "#FF5836",
  "Early Career Researchers" = "#FAB61B",
  "Citation Chains" = "#1F5DAD"
)

depth_colors <- c(
  "All Groups" = "#1f77b4", # Blue
  "Foundational" = "#ff7f0e", # Orange
  "Applied" = "#2ca02c" # Green
)


# --- Variable definitions ---
# Set desired orders for variables and names
subgroup_order <- c("All PDB - CEM", "High PDB - CEM")
depth_order <- c("All Groups")#, "Foundational", "Applied")
field_order <- c(
  "field_All Fields",
  "field_Molecular Biology",
  "field_Medicine"
)

observation_order <- c("papers", "nonecr", "ecr", "labs")

observation_labels <- c(
  "papers" = "Citation Chains",
  "labs" = "Laboratories",
  "nonecr" = "Established Researchers",
  "ecr" = "Early Career Researchers"
)


coef_order <- c(
  "strong_af:strong_ct_noai_ind", "strong_af:strong_ct_ai_ind", "strong_ct_noai_ind", "strong_ct_ai_ind", "strong_af_ind", # nolint
  "af:ct_noai_ind", "af:ct_ai_ind", "ct_noai_ind", "ct_ai_ind", "af_ind",
  "strong_af:strong_ct_noai", "strong_af:strong_ct_ai", "strong_ct_noai", "strong_ct_ai", "strong_af", # nolint
  "af:ct_noai", "af:ct_ai", "ct_noai", "ct_ai", "af"
)

# --- Function to generate coefficient plots ---
generate_coef_plots <- function(coef_table) { # nolint

  # rename vars
  coef_table <- coef_table %>% # nolint
    mutate( # nolint
      dep_var = recode(dep_var, !!!dep_var_labels), # nolint
      treat_var = ifelse(
        grepl("_ind", indep_vars), # nolint
        paste0(treat_var, "_ind"), # nolint
        treat_var
      )
    )

  unique_dep_vars <- unique(coef_table$dep_var)

  for (single_dep_var in unique_dep_vars) {
    coef_plot_data <- coef_table %>% # nolint
      filter(treat_var %in% names(coef_labels), dep_var == single_dep_var) %>% # nolint
      mutate( # nolint
        depth = factor(gsub("depth_", "", depth), levels = gsub("depth_", "", depth_order)), # nolint
        field = factor(gsub("field_", "", field), levels = gsub("field_", "", field_order)), # nolint
        subgroup = factor(gsub("subgroup_", "", subgroup), levels = gsub("subgroup_", "", subgroup_order)), # nolint
        depth_subgroup = factor(paste(field, subgroup, sep = " - "), levels = unique(paste(field, subgroup, sep = " - "))), # nolint
        source_origin = factor(source_origin, levels = observation_order), # nolint
        source_origin = recode(source_origin, !!!observation_labels), # nolint
        treat_var = factor(
          treat_var,
          levels = coef_order
        ),
        treat_var = recode(treat_var, !!!coef_labels) # nolint
      )

      coef_plot_data <- coef_plot_data %>% # nolint
        filter(!is.na(field) & !is.na(depth) & !is.na(subgroup)) # nolint

    for (depthsubgroup in unique(coef_plot_data$depth_subgroup)) {
      subgroup_coef_plot_data <- coef_plot_data %>% # nolint
        filter(depth_subgroup == depthsubgroup)

      # drop levels in column depth_subgroup if they have zero values
      vertical_levels <- levels(
        subgroup_coef_plot_data$depth_subgroup
      )[levels(subgroup_coef_plot_data$depth_subgroup) %in%
        unique(subgroup_coef_plot_data$depth_subgroup[subgroup_coef_plot_data$depth_subgroup != ""])] # nolint
      vertical_strip_colors <- strip_colors[vertical_levels]


      subgroup_coef_plot_data <- subgroup_coef_plot_data %>% # nolint
        filter(n_obs >= 650, estimate >= -10, estimate <= 10)


      horizontal_levels <- levels(
        subgroup_coef_plot_data$source_origin
      )[levels(subgroup_coef_plot_data$source_origin) %in%
        unique(subgroup_coef_plot_data$source_origin[subgroup_coef_plot_data$source_origin != ""])] # nolint
      horizontal_strip_colors <- strip_colors[horizontal_levels]

      # Check if subgroup_coef_plot_data is empty
      if (nrow(subgroup_coef_plot_data) == 0) {
        message("No data for dep_var: ", single_dep_var) # nolint
        next
      }

      # Create the plot
      coeffplot <- ggplot( # nolint
        subgroup_coef_plot_data,
        aes(x = estimate, y = treat_var) # , color = depth) # nolint
      ) +
        geom_point( # nolint
          size = 4
          # position = position_dodge(width = 0.5)
        ) +
        geom_errorbarh( # nolint
          aes(xmin = estimate - 1.645 * std_error, xmax = estimate + 1.645 * std_error), # nolint
          height = 0, linewidth = 1 # , # thicker for 10% significance
          # position = position_dodge(width = 0.5)
        ) +
        geom_errorbarh( # nolint
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.2, linewidth = 0.5 #  position = position_dodge(width = 0.5) # thinner for 5% significance
        ) +
        # scale_color_manual(values = depth_colors) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
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
        ) + # nolint
        labs( # nolint
          title = paste("Dependent Variable:", single_dep_var, " | ", depthsubgroup), # nolint
          # subtitle = paste("Field:", field_label), # nolint
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() + # nolint
        theme( # nolint
          axis.text.y = element_text(size = 40), # nolint
          axis.title.x = element_text(size = 36), # nolint
          axis.title.y = element_text(size = 36), # nolint
          strip.text = element_text(size = 30), # nolint
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"), # nolint
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
          panel.spacing = unit(2, "lines"), # nolint
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1, "cm") # nolint
        )

      # add counts of obs (subgroup all and ext.)
      coeffplot <- coeffplot + geom_text( # nolint
        data = subgroup_coef_plot_data, # nolint
        aes(label = paste0("n = ", n_obs)), # nolint
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 25.5,
        size = 12, color = "black"
      ) + theme(
        text = element_text(family = "muli", size = 40), # Base font size
        axis.text = element_text(size = 36), # Axis tick labels
        axis.title = element_text(size = 40), # Axis titles
        strip.text = element_text(size = 40, face = "bold", colour = "white"), # Facet strip text
        plot.title = element_text(size = 48, face = "bold"), # Title
        plot.subtitle = element_text(size = 44)
      )

      # Create the directory if it doesn't exist
      pathdir <- paste0("data/08_reporting/november/coefplots/", single_dep_var, "/")
      outfile <- paste0(pathdir, depthsubgroup, ".png")
      message("Saving plot to: ", outfile)
      if (!dir.exists(pathdir)) {
        message("Creating directory: ", pathdir)
        dir.create(pathdir, recursive = TRUE)
      }

      n_y_facet_rows <- length(unique(subgroup_coef_plot_data$depth_subgroup))
      plot_height <- n_y_facet_rows * 4.5 # nolint

      ggsave( # nolint
        outfile, # nolint
        coeffplot,
        width = 15,
        height = plot_height,
        dpi = 300
      )
    }
  }
}

# Run the function
generate_coef_plots(all_tables_combined)
