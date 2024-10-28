figures <- "data/05_model_output/papers/figures/"

if (!dir.exists(figures)) {
  dir.create(figures, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# FIGURE GENERATION
# ------------------------------------------------------------------------------

extract_coefficients <- function(results, dep_vars, subsets, cov_sets, fe_list, treat_vars, treat_var_interest) { # nolint
  coef_data <- list()

  # Iterate over dependent variables
  for (dep_var in dep_vars) {
    # Iterate over subsets
    for (sub in subsets) {
      # Iterate over covariate sets, fixed effects, and treatment variables # nolint
      for (cov_set in cov_sets) {
        for (fe in fe_list) {
          for (treat_var in treat_vars) {
            # Build the result name
            result_name <- paste0(
              sub, "__", dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
            )

            # Check if result exists
            if (result_name %in% names(results)) {
              message("Extracting coefficients for: ", result_name)
              # Check if the result is not null
              if (!is.null(results[[result_name]])) {
                message("Found coefficients for: ", result_name)
                coef_info <- summary(results[[result_name]])$coeftable
                n_obs <- summary(results[[result_name]])$nobs

                # Iterate over the list of treatment variables of interest
                for (treat_var_interest_item in treat_var_interest) {
                  if (treat_var_interest_item %in% rownames(coef_info)) {
                    parts <- strsplit(result_name, "__")[[1]]
                    depth <- parts[1]
                    field <- parts[2]
                    tech <- parts[3]
                    pdb <- parts[4]
                    dep_var <- parts[5]
                    cov_set <- parts[6]
                    fe <- parts[7]
                    indep_vars <- parts[8]

                    coef_data[[length(coef_data) + 1]] <- data.frame(
                      depth = depth,
                      field = field,
                      tech = tech,
                      pdb = pdb,
                      treat_var = treat_var_interest_item,
                      dep_var = dep_var,
                      indep_vars = indep_vars,
                      estimate = coef_info[treat_var_interest_item, "Estimate"],
                      std_error = coef_info[treat_var_interest_item, "Std. Error"], # nolint
                      conf_low = coef_info[treat_var_interest_item, "Estimate"] - 1.96 * coef_info[treat_var_interest_item, "Std. Error"], # nolint
                      conf_high = coef_info[treat_var_interest_item, "Estimate"] + 1.96 * coef_info[treat_var_interest_item, "Std. Error"], # nolint
                      n_obs = n_obs
                    )
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return(do.call(rbind, coef_data)) # Combine list into a data frame
}

# --------------------------------------------------------------------------
# PLOT GENERATION WITH EXTENSIVE AND INTENSIVE MARGINS
# --------------------------------------------------------------------------

# --- Variable definitions ---
# Set desired orders for variables and names
pdb_order <- c("pdb_all", "pdb_high")
tech_group_order <- c("tech_all", "tech_ct_ai", "tech_ct_noai")
depth_order <- c("depth_all", "depth_foundational", "depth_applied")

indep_vars_order <- c(
  "af_+_ct"
)

coef_order <- c(
  "ct_high", "af_high", "ct", "af"
)

dep_var_labels <- c(
  "ln1p_cited_by_count" = "ln (1 + Cited by count)",
  # "ln1p_cit_0" = "ln (1 + Citations at month < 12)",
  # "ln1p_cit_1" = "ln (1 + Citations at month < 24)",
  "ln1p_fwci" = "ln (1 + Field-Weighted Citation Impact)",
  "ln1p_cit_norm_perc" = "ln (1 + Citation percentile)",
  "ln1p_patent_count" = "ln (1 + Patent count)",
  "ln1p_patent_citation" = "ln (1 + Patent citation count)",
  "ln1p_ca_count" = "ln (1 + CA count)",
  "resolution" = "Resolution"
)

coef_labels <- c(
  "af" = "AlphaFold (All PDB)",
  "ct" = "Counterfactual (All PDB)",
  "af_high" = "AlphaFold (High PDB)",
  "ct_high" = "Counterfactual (High PDB)"
)

# --- Function to generate coefficient plots ---
generate_coef_plots <- function(coef_table) { # nolint

  # rename dep_vars with names
  coef_table <- coef_table %>% # nolint
    mutate( # nolint
      dep_var = recode(dep_var, !!!dep_var_labels), # nolint
      treat_var = ifelse(
        grepl("pdb_high", pdb), # nolint
        paste0(treat_var, "_high"), # nolint
        treat_var
      )
    ) # nolint

  # Iterate over unique field groups
  unique_fields <- unique(coef_table$field)
  unique_dep_vars <- unique(coef_table$dep_var)

  for (single_field in unique_fields) {
    for (single_dep_var in unique_dep_vars) {
      coef_plot_data <- coef_table %>% # nolint
        filter(treat_var %in% names(coef_labels), field == single_field, dep_var == single_dep_var) %>% # nolint
        mutate( # nolint
          depth = factor(gsub("depth_", "", depth), levels = gsub("depth_", "", depth_order)), # nolint
          tech = factor(gsub("tech_", "", tech), levels = gsub("tech_", "", tech_group_order)), # nolint
          tech = recode(tech, "all" = "All Technologies", "ct_ai" = "Counterfactual AI", "ct_noai" = "Counterfactual No AI"), # nolint
          depth = recode(depth, "all" = "All Fields", "foundational" = "Foundational", "applied" = "Applied"), # nolint
          treat_var = factor(
            treat_var,
            levels = coef_order
          ),
          treat_var = recode(treat_var, !!!coef_labels) # nolint
        )

      # Check if coef_plot_data is empty
      if (nrow(coef_plot_data) == 0) {
        message("No data for field: ", single_field, " and dep_var: ", single_dep_var) # nolint
        next
      }

      # Remove "field_" from the title
      field_label <- gsub("field_", "", single_field)

      # Create the plot
      coeffplot <- ggplot( # nolint
        coef_plot_data,
        aes(x = estimate, y = treat_var) # nolint
      ) +
        geom_point( # nolint
          size = 4
        ) +
        geom_errorbarh( # nolint
          aes(xmin = estimate - 1.645 * std_error, xmax = estimate + 1.645 * std_error), # nolint
          height = 0, linewidth = 1 # thicker for 10% significance
        ) +
        geom_errorbarh( # nolint
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.2, linewidth = 0.5 # thinner for 5% significance
        ) +
        geom_hline(yintercept = 2.5, color = "black", linetype = "dashed", linewidth = 1) + # nolint
        geom_hline(yintercept = 5, color = "black", linetype = "dashed", linewidth = 1) + # nolint
        geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
        ggh4x::facet_grid2(depth ~ tech, scales = "free", independent = "x", space = "fixed") + # nolint
        labs( # nolint
          title = paste("Dependent Variable:", single_dep_var),
          subtitle = paste("Field:", field_label), # nolint
          x = "Estimate (with 95% CI)",
          y = "Coefficient Variable"
        ) +
        theme_classic() + # nolint
        theme( # nolint
          axis.text.y = element_text(size = 12), # nolint
          axis.title.x = element_text(size = 14), # nolint
          axis.title.y = element_text(size = 14), # nolint
          strip.text = element_text(size = 12), # nolint
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey"), # nolint
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # nolint
          panel.spacing = unit(2, "lines"), # nolint
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1, "cm") # nolint
        )

      # add counts of obs (pdb all)
      coeffplot <- coeffplot + geom_text( # nolint
        data = coef_plot_data %>% filter(str_detect(treat_var, "All")), # nolint
        aes(label = paste0("n = ", n_obs)), # nolint
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 1.6,
        size = 3, color = "black"
      )

      # add counts of obs (pdb high)
      coeffplot <- coeffplot + geom_text( # nolint
        data = coef_plot_data %>% filter(str_detect(treat_var, "High")), # nolint
        aes(label = paste0("n = ", n_obs)), # nolint
        x = Inf, y = Inf,
        hjust = 1.1, vjust = 26, # Adjusted vjust to move the text down
        size = 3, color = "black"
      )

      # Create the directory if it doesn't exist
      pathdir <- paste0(
        figures, "coef_plot/", single_field, "/"
      )
      outfile <- paste0(pathdir, single_dep_var, ".png")
      message("Saving plot to: ", outfile)
      if (!dir.exists(pathdir)) {
        message("Creating directory: ", pathdir)
        dir.create(pathdir, recursive = TRUE)
      }

      ggsave( # nolint
        outfile, # nolint
        coeffplot,
        width = 20,
        height = 10,
        dpi = 300
      )
    }
  }
}
