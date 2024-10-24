figures <- "data/05_model_output/figures/ecr/"

if (!dir.exists(figures)) {
  dir.create(figures, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# FIGURE GENERATION
# ------------------------------------------------------------------------------

extract_coefficients <- function(results, dep_vars, subsets, cov_sets, fe_list, treat_vars, treat_var_interest = c("is_af")) { # nolint
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

                # Iterate over the list of treatment variables of interest
                for (treat_var_interest_item in treat_var_interest) {
                  if (treat_var_interest_item %in% rownames(coef_info)) {
                    parts <- strsplit(result_name, "__")[[1]]
                    depth <- parts[1]
                    field <- parts[2]
                    dep_var <- parts[3]
                    cov_set <- parts[4]
                    fe <- parts[5]
                    indep_vars <- parts[6]

                    coef_data[[length(coef_data) + 1]] <- data.frame(
                      depth = depth,
                      field = field,
                      treat_var = treat_var_interest_item,
                      dep_var = dep_var,
                      indep_vars = indep_vars,
                      estimate = coef_info[treat_var_interest_item, "Estimate"],
                      std_error = coef_info[treat_var_interest_item, "Std. Error"], # nolint
                      conf_low = coef_info[treat_var_interest_item, "Estimate"] - 1.96 * coef_info[treat_var_interest_item, "Std. Error"], # nolint
                      conf_high = coef_info[treat_var_interest_item, "Estimate"] + 1.96 * coef_info[treat_var_interest_item, "Std. Error"] # nolint
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

# ------------------------------------------------------------------------------
# PLOT GENERATION
# ------------------------------------------------------------------------------

generate_coef_plots <- function(coef_table, interaction_bool = TRUE) { # nolint
  # Filter the coef_table based on the interaction argument
  if (interaction_bool) {
    coef_table <- coef_table[
      coef_table$indep_vars != "is_af_+_is_ct_+_is_af_ct",
    ]
    interaction_folder <- "interacted"
  } else {
    coef_table <- coef_table[
      coef_table$indep_vars == "is_af_+_is_ct_+_is_af_ct",
    ]
    interaction_folder <- "not_interacted"
  }

  # Set the desired order for dependent variables (x-axis) and depth (y-axis)
  dep_var_order <- c(
    "ln1p_cited_by_count", "ln1p_cit_0", "ln1p_cit_1", "fwci",
    "citation_normalized_percentile_value", "ln1p_patent_count",
    "ln1p_patent_citation", "ln1p_ca_count", "resolution", "R_free"
  )

  depth_order <- c("depth_all", "depth_foundational", "depth_applied")

  coef_order <- c(
    "AlphaFold + Counterfactual (int.)",
    "Counterfactual (int. x ext.)", "Counterfactual (int.)",
    "AlphaFold (int. x ext.)", "AlphaFold (int.)"
  )

  # Change coefficient names to more readable labels
  coef_labels <- c(
    "is_afTRUE" = "AlphaFold (int.)",
    "is_afTRUE:af" = "AlphaFold (int. x ext.)",
    "is_ctTRUE" = "Counterfactual (int.)",
    "is_ctTRUE:ct" = "Counterfactual (int. x ext.)",
    "is_af_ctTRUE" = "AlphaFold + Counterfactual (int.)"
  )

  # Split dependent variables into chunks of three
  dep_var_subset <- dep_var_order[dep_var_order %in% coef_table$dep_var]
  dep_var_chunks <- split(
    dep_var_subset, ceiling(seq_along(dep_var_subset) / 3)
  )

  # Iterate over unique field groups
  unique_fields <- unique(coef_table$field)
  for (single_field in unique_fields) {
    for (dep_vars in dep_var_chunks) {
      coef_plot_data <- coef_table %>% # nolint
        filter(treat_var %in% names(coef_labels)) %>% # nolint
        filter(field == single_field) %>% # nolint
        filter(dep_var %in% dep_vars) %>% # nolint
        mutate( # nolint
          dep_var = factor(dep_var, levels = dep_var_order), # nolint
          depth = factor(depth, levels = depth_order), # nolint
          treat_var = factor(
            treat_var,
            levels = names(coef_labels), labels = coef_labels
          )
        )

      # Reorder the treat_var to match the desired order
      coef_plot_data$treat_var <- factor(
        coef_plot_data$treat_var,
        levels = coef_order
      )

      # Create the coefficient plot for the current field group
      coeffplot <- ggplot( # nolint
        coef_plot_data,
        aes(x = estimate, y = treat_var) # nolint
      ) +
        geom_point( # nolint
          aes(color = treat_var),
          size = 3
        ) +
        geom_errorbarh( # nolint
          aes(xmin = conf_low, xmax = conf_high), # nolint
          height = 0.2
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
        facet_grid(depth ~ dep_var, scales = "free", space = "free_x") + # nolint
        labs( # nolint
          title = paste("Coefficient plot for field:", single_field), # nolint
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

      # Create the directory if it doesn't exist
      pathdir <- paste0(
        figures, "coef_plot/", interaction_folder, "/", single_field, "/"
      )
      dep_var_names <- paste(dep_vars, collapse = "_")
      outfile <- paste0(pathdir, dep_var_names, ".png")
      message("Saving plot to: ", outfile)
      if (!dir.exists(pathdir)) {
        message("Creating directory: ", pathdir)
        dir.create(pathdir, recursive = TRUE)
      }

      # if error occurs, skip and continue loop
      tryCatch({
        ggsave( # nolint
          outfile, # nolint
          coeffplot,
          width = 20,
          height = 10,
          dpi = 300
        )
      }, error = function(e) {
        message("Error occurred: ", e)
      })
    }
  }
}