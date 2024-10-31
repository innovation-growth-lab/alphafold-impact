figures <- "data/05_model_output/ecr/figures/"

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
                n_obs <- summary(results[[result_name]])$nobs

                # Iterate over the list of treatment variables of interest
                for (treat_var_interest_item in treat_var_interest) {
                  if (treat_var_interest_item %in% rownames(coef_info)) {
                    parts <- strsplit(result_name, "__")[[1]]
                    depth <- parts[1]
                    field <- parts[2]
                    pdb <- parts[3]
                    dep_var <- parts[4]
                    cov_set <- parts[5]
                    fe <- parts[6]
                    indep_vars <- parts[7]

                    coef_data[[length(coef_data) + 1]] <- data.frame(
                      depth = depth,
                      field = field,
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
pdb_group_order <- c("pdb_All PDB", "pdb_High PDB")
depth_order <- c("All Groups", "Foundational", "Applied")
field_order <- c(
  "field_All Fields",
  "field_Molecular Biology",
  "field_Medicine"
)


coef_order <- c(
  "af:ct_noai", "af:ct_ai", "ct_noai", "ct_ai", "af",
  "af:ct_noai_ind", "af:ct_ai_ind", "ct_noai_ind", "ct_ai_ind", "af_ind",
  "af:ct_noai_high", "af:ct_ai_high", "ct_noai_high", "ct_ai_high", "af_high",
  "af:ct_noai_ind_high", "af:ct_ai_ind_high", "ct_noai_ind_high",
  "ct_ai_ind_high", "af_ind_high"
)

dep_var_labels <- c(
  "ln1p_cited_by_count" = "ln (1 + Cited by count)",
  "ln1p_cit_0" = "ln (1 + Citations at month < 12)",
  "ln1p_cit_1" = "ln (1 + Citations at month < 24)",
  "ln1p_fwci" = "ln (1 + Field-Weighted Citation Impact)",
  "logit_cit_norm_perc" = "logit (Citation percentile)",
  "ln1p_patent_count" = "ln (1 + Patent count)",
  "ln1p_patent_citation" = "ln (1 + Patent citation count)",
  "ln1p_ca_count" = "ln (1 + CA count)",
  "resolution" = "Resolution",
  "R_free" = "R free",
  "pdb_submission" = "PDB submission"
)

coef_labels <- c(
  "af_ind" = "AlphaFold (ext.)",
  "ct_ai_ind" = "Counterfactual AI (ext.)",
  "ct_noai_ind" = "Counterfactual no AI (ext.)",
  "af" = "AlphaFold (int.)",
  "ct_ai" = "Counterfactual AI (int.)",
  "ct_noai" = "Counterfactual no AI (int.)"
)

strip_colors <- c(
  "All Groups - All PDB" = "lightyellow",
  "All Groups - High PDB" = "lightyellow",
  "Foundational - All PDB" = "lightcoral",
  "Foundational - High PDB" = "lightcoral",
  "Applied - All PDB" = "lightblue",
  "Applied - High PDB" = "lightblue"
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
    tryCatch(
      {
        coef_plot_data <- coef_table %>% # nolint
          filter(treat_var %in% names(coef_labels), dep_var == single_dep_var) %>% # nolint
          mutate( # nolint
            depth = factor(gsub("depth_", "", depth), levels = gsub("depth_", "", depth_order)), # nolint
            pdb = factor(gsub("pdb_", "", pdb), levels = gsub("pdb_", "", pdb_group_order)), # nolint
            field = factor(gsub("field_", "", field), levels = gsub("field_", "", field_order)), # nolint
            depth_pdb = factor(paste(depth, pdb, sep = " - "), levels = unique(paste(depth, pdb, sep = " - "))), # nolint
            treat_var = factor(
              treat_var,
              levels = coef_order
            ),
            treat_var = recode(treat_var, !!!coef_labels) # nolint
          )

        # drop depth_pdb with Foundational and High PDB
        coef_plot_data <- coef_plot_data %>% # nolint
          filter(!(pdb == "High PDB")) # nolint

        # Check if coef_plot_data is empty
        if (nrow(coef_plot_data) == 0) {
          message("No data for dep_var: ", single_dep_var) # nolint
          next
        }

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
          geom_hline(yintercept = 3.5, color = "gray", linetype = "dashed", linewidth = 1) + # nolint
          geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
          ggh4x::facet_grid2(
            depth_pdb ~ field,
            scales = "free",
            independent = "x",
            space = "fixed",
            strip = ggh4x::strip_themed(
              background_y = ggh4x::elem_list_rect(
                fill = strip_colors[levels(coef_plot_data$depth_pdb)]
              )
            )
          ) + # nolint
          labs( # nolint
            title = paste("Dependent Variable:", single_dep_var),
            # subtitle = paste("Field:", field_label), # nolint
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

        # add counts of obs (pdb all and ext.)
        coeffplot <- coeffplot + geom_text( # nolint
          data = coef_plot_data %>% filter(str_detect(treat_var, "ext\\.")), # nolint
          aes(label = paste0("n = ", n_obs)), # nolint
          x = Inf, y = Inf,
          hjust = 1.1, vjust = 29.5,
          size = 3, color = "black"
        )

        # Create the directory if it doesn't exist
        pathdir <- paste0(
          figures, "coef_plot/"
        )
        outfile <- paste0(pathdir, single_dep_var, ".png")
        message("Saving plot to: ", outfile)
        if (!dir.exists(pathdir)) {
          message("Creating directory: ", pathdir)
          dir.create(pathdir, recursive = TRUE)
        }

        n_y_facet_rows <- length(unique(coef_plot_data$depth_pdb))
        plot_height <- n_y_facet_rows * 3.5 # nolint

        ggsave( # nolint
          outfile, # nolint
          coeffplot,
          width = 15,
          height = plot_height,
          dpi = 300
        )
      },
      error = function(e) {
        message("Error occurred: ", e)
      }
    )
  }
}
