figures <- "data/05_model_output/authors/nonecr/quarterly/figures/"

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
                    parts <- strsplit(result_name, "__")[[1]] # nolint
                    depth_c <- parts[1]
                    field_c <- parts[2]
                    subgroup_c <- parts[3]
                    dep_var_c <- parts[4]
                    cov_set_c <- parts[5] # nolint
                    fe_c <- parts[6] # nolint
                    indep_vars_c <- parts[7]

                    coef_data[[length(coef_data) + 1]] <- data.frame(
                      depth = depth_c,
                      field = field_c,
                      subgroup = subgroup_c,
                      treat_var = treat_var_interest_item,
                      dep_var = dep_var_c,
                      indep_vars = indep_vars_c,
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
subgroup_order <- c("All PDB", "High PDB", "CEM")
depth_order <- c("All Groups", "Foundational", "Applied")
field_order <- c(
  "field_All Fields",
  "field_Molecular Biology",
  "field_Medicine"
)


coef_order <- c(
  "strong_af:strong_ct_noai", "strong_af:strong_ct_ai", "strong_ct_noai", "strong_ct_ai", "strong_af", # nolint
  "af:ct_noai", "af:ct_ai", "ct_noai", "ct_ai", "af",
  "strong_af:strong_ct_noai_ind", "strong_af:strong_ct_ai_ind", "strong_ct_noai_ind", "strong_ct_ai_ind", "strong_af_ind", # nolint
  "af:ct_noai_ind", "af:ct_ai_ind", "ct_noai_ind", "ct_ai_ind", "af_ind"
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
  "strong_af_ind" = "AlphaFold - Method (ext.)",
  "strong_ct_ai_ind" = "Counterfactual AI - Method (ext.)",
  "strong_ct_noai_ind" = "Counterfactual no AI - Method (ext.)",
  "af" = "AlphaFold (int.)",
  "ct_ai" = "Counterfactual AI (int.)",
  "ct_noai" = "Counterfactual no AI (int.)",
  "strong_af" = "AlphaFold - Method (int.)",
  "strong_ct_ai" = "Counterfactual AI - Method (int.)",
  "strong_ct_noai" = "Counterfactual no AI - Method (int.)"
)

strip_colors <- c(
  "All Groups - All PDB" = "lightyellow",
  "All Groups - High PDB" = "lightyellow",
  "All Groups - CEM" = "lightyellow",
  "Foundational - All PDB" = "lightcoral",
  "Foundational - High PDB" = "lightcoral",
  "Foundational - CEM" = "lightcoral",
  "Applied - All PDB" = "lightblue",
  "Applied - High PDB" = "lightblue",
  "Applied - CEM" = "lightblue"
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
    for (single_subgroup in subgroup_order) {
      tryCatch(
        {
          coef_plot_data <- coef_table %>% # nolint
            filter(treat_var %in% names(coef_labels), dep_var == single_dep_var, subgroup == paste0("subgroup_", single_subgroup, sep="")) %>% # nolint
            mutate( # nolint
              depth = factor(gsub("depth_", "", depth), levels = gsub("depth_", "", depth_order)), # nolint
              field = factor(gsub("field_", "", field), levels = gsub("field_", "", field_order)), # nolint
              subgroup = factor(gsub("subgroup_", "", subgroup), levels = gsub("subgroup_", "", subgroup_order)), # nolint
              depth_subgroup = factor(paste(depth, subgroup, sep = " - "), levels = unique(paste(depth, subgroup, sep = " - "))), # nolint
              treat_var = factor(
                treat_var,
                levels = coef_order
              ),
              treat_var = recode(treat_var, !!!coef_labels) # nolint
            )

          # drop levels in column depth_subgroup if they have zero values
          present_levels <- levels(
            coef_plot_data$depth_subgroup
          )[levels(coef_plot_data$depth_subgroup) %in%
            unique(coef_plot_data$depth_subgroup[coef_plot_data$depth_subgroup != ""])] # nolint
          present_strip_colors <- strip_colors[present_levels]

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
            geom_hline(yintercept = 6.5, color = "gray", linetype = "dashed", linewidth = 1) + # nolint
            geom_vline(xintercept = 0, color = "black", linewidth = 1) + # nolint
            ggh4x::facet_grid2(
              depth_subgroup ~ field,
              scales = "free",
              independent = "x",
              space = "fixed",
              strip = ggh4x::strip_themed(
                background_y = ggh4x::elem_list_rect(
                  fill = present_strip_colors
                )
              )
            ) + # nolint
            labs( # nolint
              title = paste("Dependent Variable:", single_dep_var, " - ", single_subgroup), # nolint
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

          # add counts of obs (subgroup all and ext.)
          coeffplot <- coeffplot + geom_text( # nolint
            data = coef_plot_data %>% filter(str_detect(treat_var, "ext\\.")), # nolint
            aes(label = paste0("n = ", n_obs)), # nolint
            x = Inf, y = Inf,
            hjust = 1.1, vjust = 39.5,
            size = 3, color = "black"
          )

          # Create the directory if it doesn't exist
          pathdir <- paste0(
            figures, "coef_plot/", single_subgroup, "/"
          )
          outfile <- paste0(pathdir, single_dep_var, ".png")
          message("Saving plot to: ", outfile)
          if (!dir.exists(pathdir)) {
            message("Creating directory: ", pathdir)
            dir.create(pathdir, recursive = TRUE)
          }

          n_y_facet_rows <- length(unique(coef_plot_data$depth_subgroup))
          plot_height <- n_y_facet_rows * 4.5 # nolint

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
}
