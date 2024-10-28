tables <- "data/05_model_output/papers/tables/"
if (!dir.exists(tables)) {
  dir.create(tables, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------

# Define variables of interest
variable_interest <- c(
  "treatment_af_dyn",
  "treatment_af_dyn:strong",
  "treatment_af_dyn:high_pdb",
  "treatment_af_dyn:strong:high_pdb"
)

# Define mapping from dependent variables to variables of interest and file names # nolint
table_info <- list(
  "cited_by_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "03_cited_by_count_ln_productivity_lvl0.tex"
  ),
  "patent_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "01_patent_count_ln_productivity_lvl0.tex"
  ),
  "patent_citation_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "02_patent_citation_ln_productivity_lvl0.tex"
  ),
  "ca_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "04_ca_count_ln_productivity_lvl0.tex"
  )
)

# Define subsets, covariate sets, fixed effects
subsets_lvl0 <- c(
  "all_lvl0", "all_lvl0_ct", "all_lvl0_w_high_pdb", "af_ct_lvl0_w_high_pdb"
)

# Define mapping from dependent variables to variables
table_info_lvl12 <- list(
  "patent_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "01_patent_count_ln_productivity_lvl12.tex"
  ),
  "cited_by_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "03_cited_by_count_ln_productivity_lvl12.tex"
  ),
  "patent_citation_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "02_patent_citation_ln_productivity_lvl12.tex"
  ),
  "ca_count_ln" = list(
    vars_to_keep = variable_interest,
    file_name = "04_ca_count_ln_productivity_lvl12.tex"
  )
)

# Define subsets, covariate sets, fixed effects
subsets_lvl12 <- c(
  "all_lvl12", "all_lvl12_ct", "all_lvl12_w_high_pdb", "af_ct_lvl12_w_high_pdb"
)

cov_sets <- c("base2")
fe_list <- c("fe1")
treat_vars <- c(
  "treatment_af_dyn",
  "treatment_af_dyn + treatment_af_dyn:strong + treatment_af_dyn:high_pdb + treatment_af_dyn:strong:high_pdb" # nolint
)


# Function to generate tables
generate_tables <- function(results, dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint
  for (dep_var in dep_vars) {
    vars_to_keep <- table_info[[dep_var]]$vars_to_keep
    file_name <- table_info[[dep_var]]$file_name

    # Iterate over subsets
    for (sub in subsets) {
      result_names <- c()

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
              result_names <- c(result_names, result_name)
            }
          }
        }
      }

      # Use etable to output the table for each subset
      if (length(result_names) > 0) {
        # split sub on __
        parts <- strsplit(sub, "__")[[1]]
        depth_path <- parts[1]
        field_path <- parts[2]
        pathdir <- paste0(tables, depth_path, "/", field_path, "/")
        if (!dir.exists(pathdir)) {
          dir.create(pathdir, recursive = TRUE)
        }
        fixest::etable(
          results[result_names],
          keep = vars_to_keep,
          file = paste0(pathdir, file_name)
        )
      }
    }
  }
}