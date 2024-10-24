tables <- "data/05_model_output/ecr/tables/"
if (!dir.exists(tables)) {
  dir.create(tables, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------

variable_interest <- c(
  "is_af", "is_af:af", "is_ct", "is_ct:ct", "is_af_ct", "is_af_ct:af_ct"
)

table_info <- list(
  "ln1p_cited_by_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cited_by_count.tex"
  ),
  "ln1p_cit_0" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cit_0.tex"
  ),
  "ln1p_cit_1" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_cit_1.tex"
  ),
  "fwci" = list(
    vars_to_keep = variable_interest,
    file_name = "fwci.tex"
  ),
  "citation_normalized_percentile_value" = list(
    vars_to_keep = variable_interest,
    file_name = "citation_normalized_percentile_value.tex"
  ),
  "ln1p_patent_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_patent_count.tex"
  ),
  "ln1p_patent_citation" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_patent_citation.tex"
  ),
  "ln1p_ca_count" = list(
    vars_to_keep = variable_interest,
    file_name = "ln1p_ca_count.tex"
  ),
  "resolution" = list(
    vars_to_keep = variable_interest,
    file_name = "resolution.tex"
  ),
  "R_free" = list(
    vars_to_keep = variable_interest,
    file_name = "R_free.tex"
  )
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