options(scipen = 999)
tables <- "data/05_model_output/nonecr/tables/"
if (!dir.exists(tables)) {
  dir.create(tables, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------

table_info <- list(
  "ln1p_cited_by_count" = list(
    file_name = "ln1p_cited_by_count.tex"
  ),
  "ln1p_cit_0" = list(
    file_name = "ln1p_cit_0.tex"
  ),
  "ln1p_cit_1" = list(
    file_name = "ln1p_cit_1.tex"
  ),
  "ln1p_fwci" = list(
    file_name = "ln1p_fwci.tex"
  ),
  "ln1p_cit_norm_perc" = list(
    file_name = "ln1p_cit_norm_perc.tex"
  ),
  "ln1p_patent_count" = list(
    file_name = "ln1p_patent_count.tex"
  ),
  "ln1p_patent_citation" = list(
    file_name = "ln1p_patent_citation.tex"
  ),
  "ln1p_ca_count" = list(
    file_name = "ln1p_ca_count.tex"
  ),
  "resolution" = list(
    file_name = "resolution.tex"
  ),
  "R_free" = list(
    file_name = "R_free.tex"
  )
)

dict_vars <- c(
  "af" = "AlphaFold",
  "ct" = "Counterfactual"
)

# Function to generate tables
generate_tables <- function(results, dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint

  depths <- c()
  fields <- c()

  # Iterate over subsets to extract depths and fields
  for (sub in subsets) {
    parts <- strsplit(sub, "__")[[1]]
    depths <- c(depths, parts[1])
    fields <- c(fields, parts[2])
  }

  # Get unique depths and fields
  unique_depths <- unique(depths)
  unique_fields <- unique(fields)

  # Generate all possible depth-field combinations
  depth_field_pairs <- expand.grid(unique_depths, unique_fields)
  pdb_groups <- c("pdb_all", "pdb_high")

  for (dep_var in dep_vars) {
    file_name <- table_info[[dep_var]]$file_name

    # Iterate over subsets
    for (i in seq_len(nrow(depth_field_pairs))) {
      depth <- depth_field_pairs[i, 1]
      field <- depth_field_pairs[i, 2]

      result_names <- c()
      for (pdb_group in pdb_groups) {
        # Iterate over covariate sets, fixed effects, and treatment variables # nolint
        for (cov_set in cov_sets) {
          for (fe in fe_list) {
            for (treat_var in treat_vars) {
              # Build the result name
              result_name <- paste0(
                depth, "__", field, "__", pdb_group, "__",
                dep_var, "__", cov_set, "__", fe, "__", gsub(" ", "_", treat_var) # nolint
              )
              # Check if result exists
              if (result_name %in% names(results)) {
                result_names <- c(result_names, result_name)
              }
            }
          }
        }
      }

      # Use etable to output the table for each depth-field combination
      if (length(result_names) > 0) {
        message(
          paste0(
            "Generating table for depth: ", depth,
            ", field: ", field,
            " and ", dep_var
          )
        )
        pathdir <- paste0(
          tables,
          depth, "/",
          field, "/"
        )
        if (!dir.exists(pathdir)) {
          dir.create(pathdir, recursive = TRUE)
        }

        # Generate the etable output
        etable_output <- fixest::etable(
          results[result_names],
          drop = c("num_publications"),
          tex = TRUE,
          dict = dict_vars,
          digits = 3,
          digits.stats = 2,
          powerBelow = -20,
          fitstat = c("n", "r2")
        )

        # force mean y to appear
        mean_y_values <- sapply(results[result_names], function(model) {
          mean(model$fitted.values + model$residuals)
        })

        mean_y_row <- paste0(
          "Mean(Dep. Var.) & ", paste(sprintf("%.3f", mean_y_values), collapse = " & "), " \\\\" # nolint
        )

        # tech_groups_latex <- c(
        #   "All Technologies",
        #   "Counterfactual AI",
        #   "Counterfactual No AI"
        # )

        # # Add tech_group headers and \cmidrule after row 5
        # tech_group_headers <- paste0(
        #   "\\multicolumn{4}{c}{", tech_groups_latex, "}"
        # )
        # tech_group_headers <- paste0(
        #   " & ", paste(tech_group_headers, collapse = " & "), " \\\\"
        # )

        # tech_group_cmidrules <- paste0(
        #   "\\cmidrule(lr){",
        #   seq(2, length(tech_groups_latex) * 4 + 1, by = 4), "-",
        #   seq(5, length(tech_groups_latex) * 4 + 1, by = 4), "}"
        # )
        # tech_group_cmidrules <- paste0(
        #   paste(tech_group_cmidrules, collapse = " ")
        # )

        # pdb headers
        pdb_headers <- "\\multicolumn{2}{c}{All PDB} & \\multicolumn{2}{c}{High PDB} \\\\"
        pdb_cmidrules <- "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}"

        coefficient_headers <- paste0(
          rep(
            "\\multicolumn{1}{c}{Extensive} & \\multicolumn{1}{c}{Intensive}",
            2
          )
        )

        coefficient_headers <- paste0(
          "Variables & ", paste(coefficient_headers, collapse = " & "), " \\\\"
        )

        coefficient_cmidrules <- paste0(
          "\\cmidrule(lr){", seq(2, 5, by = 2), "-",
          seq(2, 5, by = 2), "} \\cmidrule(lr){",
          seq(3, 5, by = 2), "-",
          seq(3, 5, by = 2), "}"
        )
        coefficient_cmidrules <- paste0(
          "\\cmidrule(lr){1-1}",
          paste(coefficient_cmidrules, collapse = " ")
        )

        # Split the etable output into lines
        etable_lines <- unlist(strsplit(etable_output, "\n"))

        # # Insert tech_group headers after row 5
        # etable_lines <- append(etable_lines, tech_group_headers, after = 5)
        # etable_lines <- append(etable_lines, tech_group_cmidrules, after = 6)

        # Insert pdb headers after row 6
        etable_lines <- append(etable_lines, pdb_headers, after = 5)
        etable_lines <- append(etable_lines, pdb_cmidrules, after = 6)

        # Insert coefficient headers after row 8
        etable_lines <- append(etable_lines, coefficient_headers, after = 7)
        etable_lines <- append(etable_lines, coefficient_cmidrules, after = 8)

        # drop lines 10-11
        etable_lines <- etable_lines[-c(10, 11, 12)]

        # add mean y row in 6th last row
        etable_lines <- append(
          etable_lines, mean_y_row,
          after = length(etable_lines) - 5
        )

        # Combine the lines back into a single string
        etable_output <- paste(etable_lines, collapse = "\n")

        # Write the table to a file
        writeLines(
          etable_output,
          con = paste0(pathdir, file_name)
        )
      }
    }
  }
}
