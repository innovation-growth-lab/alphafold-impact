options(scipen = 999)
tables <- "data/05_model_output/papers/tables/"
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
  "logit_cit_norm_perc" = list(
    file_name = "logit_cit_norm_perc.tex"
  ),
  "patent_count" = list(
    file_name = "patent_count.tex"
  ),
  "patent_citation" = list(
    file_name = "patent_citation.tex"
  ),
  "ca_count" = list(
    file_name = "ca_count.tex"
  ),
  "ln1p_resolution" = list(
    file_name = "ln1p_resolution.tex"
  ),
  "ln1p_R_free" = list(
    file_name = "ln1p_R_free.tex"
  ),
  "num_pdb_ids" = list(
    file_name = "num_pdb_ids.tex"
  ),
  "num_publications" = list(
    file_name = "num_publications.tex"
  ),
  "num_uniprot_structures" = list(
    file_name = "num_uniprot_structures.tex"
  ),
  "num_primary_submissions" = list(
    file_name = "num_primary_submissions.tex"
  ),
  "organism_rarity_mean" = list(
    file_name = "organism_rarity_mean.tex"
  ),
  "mean_tmscore" = list(
    file_name = "mean_tmscore.tex"
  ),
  "num_uniprot_structures_w_disease" = list(
    file_name = "num_uniprot_structures_w_disease.tex"
  ),
  "num_primary_submissions_w_disease" = list(
    file_name = "num_primary_submissions_w_disease.tex"
  ),
  "num_uniprot_structures_w_rare_organisms" = list(
    file_name = "num_uniprot_structures_w_rare_organisms.tex"
  ),
  "num_primary_submissions_w_rare_organisms" = list(
    file_name = "num_primary_submissions_w_rare_organisms.tex"
  ),
  "num_uniprot_structures_w_low_similarity" = list(
    file_name = "num_uniprot_structures_w_low_similarity.tex"
  ),
  "num_primary_submissions_w_low_similarity" = list(
    file_name = "num_primary_submissions_w_low_similarity.tex"
  ),
  "num_diseases" = list(
    file_name = "num_diseases.tex"
  ),
  "mesh_C" = list(
    file_name = "mesh_C.tex"
  )
)

dict_vars <- c(
  "af" = "AlphaFold",
  "ct_ai" = "Counterfactual AI",
  "ct_noai" = "Counterfactual No AI",
  "strong1" = "Method"
)

# Function to generate tables

# Function to generate tables
generate_tables <- function(results, dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint

  depths <- c("depth_All Groups", "depth_Foundational", "depth_Applied")
  fields <- c(
    "field_All Fields",
    "field_Molecular Biology",
    "field_Medicine"
  )
  field_labels <- gsub("field_", "", fields)
  subgroups <- c("subgroup_All PDB", "subgroup_High PDB", "subgroup_CEM")

  for (dep_var in dep_vars) {
    file_name <- table_info[[dep_var]]$file_name

    # Iterate over subsets
    for (depth in depths) {
      result_names <- c()
      for (field in fields) {
        for (subgroup in subgroups) {
          # Iterate over covariate sets, fixed effects, and treatment variables # nolint
          for (cov_set in cov_sets) {
            for (fe in fe_list) {
              for (treat_var in treat_vars) {
                # Build the result name
                result_name <- paste0(
                  depth, "__", field, "__", subgroup, "__", # nolint
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
      }

      # Use etable to output the table for each depth-field combination
      if (length(result_names) > 0) {
        message(
          paste0(
            "Generating table for depth: ", depth,
            ", pdb group: ", subgroup,
            " and ", dep_var
          )
        )
        pathdir <- paste0(
          tables,
          depth, "/"
        )
        if (!dir.exists(pathdir)) {
          dir.create(pathdir, recursive = TRUE)
        }

        # Generate the etable output
        etable_output <- fixest::etable(
          results[result_names],
          drop = c("num_publications", "Constant"),
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

        # Add tech_group headers and \cmidrule after row 5
        field_headers <- paste0(
          "\\multicolumn{3}{c}{", field_labels, "}"
        )
        field_headers <- paste0(
          " & ", paste(field_headers, collapse = " & "), " \\\\"
        )

        field_cmidrules <- paste0(
          "\\cmidrule(lr){",
          seq(2, length(field_labels) * 3 + 1, by = 3), "-",
          seq(4, length(field_labels) * 3 + 1, by = 3), "}"
        )
        field_cmidrules <- paste0(
          paste(field_cmidrules, collapse = " ")
        )

        # subgroup headers
        subgroup_headers <- paste0(
          rep(
            "\\multicolumn{1}{c}{All PDB} & \\multicolumn{1}{c}{High PDB} & \\multicolumn{1}{c}{CEM}", # nolint
            length(field_labels)
          )
        )

        subgroup_headers <- paste0(
          " & ", paste(subgroup_headers, collapse = " & "), " \\\\"
        )

        subgroup_cmidrules <- paste0(
          "\\cmidrule(lr){", seq(2, length(field_labels) * 3 + 1, by = 2), "-", # nolint
          seq(2, length(field_labels) * 3 + 1, by = 2), "} \\cmidrule(lr){",
          seq(3, length(field_labels) * 3 + 1, by = 2), "-",
          seq(3, length(field_labels) * 3 + 1, by = 2), "}"
        )
        subgroup_cmidrules <- paste0(
          "\\cmidrule(lr){1-1}",
          paste(subgroup_cmidrules, collapse = " ")
        )

        etable_lines <- unlist(strsplit(etable_output, "\n"))

        # Insert tech_group headers after row 5
        etable_lines <- append(etable_lines, field_headers, after = 5)
        etable_lines <- append(etable_lines, field_cmidrules, after = 6)

        # Insert subgroup headers after row 6
        etable_lines <- append(etable_lines, subgroup_headers, after = 7)
        etable_lines <- append(etable_lines, subgroup_cmidrules, after = 8)

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
