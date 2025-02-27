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
  "ct_ai" = "AI Frontiers",
  "ct_noai" = "No AI Frontiers",
  "af:strong0" = "AlphaFold - Background",
  "af:strong1" = "AlphaFold - Method",
  "ct_ai:strong0" = "AI Frontiers - Background",
  "ct_ai:strong1" = "AI Frontiers - Method",
  "strong0:ct_ai" = "AI Frontiers - Background",
  "strong1:ct_ai" = "AI Frontiers - Method",
  "ct_noai:strong0" = "No AI Frontiers - Background",
  "ct_noai:strong1" = "No AI Frontiers - Method",
  "strong0:ct_noai" = "No AI Frontiers - Background",
  "strong1:ct_noai" = "No AI Frontiers - Method",
  "strong1" = "Method",
  "strong0" = "Background"
)

# Function to generate tables
generate_tables <- function(results, dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint

  scopes <- c("scope_All", "scope_Intent")
  fields <- c(
    "field_All Fields",
    "field_Molecular Biology",
    "field_Medicine"
  )
  field_labels <- gsub("field_", "", fields)
  subgroups <- c("subgroup_All PDB", "subgroup_High PDB")

  for (dep_var in dep_vars) {
    file_name <- table_info[[dep_var]]$file_name

    # Iterate over subsets
    for (field in fields) {
      result_names <- c()
      for (subgroup in subgroups) {
        for (scope in scopes) {
          # Iterate over covariate sets, fixed effects, and treatment variables # nolint
          for (cov_set in cov_sets) {
            for (fe in fe_list) {
              for (treat_var in treat_vars) {
                # Build the result name
                result_name <- paste0(
                  scope, "__", field, "__", subgroup, "__", # nolint
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


      # Use etable to output the table for each scope-field combination
      if (length(result_names) > 0) {
        message(
          paste0(
            "Generating table for dep_var: ", dep_var
          )
        )
        pathdir <- paste0(
          tables,
          "/", field, "/"
        )
        if (!dir.exists(pathdir)) {
          dir.create(pathdir, recursive = TRUE)
        }

        # Generate the etable output
        etable_output <- fixest::etable(
          results[result_names],
          drop = c(
            "num_publications", "Constant",
            "^covid_", "^field_", "^mesh_", "^institution_"
          ),
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

        # Add headers and cmidrules
        # First level: PDB Experience headers
        pdb_headers <- c(
          "\\multicolumn{3}{c}{All Authors} & \\multicolumn{3}{c}{Experienced Authors}"
        )
        pdb_headers <- paste0(
          " & ", paste(pdb_headers, collapse = " & "), " \\\\"
        )
        pdb_cmidrules <- "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}"

        # Second level: Intent headers
        intent_headers <- paste0(
          rep(
            "\\multicolumn{1}{c}{All} & \\multicolumn{2}{c}{With Intent Data}",
            2
          )
        )
        intent_headers <- paste0(
          " & ", paste(intent_headers, collapse = " & "), " \\\\"
        )
        intent_cmidrules <- paste0(
          "\\cmidrule(lr){2-2} \\cmidrule(lr){3-4}",
          " \\cmidrule(lr){5-5} \\cmidrule(lr){6-7}"
        )

        etable_lines <- unlist(strsplit(etable_output, "\n"))

        # Insert headers and cmidrules
        etable_lines <- append(etable_lines, pdb_headers, after = 5)
        etable_lines <- append(etable_lines, pdb_cmidrules, after = 6)
        etable_lines <- append(etable_lines, intent_headers, after = 7)
        etable_lines <- append(etable_lines, intent_cmidrules, after = 8)

        # drop old header lines
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
