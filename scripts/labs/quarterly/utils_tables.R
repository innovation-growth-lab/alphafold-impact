options(scipen = 999)
tables <- "data/05_model_output/labs/quarterly/tables/"
if (!dir.exists(tables)) {
  dir.create(tables, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# TABLE GENERATION
# ------------------------------------------------------------------------------

table_info <- list(
  "cited_by_count" = list(
    file_name = "cited_by_count.tex"
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
  "num_pdb_submissions" = list(
    file_name = "num_pdb_submissions.tex"
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
  "ln1p_organism_rarity_mean" = list(
    file_name = "ln1p_organism_rarity_mean.tex"
  ),
  "ln1p_organism_rarity_max" = list(
    file_name = "ln1p_organism_rarity_max.tex"
  ),
  "ln1p_max_tmscore" = list(
    file_name = "ln1p_max_tmscore.tex"
  ),
  "ln1p_max_fident" = list(
    file_name = "ln1p_max_fident.tex"
  ),
  "ln1p_max_score" = list(
    file_name = "ln1p_max_score.tex"
  ),
  "normalised_max_tmscore" = list(
    file_name = "normalised_max_tmscore.tex"
  ),
  "normalised_max_fident" = list(
    file_name = "normalised_max_fident.tex"
  ),
  "normalised_max_score" = list(
    file_name = "normalised_max_score.tex"
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
  "ln1p_mesh_C" = list(
    file_name = "ln1p_mesh_C.tex"
  )
)

dict_vars <- c(
  "af" = "AlphaFold",
  "ct_ai" = "AI Frontier",
  "ct_pp" = "PP Frontier",
  "ct_sb" = "SB Frontier",
  "af_intent_strong" = "AlphaFold - Method",
  "af_intent_weak" = "AlphaFold - Bkg.",
  "af_intent_mixed" = "AlphaFold - Mixed",
  "ct_ai_intent_strong" = "AI Frontier - Method",
  "ct_ai_intent_weak" = "AI Frontier - Bkg.",
  "ct_ai_intent_mixed" = "AI Frontier - Mixed",
  "ct_pp_intent_strong" = "PP Frontier - Method",
  "ct_pp_intent_weak" = "PP Frontier - Bkg.",
  "ct_pp_intent_mixed" = "PP Frontier - Mixed",
  "ct_sb_intent_strong" = "SB Frontier - Method",
  "ct_sb_intent_weak" = "SB Frontier - Bkg.",
  "ct_sb_intent_mixed" = "SB Frontier - Mixed"
)

# Function to generate tables
generate_tables <- function(results, dep_vars, table_info, subsets, cov_sets, fe_list, treat_vars) { # nolint

  scopes <- c("scope_All", "scope_Intent")
  fields <- c(
    "field_All Fields",
    "field_Molecular Biology",
    "field_Medicine"
  )
  subgroups <- c("subgroup_All PDB", "subgroup_High PDB")

  for (dep_var in dep_vars) {
    file_name <- table_info[[dep_var]]$file_name

    # Iterate over subsets
    latex_output <- list()
    for (field in fields) {
      field_label <- gsub("field_", "", field)
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

        # Generate the etable output
        etable_output <- fixest::etable(
          results[result_names],
          drop = c(
            # Drop all interaction terms (both colon and underscore format)
            ".*:.*",
            ".*_ct_.*",
            # Drop specific patterns we don't want
            "^af_ct.*",
            "^ct_ai_ct.*",
            # Drop renamed interaction terms
            ".*\\$\\\\times\\$.*", # matches the LaTeX formatted interactions
            "AlphaFold.*AI Frontier.*",
            "AlphaFold.*No AI Frontier.*",
            "AI Frontier.*No AI Frontier.*",
            # Drop any other controls
            "num_publications",
            "Constant",
            "^covid_",
            "^field_",
            "^mesh_",
            "^institution_"
          ),
          tex = TRUE,
          dict = dict_vars,
          digits = 3,
          digits.stats = 2,
          powerBelow = -20,
          fitstat = c("n", "pr2", "r2")
        )

        # force mean y to appear
        mean_y_values <- sapply(results[result_names], function(model) {
          mean(model$fitted.values + model$residuals)
        })

        mean_y_row <- paste0(
          "Mean(Dep. Var.) & ", paste(sprintf("%.3f", mean_y_values), collapse = " & "), " \\\\" # nolint
        )

        etable_lines <- unlist(strsplit(etable_output, "\n"))

        # Top Field: All Fields
        if (field_label == "All Fields") {
          pdb_headers <- c(
            "\\multicolumn{3}{c}{All Authors} & \\multicolumn{3}{c}{Experienced Authors}" # nolint
          )
          pdb_headers <- paste0(
            " & ", paste(pdb_headers, collapse = " & "), " \\\\"
          )
          pdb_cmidrules <- "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}"

          # Second level: Intent headers
          intent_headers <- paste0(
            rep(
              "\\multicolumn{1}{c}{All} & \\multicolumn{2}{c}{With Intent Data}", # nolint
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

          # add field label in italics and with midrules
          field_label_line <- paste0(
            "\\multicolumn{6}{c}{\\textit{", field_label, "}} \\\\"
          )
          field_label_line <- paste0(
            " & ", paste(field_label_line, collapse = " & "), " \\\\"
          )

          # Insert headers and cmidrules
          etable_lines <- append(etable_lines, pdb_headers, after = 5)
          etable_lines <- append(etable_lines, pdb_cmidrules, after = 6)
          etable_lines <- append(etable_lines, intent_headers, after = 7)
          etable_lines <- append(etable_lines, intent_cmidrules, after = 8)
          etable_lines <- append(etable_lines, field_label_line, after = 9)
          # drop old header lines
          etable_lines <- etable_lines[-c(11, 12, 13)]


          # add mean y row in 6th last row
          etable_lines <- append(
            etable_lines, mean_y_row,
            after = length(etable_lines) - 5
          )

          chunks <- unlist(
            strsplit(paste(etable_lines, collapse = "\n"), "\\\\midrule")
          )
          etable_lines <- unlist(
            strsplit(paste(chunks[c(1, 2, 3, 5)], collapse = "\\midrule"), "\n")
          )
          latex_output <- c(latex_output, etable_lines)


          # Second Field: Molecular Biology
        } else if (field_label == "Molecular Biology") {
          chunks <- unlist(
            strsplit(paste(etable_lines, collapse = "\n"), "\\\\midrule")
          )
          etable_lines <- unlist(
            strsplit(paste(chunks[c(4, 6)], collapse = "\\midrule"), "\n")
          )

          # add field label in italics and with midrules
          field_label_line <- paste0(
            "\\multicolumn{6}{c}{\\textit{", field_label, "}} \\\\"
          )
          field_label_line <- paste0(
            " & ", paste(field_label_line, collapse = " & "), " \\\\"
          )
          # add as new header the field label
          etable_lines <- append(etable_lines, field_label_line, after = 1)

          # remove "Variables" row (3)
          etable_lines <- etable_lines[-c(1, 3)]

          # add mean y row in the last row
          etable_lines <- append(
            etable_lines, mean_y_row,
            after = length(etable_lines)
          )
          latex_output <- c(latex_output, etable_lines)
        } else if (field_label == "Medicine") {
          chunks <- unlist(
            strsplit(paste(etable_lines, collapse = "\n"), "\\\\midrule")
          )
          etable_lines <- unlist(
            strsplit(paste(chunks[c(4, 6, 7, 8)], collapse = "\\midrule"), "\n")
          )

          # add field label in italics and with midrules
          field_label_line <- paste0(
            "\\multicolumn{6}{c}{\\textit{", field_label, "}} \\\\"
          )
          field_label_line <- paste0(
            " & ", paste(field_label_line, collapse = " & "), " \\\\"
          )
          # add as new header the field label
          etable_lines <- append(etable_lines, field_label_line, after = 1)

          # remove "Variables" row (3)
          etable_lines <- etable_lines[-c(1, 3)]

          # add mean y row in the 6th last row
          etable_lines <- append(
            etable_lines, mean_y_row,
            after = length(etable_lines) - 5
          )

          latex_output <- c(latex_output, etable_lines)
        } else {
          # raise error
          stop(paste0("Field label not found: ", field_label))
        }
      }
    }
    # join the latex output into a single file
    latex_output <- paste(latex_output, collapse = "\n")

    pathdir <- paste0(
      tables,
      "/"
    )
    if (!dir.exists(pathdir)) {
      dir.create(pathdir, recursive = TRUE)
    }

    # Write the table to a file
    writeLines(
      latex_output,
      con = paste0(pathdir, file_name)
    )
  }
}
