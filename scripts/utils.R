
generate_summary_table <- function(df, output_label) {
  sb_data_qtly_t0 <- df %>%
    filter(time == 0)
  # Specify columns for which to calculate summary statistics
  cols <- c(
    "num_publications",
    "cited_by_count",
    "ct0",
    "ct1",
    # "ct2",
    "pdb_share",
    "strong",
    "protein_share",
    "experimental_share",
    "patent_count",
    "field_biochemistry_genetics_and_molecular_biology",
    "mesh_C"
  )

  # Convert 'strong' column to numeric
  sb_data_qtly_t0$strong <- as.integer(sb_data_qtly_t0$strong) - 1
  options(scipen = 999)

  # Calculate summary statistics for each pi_id and seed_type
  summary_means <- sb_data_qtly_t0 %>%
    group_by(seed) %>%
    summarise(
      across(all_of(cols), list(mean = ~ mean(., na.rm = TRUE)),
        .names = "{.col}"
      ),
      nobs = n()
    ) %>%
    mutate_at(vars(
      num_publications, ct0, ct1, strong,
      protein_share, experimental_share, patent_count
    ), ~ round(., 2)) %>%
    mutate_at(vars(pdb_share), ~ format(round(., 5), scientific = TRUE)) %>%
    mutate_at(
      vars(
        field_biochemistry_genetics_and_molecular_biology,
        mesh_C
      ),
      ~ round(., 3)
    ) %>%
    mutate_at(vars(cited_by_count, nobs), ~ round(., 0)) %>%
    mutate(across(where(is.numeric), ~ format(., scientific = FALSE))) %>%
    pivot_longer(cols = -seed, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = seed, values_from = value) %>%
    rename_with(~ paste0(gsub(",", "_", .), "_mean"))

  summary_sds <- sb_data_qtly_t0 %>%
    group_by(seed) %>%
    summarise(
      across(all_of(cols), list(sd = ~ sd(., na.rm = TRUE)),
        .names = "{.col}"
      ),
      nobs = n()
    ) %>%
    mutate_at(vars(
      num_publications, ct0, ct1, strong,
      protein_share, experimental_share, patent_count
    ), ~ round(., 2)) %>%
    mutate_at(vars(pdb_share), ~ format(round(., 5), scientific = TRUE)) %>%
    mutate_at(
      vars(
        field_biochemistry_genetics_and_molecular_biology,
        mesh_C
      ),
      ~ round(., 3)
    ) %>%
    mutate_at(vars(cited_by_count, nobs), ~ round(., 0)) %>%
    mutate(across(where(is.numeric), ~ format(., scientific = FALSE))) %>%
    pivot_longer(cols = -seed, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = seed, values_from = value) %>%
    rename_with(~ paste0(gsub(",", "_", .), "_sd"))





  # Bind the columns
  tab1 <- bind_cols(
    summary_means$alphafold_mean, summary_sds$alphafold_sd,
    summary_means$alphafold_ct_mean, summary_sds$alphafold_ct_sd,
    summary_means$ct_mean, summary_sds$ct_sd,
    summary_means$other_mean, summary_sds$other_sd
  )

  # Add the 'names' column
  tab1 <- tab1 %>%
    mutate(names = c(
      "Number of Publications", "Number of Citations",
      "Citations - 12 months", "Citations - 24 months", "PDB Share",
      "Strong Chain", "Protein Prediction Share",
      "Experimental Share", "Patent Citation Count",
      "Biochem, Genetics, Molecular Share",
      "C-class MeSH Share", "Number of Observations"
    ))

  # Select the columns
  tab1 <- tab1 %>%
    select(names, everything())

  otab1 <- tab1 %>%
    kable("latex",
      align = "lcccccccccc", booktabs = TRUE, linesep = c("", "", ""),
      col.names = NULL,
      escape = FALSE,
      label = "table1",
      caption = "Descriptive Statistics"
    ) %>%
    kable_styling(
      position = "center", latex_options = c("HOLD_position", "scale_down")
    ) %>%
    add_header_above(
      c(" ", "Mean", "SD", "Mean", "SD", "Mean", "SD", "Mean", "SD")
    ) %>%
    add_header_above(
      c(
        "\\\\textit{Variable}" = 1,
        "Alphafold" = 2,
        "Alphafold & Other Tech" = 2,
        "Other Tech" = 2,
        "Broader SB" = 2
      ),
      escape = FALSE
    ) %>%
    row_spec(12, extra_latex_after = "\\midrule") %>%
    footnote(
      general_title = "",
      general = c(
        "\\\\footnotesize \\\\textit{Note:} Lorem Ipsum"
      ),
      threeparttable = TRUE,
      footnote_as_chunk = TRUE,
      escape = FALSE
    )

  write_lines(otab1, file = paste(tables, output_label, sep = ""))
}