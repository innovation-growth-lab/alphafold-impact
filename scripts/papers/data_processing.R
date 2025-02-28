# Clean out the workspace
rm(list = ls())
options(max.print = 1000)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "aws.s3", "yaml", "zoo", "MatchIt"
)
new_packages <- list_of_packages[
  !(list_of_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(
    new_packages,
    repos = "http://cran.us.r-project.org"
  )
}
invisible(lapply(list_of_packages, library, character.only = TRUE))

# Set working directory
setwd("~/projects/alphafold-impact/")
pathdir <- "data/05_model_output/papers/data/"

# Create directories if they do not exist
if (!dir.exists(pathdir)) {
  dir.create(pathdir, recursive = TRUE)
}

# Assign commonly used dplyr functions
select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize
bind_rows <- dplyr::bind_rows

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

credentials <- yaml.load_file("conf/base/credentials.yml")

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key,
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret,
  "AWS_DEFAULT_REGION" = "eu-west-2"
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
path <- "oct/04_output/publications/regression_inputs.parquet" # nolint


# Fetch the data from the S3 bucket
papers <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# Replace commas in column names
colnames(papers) <- gsub(",", "", colnames(papers))

# creating quantiles for institution controls and translational
papers <- papers %>%
  mutate(
    institution_2yr_mean_citedness = factor(
      ntile("2yr_mean_citedness", 4)
    ),
    institution_h_index = factor(ntile(h_index, 4)),
    institution_i10_index = factor(ntile(i10_index, 4)),
    institution_cited_by_count = factor(
      ntile(cited_by_count_institution, 4)
    ),
    organism_rarity_mean_quantile = factor(ntile(organism_rarity_mean, 4)),
    organism_rarity_max_quantile = factor(ntile(organism_rarity_max, 4)),
    mean_tmscore_quantile = factor(ntile(mean_tmscore, 4)),
    max_tmscore_quantile = factor(ntile(max_tmscore, 4))
  )

# for duplicate ids, keep first. fill few nas
papers <- papers %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    institution_type = ifelse(
      is.na(type), "unknown", type
    ),
    institution_country_code = ifelse(
      is.na(country_code), "unknown", country_code
    ),
    institution_2yr_mean_citedness = ifelse(
      is.na(institution_2yr_mean_citedness), "unknown",
      institution_2yr_mean_citedness
    ),
    institution_h_index = ifelse(
      is.na(institution_h_index), "unknown", institution_h_index
    ),
    institution_i10_index = ifelse(
      is.na(institution_i10_index), "unknown", institution_i10_index
    ),
    institution_cited_by_count = ifelse(
      is.na(institution_cited_by_count), "unknown",
      institution_cited_by_count
    ),
    pdb_submission = ifelse(pdb_submission > 0, 1, 0),
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    fwci = ifelse(is.na(fwci), 0, fwci),
    num_uniprot_structures = ifelse(
      is.na(num_uniprot_structures), 0, num_uniprot_structures
    ),
    num_pdb_ids = ifelse(is.na(num_pdb_ids), 0, num_pdb_ids),
    num_primary_submissions = ifelse(is.na(num_primary_submissions), 0, num_primary_submissions), # nolint
  )

# adding translational variables
papers <- papers %>%
  mutate(
    num_uniprot_structures_w_disease = ifelse(
      num_uniprot_structures > 0 & num_diseases > 0, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_disease = ifelse(
      num_primary_submissions > 0 & num_diseases > 0, num_primary_submissions, 0 # nolint
    ),
    num_uniprot_structures_w_rare_organisms = ifelse(
      num_uniprot_structures > 0 & organism_rarity_mean_quantile == 4, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_rare_organisms = ifelse(
      num_primary_submissions > 0 & organism_rarity_mean_quantile == 4, num_primary_submissions, 0 # nolint
    ),
    num_uniprot_structures_w_low_similarity = ifelse(
      num_uniprot_structures > 0 & mean_tmscore_quantile == 1, num_uniprot_structures, 0 # nolint
    ),
    num_primary_submissions_w_low_similarity = ifelse(
      num_primary_submissions > 0 & mean_tmscore_quantile == 1, num_primary_submissions, 0 # nolint
    )
  )

# Create 'strong' variable
papers <- papers %>%
  mutate(
    strong = as.factor(case_when(
      chain_label %in% c("strong", "partial_strong") ~ 1,
      chain_label %in% c("weak", "partial_weak") ~ 0,
      TRUE ~ NA_real_
    )),
    depth = as.factor(if_else(level == 0, "foundational", "applied")),
    af = if_else(source == "af", 1, 0),
    ct = if_else(source %in% c("ct_ai", "ct_noai"), 1, 0),
    ct_ai = if_else(source == "ct_ai", 1, 0),
    ct_noai = if_else(source == "ct_noai", 1, 0),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci), # nolint
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    logit_cit_norm_perc = log(
      citation_normalized_percentile_value /
        (1 - citation_normalized_percentile_value)
    ),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field), # nolint
    publication_date = as.Date(publication_date),
    ln1p_resolution = log1p(as.numeric(resolution_mean)),
    ln1p_R_free = log1p(as.numeric(R_free_mean)),
    ln1p_score = log1p(as.numeric(score_mean)),
    year = as.integer(str_sub(publication_date, 1, 4)),
    quarter_year = paste0(
      year(publication_date), " Q", quarter(publication_date)
    ),
  ) %>%
  rename(author = last_author)


pubs_per_quarter <- papers %>%
  group_by(author, quarter_year) %>%
  summarise(num_publications = n()) %>%
  ungroup()


# Merge the summary back into the original papers DataFrame
papers <- papers %>%
  left_join(pubs_per_quarter, by = c("author", "quarter_year"))

# Fill NA values in 'field_' and 'mesh_' prefix columns with 0
papers <- papers %>%
  mutate(across(starts_with("field_"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("mesh_"), ~ replace_na(., 0)))

# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

papers$primary_field <- recode(papers$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# SUBSET
# ------------------------------------------------------------------------------

# Define sub_samples as a list of samples
sub_samples <- list()
sub_groups <- c("All PDB", "High PDB")
unique_scopes <- c("All", "Intent") # Changed to All/Intent distinction
unique_fields <- c(
  "All Fields",
  "Molecular Biology",
  "Medicine"
)

papers <- papers %>%
  select(c(
    "quarter_year",
    "author",
    "institution_type",
    "institution_country_code",
    "institution_cited_by_count",
    "institution_2yr_mean_citedness",
    "institution_h_index",
    "institution_i10_index",
    "num_publications",
    "ln1p_cited_by_count",
    "ln1p_fwci",
    "logit_cit_norm_perc",
    "patent_count",
    "patent_citation",
    "ca_count",
    "ln1p_resolution",
    "ln1p_R_free",
    "pdb_submission",
    "af",
    "ct_ai",
    "ct_noai",
    "strong",
    "depth",
    "primary_field",
    "q4_pdb_pre2021_any",
    "mesh_C",
    grep("^field_", names(papers), value = TRUE),
    "num_uniprot_structures",
    "num_pdb_ids",
    "num_primary_submissions",
    "num_diseases",
    "organism_rarity_mean",
    "mean_tmscore",
    "ln1p_score",
    "num_pdb_ids",
    "num_uniprot_structures_w_disease",
    "num_primary_submissions_w_disease",
    "num_uniprot_structures_w_rare_organisms",
    "num_primary_submissions_w_rare_organisms",
    "num_uniprot_structures_w_low_similarity",
    "num_primary_submissions_w_low_similarity"
  ))

# Create subsets for all combinations of depth, field, and sub_group
for (scope_lvl in unique_scopes) {
  for (field in unique_fields) {
    for (sub_group in sub_groups) {
      sample_name <- paste0(
        "scope_", scope_lvl, "__field_",
        field, "__subgroup_", sub_group
      )
      message("Creating sample: ", sample_name)

      # Start with the full dataset
      sub_sample <- papers

      # # Apply depth filter
      if (scope_lvl == "Intent") {
        sub_sample <- subset(sub_sample, !is.na(strong))
      } else {
        # drop strong
        sub_sample <- subset(sub_sample, select = -strong)
      }

      # Apply field filter
      if (field != "All Fields") {
        sub_sample <- subset(sub_sample, primary_field == field)
      }

      # Apply sub_group filter
      if (grepl("High PDB", sub_group)) {
        sub_sample <- subset(sub_sample, q4_pdb_pre2021_any == TRUE)
      }

      # Store the subset
      sub_samples[[sample_name]] <- sub_sample
    }
  }
}
# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
