# Clean out the workspace
rm(list = ls())
options(max.print = 1000)

# Check installation & load required packages
list_of_packages <- c(
  "arrow", "tidyverse", "MatchIt", "fastDummies", "aws.s3",
  "yaml", "zoo"
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
pathdir <- "data/05_model_output/papers/"

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
path <- "oct/04_output/publications/outputs.parquet" # nolint


# Fetch the data from the S3 bucket
papers <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# Replace commas in column names
colnames(papers) <- gsub(",", "", colnames(papers))

# for duplicate ids, keep first
papers <- papers %>% distinct(id, .keep_all = TRUE)

# Create 'strong' variable
papers <- papers %>%
  mutate(
    strong = if_else(chain_label %in% c("strong", "partial_strong"), 1, 0),
    strong = replace_na(strong, 0),
    strong = as.factor(strong),
    depth = as.factor(if_else(level == 0, "foundational", "applied")),
    treatment_af_dyn = if_else(source == "af", 1, 0),
    treatment_ct_dyn = if_else(source %in% c("ct_ai", "ct_noai"), 1, 0),
    treatment_ct_ai_dyn = if_else(source == "ct_ai", 1, 0),
    treatment_ct_noai_dyn = if_else(source == "ct_noai", 1, 0),
    group_pdb_count = group_pdb_count / max(group_pdb_count, na.rm = TRUE),
    institution = as.factor(last_author_institution),
    institution_type = as.factor(type),
    institution_country_code = as.factor(country_code),
    # ln1p_cit_0 = log1p(cit_0), # nolint
    # ln1p_cit_1 = log1p(cit_1), # nolint
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_norm_perc = log1p(citation_normalized_percentile_value),
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    publication_date = as.Date(publication_date),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free),
    time_qtly = (year(publication_date) - 2018) * 4 + quarter(publication_date),
    time_qtly = as.factor(time_qtly)
  )

pubs_per_quarter <- papers %>%
  group_by(last_author, time_qtly) %>%
  summarize(num_publications = n()) %>%
  ungroup()

# Merge the summary back into the original papers DataFrame
papers <- papers %>%
  left_join(pubs_per_quarter, by = c("last_author", "time_qtly"))

# Fill NA values in 'field_' and 'mesh_' prefix columns with 0
papers <- papers %>%
  mutate(across(starts_with("field_"), ~ replace_na(., 0))) %>%
  mutate(across(starts_with("mesh_"), ~ replace_na(., 0)))

# Create 'high_pdb' variable
papers <- papers %>%
  mutate(high_pdb = if_else(group_pdb_count >= quantile(group_pdb_count, 0.75, na.rm = TRUE), 1, 0)) # nolint

# improve field names in primary_field
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "biochem_genetics_molecular_biology", # nolint
  "Medicine" = "medicine",
  "Chemistry" = "chemistry",
  "Immunology and Microbiology" = "immunology_microbiology"
)

# Apply the mapping to the primary_field column
papers$primary_field <- recode(papers$primary_field, !!!field_mapping)

# ------------------------------------------------------------------------------
# SUBSET
# ------------------------------------------------------------------------------

sub_samples <- list(
  depth_all__field_all = papers,
  depth_foundational__field_all = subset(papers, depth == "foundational"),
  depth_applied__field_all = subset(papers, depth == "applied")
)

tech_groups <- c(
  "all", "ct", "ct_ai", "ct_noai", "w_high_pdb",
  "ct_w_high_pdb", "ct_ai_w_high_pdb", "ct_noai_w_high_pdb"
)
unique_depths <- c("all", "foundational", "applied")
unique_fields <- c(
  "biochem_genetics_molecular_biology", "medicine", 
  "chemistry", "immunology_microbiology"
)

percentile_75_pdb <- quantile(papers$group_pdb_count, 0.75)

# function to create subsets based on tech group
create_tech_group_subset <- function(data, tech_group) {
  if (tech_group == "all") {
    return(data)
  } else if (tech_group == "ct") {
    return(
      data %>% # nolint
        filter(
          stringr::str_detect(source, "af") |
            stringr::str_detect(source, "ct_ai") |
            stringr::str_detect(source, "ct_noai")
        )
    )
  } else if (tech_group == "ct_ai") {
    return(
      data %>% filter( # nolint
        stringr::str_detect(source, "af") |
          stringr::str_detect(source, "ct_ai")
      )
    )
  } else if (tech_group == "ct_noai") {
    return(
      data %>% filter( # nolint
        stringr::str_detect(source, "af") |
          stringr::str_detect(source, "ct_noai")
      )
    )
  } else if (tech_group == "w_high_pdb") {
    return(
      data %>% filter(group_pdb_count >= percentile_75_pdb) # nolint
    )
  } else if (tech_group == "ct_w_high_pdb") {
    return(
      data %>% # nolint
        filter(group_pdb_count >= percentile_75_pdb) %>% # nolint
        filter(
          stringr::str_detect(source, "af") |
            stringr::str_detect(source, "ct_ai") |
            stringr::str_detect(source, "ct_noai")
        )
    )
  } else if (tech_group == "ct_ai_w_high_pdb") {
    return(
      data %>% # nolint
        filter(group_pdb_count >= percentile_75_pdb) %>% # nolint
        filter(
          stringr::str_detect(source, "af") |
            stringr::str_detect(source, "ct_ai")
        )
    )
  } else if (tech_group == "ct_noai_w_high_pdb") {
    return(
      data %>% # nolint
        filter(group_pdb_count >= percentile_75_pdb) %>% # nolint
        filter(
          stringr::str_detect(source, "af") |
            stringr::str_detect(source, "ct_noai")
        )
    )
  }
}

for (depth_lvl in unique_depths) {
  for (field in unique_fields) {
    for (tech_group in tech_groups) {
      sample_name <- paste0(
        "depth_", depth_lvl, "__field_", field, "__tech_", tech_group
      )
      if (depth_lvl == "all") {
        sub_samples[[sample_name]] <- create_tech_group_subset(
          subset(papers, primary_field == field), tech_group
        )
      } else {
        sub_samples[[sample_name]] <- create_tech_group_subset(
          subset(papers, depth == depth_lvl & primary_field == field),
          tech_group
        )
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Save data
# ------------------------------------------------------------------------------
saveRDS(sub_samples, paste0(pathdir, "sub_samples.rds"))
