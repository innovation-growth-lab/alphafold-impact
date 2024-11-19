# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
options(width = 100)

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
pathdir <- "data/05_model_output/authors/nonecr/quarterly/data/"

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
  "AWS_ACCESS_KEY_ID" = credentials$s3_credentials$key, # nolint
  "AWS_SECRET_ACCESS_KEY" = credentials$s3_credentials$secret, # nolint
  "AWS_DEFAULT_REGION" = "eu-west-2" # nolint
)

# Define the S3 bucket and path
bucket <- "igl-alphafold"
path <- "oct/03_primary/nonecr/publications_quarterly.parquet" # nolint

# Fetch the data from the S3 bucket
nonecr_data <- s3read_using(
  FUN = arrow::read_parquet,
  object = path,
  bucket = bucket
)

# ------------------------------------------------------------------------------
# Data Prep
# ------------------------------------------------------------------------------

# create a new factor variable
nonecr_data <- nonecr_data %>%
  mutate(
    strong_af = ifelse(is.na(strong_af), 0, strong_af),
    strong_ct_ai = ifelse(is.na(strong_ct_ai), 0, strong_ct_ai),
    strong_ct_noai = ifelse(is.na(strong_ct_noai), 0, strong_ct_noai),
    ct = ct_ai + ct_noai,
    af_ind = ifelse(af > 0, 1, 0),
    ct_ind = ifelse(ct > 0, 1, 0),
    ct_ai_ind = ifelse(ct_ai > 0, 1, 0),
    ct_noai_ind = ifelse(ct_noai > 0, 1, 0),
    "af:ct_ai_ind" = ifelse(af > 0 & ct_ai > 0, 1, 0),
    "af:ct_noai_ind" = ifelse(af > 0 & ct_noai > 0, 1, 0),
    strong_af_ind = ifelse(strong_af > 0, 1, 0),
    strong_ct_ai_ind = ifelse(strong_ct_ai > 0, 1, 0),
    strong_ct_noai_ind = ifelse(strong_ct_noai > 0, 1, 0),
    "strong_af:strong_ct_ai_ind" = ifelse(strong_af > 0 & strong_ct_ai > 0, 1, 0), # nolint
    "strong_af:strong_ct_noai_ind" = ifelse(strong_af > 0 & strong_ct_noai > 0, 1, 0) # nolint
  )

# creating quantiles for institution controls
nonecr_data <- nonecr_data %>%
  mutate(
    institution_2yr_mean_citedness = factor(
      ntile(institution_2yr_mean_citedness, 4)
    ),
    institution_h_index = factor(ntile(institution_h_index, 4)),
    institution_i10_index = factor(ntile(institution_i10_index, 4)),
    institution_cited_by_count = factor(
      ntile(institution_cited_by_count, 4)
    )
  )

# fill with nan
nonecr_data <- nonecr_data %>%
  mutate(
    institution_type = ifelse(
      is.na(institution_type), "unknown", institution_type
    ),
    institution_country_code = ifelse(
      is.na(institution_country_code), "unknown", institution_country_code
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
    institution = ifelse(is.na(institution_works_count), 0, institution_works_count), # nolint
    ca_count = ifelse(is.na(ca_count), 0, ca_count),
    high_pdb = as.factor(ifelse(is.na(high_pdb), 0, high_pdb)),
    covid_share_2020 = ifelse(is.na(covid_share_2020), 0, covid_share_2020),
    across(starts_with("field_"), ~ ifelse(is.na(.), 0, .)),
    across(starts_with("mesh_"), ~ ifelse(is.na(.), 0, .))
  )

# create factors, log transforms
nonecr_data <- nonecr_data %>%
  mutate(
    author = as.factor(author),
    author_position = as.factor(author_position),
    depth = as.factor(depth),
    institution = as.factor(institution),
    institution_type = as.factor(institution_type),
    institution_country_code = as.factor(institution_country_code),
    ln1p_cited_by_count = log1p(cited_by_count),
    ln1p_fwci = log1p(fwci),
    ln1p_cit_0 = log1p(cit_0),
    ln1p_cit_1 = log1p(cit_1),
    ln1p_cit_norm_perc = log1p(percentile_value),
    logit_cit_norm_perc = log(
      percentile_value /
        (1 - percentile_value)
    ),
    num_pdb_submissions = num_publications_pdb,
    ln1p_ca_count = log1p(ca_count),
    ln1p_patent_count = log1p(patent_count),
    ln1p_patent_citation = log1p(patent_citation),
    primary_field = as.factor(primary_field),
    resolution = as.numeric(resolution),
    R_free = as.numeric(R_free)
  )
# Define the mapping of old values to new values
field_mapping <- c(
  "Biochemistry, Genetics and Molecular Biology" = "Molecular Biology"
)

# Apply the mapping to the primary_field column
nonecr_data$primary_field <- recode(nonecr_data$primary_field, !!!field_mapping)

# filter out non complete cases of columns starting with institution_
nonecr_data <- nonecr_data %>%
  filter(complete.cases(across(starts_with("institution_"))))

# ------------------------------------------------------------------------------
# CEM (Coarsened Exact Matching)
# ------------------------------------------------------------------------------

# Define the columns to be used for matching
coarse_cols <- c(
  "ln1p_cited_by_count", "num_publications", "covid_share_2020"
)

exact_cols <- c("institution_type", "institution_country_code", "high_pdb")

mode_function <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Filter and prepare data for collapsing
quarterly_cem <- nonecr_data %>%
  group_by(author) %>%
  mutate(af_ind = max(af_ind)) %>%
  ungroup() %>%
  filter(complete.cases(across(coarse_cols))) %>%
  filter(quarter %in% 192:196) %>%
  select(af_ind, author, all_of(coarse_cols), all_of(exact_cols)) %>%
  group_by(author) %>%
  summarise(
    af_ind = max(af_ind),
    across(coarse_cols, \(x) mean(x, na.rm = TRUE)),
    across(exact_cols, mode_function)
  )

# ------------------------------------------------------------------------------
# Sample Prep
# ------------------------------------------------------------------------------
# Define sub_samples as a list of samples

nonecr_data <- nonecr_data %>%
  select(
    "quarter",
    "author",
    "institution",
    "institution_type",
    "institution_country_code",
    "institution_cited_by_count",
    "institution_2yr_mean_citedness",
    "institution_h_index",
    "institution_i10_index",
    "num_publications",
    "ln1p_cited_by_count",
    "ln1p_cit_0",
    "ln1p_cit_1",
    "ln1p_fwci",
    "logit_cit_norm_perc",
    "ln1p_patent_count",
    "ln1p_patent_citation",
    "ln1p_ca_count",
    "resolution",
    "R_free",
    "num_pdb_submissions",
    "af_ind",
    "ct_ai_ind",
    "ct_noai_ind",
    "af",
    "ct_ai",
    "ct_noai",
    "strong_af_ind",
    "strong_ct_ai_ind",
    "strong_ct_noai_ind",
    "strong_af",
    "strong_ct_ai",
    "strong_ct_noai",
    "depth",
    "primary_field",
    "high_pdb",
    "covid_share_2020",
    grep("^field_", names(nonecr_data), value = TRUE),
    grep("^mesh_", names(nonecr_data), value = TRUE)
  )

colnames(nonecr_data) <- gsub(",", "", colnames(nonecr_data))


# Define sub_samples as a list of samples
sub_samples <- list()
sub_groups <- c(
  # "All PDB",
  "All PDB - CEM",
  # "High PDB"
  "High PDB - CEM"
)
unique_depths <- c("All Groups", "Foundational", "Applied")
unique_fields <- c(
  "All Fields",
  "Molecular Biology",
  "Medicine"
)

# Create subsets for all combinations of depth, field, and sub_group
for (depth_lvl in unique_depths) { # nolint
  for (field in unique_fields) {
    for (sub_group in sub_groups) {
      sample_name <- paste0(
        "depth_", depth_lvl, "__field_",
        field, "__subgroup_", sub_group
      )
      message("Creating sample: ", sample_name)

      # Start with the full dataset
      sub_sample <- nonecr_data

      # Apply depth filter
      if (depth_lvl == "Foundational") {
        sub_sample <- subset(sub_sample, depth == "foundational")
      } else if (depth_lvl == "Applied") {
        sub_sample <- subset(sub_sample, depth == "applied")
      }

      # Apply field filter
      if (field != "All Fields") {
        sub_sample <- subset(sub_sample, primary_field == field)
      }

      # Apply sub_group filter
      if (grepl("High PDB", sub_group)) {
        sub_sample <- subset(sub_sample, high_pdb == 1)
      }

      if (grepl("CEM", sub_group)) {
        # CEM matching on the collapsed data using coarse_cols
        match_out_af_coarse <- matchit(
          as.formula(paste0(
            "af_ind ~ ",
            paste(coarse_cols, collapse = " + ")
          )),
          data = quarterly_cem, method = "cem", k2k = FALSE
        )

        # Store matched data for coarse_cols
        cem_data_coarse <- match.data(match_out_af_coarse)

        # Exact matching on the collapsed data using exact_cols
        match_out_af_exact <- matchit(
          as.formula(paste0("af_ind ~ ", paste0(exact_cols, collapse = " + "))),
          data = quarterly_cem, method = "exact"
        )

        # Store matched data for exact_cols
        cem_data_exact <- match.data(match_out_af_exact)

        # Combine the matched results
        combined_cem_data <- intersect(
          cem_data_coarse$author,
          cem_data_exact$author
        )

        # Sample based on the combined matched group
        sub_sample <- sub_sample %>% filter(author %in% combined_cem_data)
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
