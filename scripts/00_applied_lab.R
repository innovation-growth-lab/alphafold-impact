# Clean out the workspace
rm(list = ls())
options(max.print = 1000)
source("scripts/utils.R")

# Check installation & LOAD PACKAGES
list_of_packages <- c(
  "arrow", "tidyverse", "ggplot2", "data.table", "bacondecomp",
  "fixest", "broom", "stargazer", "kableExtra", "patchwork",
  "extrafont", "RColorBrewer", "plotrix", "MatchIt", "cem"
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

setwd("~/projects/alphafold-impact/")

figures <- "data/05_model_output/figures/"
tables <- "data/05_model_output/tables/"
figs <- list()

getwd()

select <- dplyr::select
summarise <- dplyr::summarise
summarize <- dplyr::summarize

################################################################################
################################################################################
############################### DATA PREP ######################################
################################################################################
################################################################################

sb_data_qtly <- read_parquet(
  "data/04_outputs/applied_labs/staggered/outputs_quarterly.parquet"
)

# print columns
colnames(sb_data_qtly)

# replace columns with comma by ""
colnames(sb_data_qtly) <- gsub(",", "", colnames(sb_data_qtly))

# print parent_time values
unique(sb_data_qtly$parent_time)

# count seed values
sb_data_qtly %>%
  group_by(seed) %>%
  count()


# make strong be 1 or 0 from TRUE FALSE
sb_data_qtly$strong <- as.factor(sb_data_qtly$strong)
sb_data_qtly$time_qtly <- as.factor(sb_data_qtly$time)
sb_data_qtly$country <- as.factor(sb_data_qtly$institution_country_code)
sb_data_qtly$institution_type <- as.factor(sb_data_qtly$institution_type)


# create missing rows (time should go 1 - 25) for each pi_id, fill with NA
sb_data_qtly <- sb_data_qtly %>%
  group_by(pi_id, parent_time) %>%
  complete(time = 1:25) %>%
  ungroup()

# fill seed, strong NA with bfill, other columns with 0
sb_data_qtly <- sb_data_qtly %>%
  group_by(pi_id) %>%
  fill(seed, strong, country, starts_with("institution_"), .direction = "downup") %>%
  replace_na(list(
    num_publications = 0,
    cited_by_count = 0,
    ct0 = 0,
    ct1 = 0,
    ca_count = 0,
    patent_count = 0,
    grant_count = 0,
    pdb_share = 0,
    resolution = 0,
    R_free = 0,
    covid_share_2020 = 0,
  )) %>%
  ungroup()

# Replace NA values in all columns that start with "mesh_" with 0
sb_data_qtly <- sb_data_qtly %>%
  mutate_at(vars(starts_with("mesh_")), ~replace_na(., 0))


sb_data_qtly <- sb_data_qtly %>%
  rename(g = parent_time) %>%
  mutate(
    g = ifelse(seed == "other", NA, g),
    g_af = ifelse(str_detect(seed, "af"), g, NA),
    g_ct = ifelse(str_detect(seed, "ct"), g, NA),
    treatment_af_dyn = ifelse(time > g_af, 1, 0),
    treatment_af = ifelse(time > 15 & !is.na(g_af), 1, 0),
    treatment_ct = ifelse(time > g_ct, 1, 0),
    experimental_share = experimental / num_publications,
    protein_share = protein_concept / num_publications
  ) %>%
  replace_na(
    list(
      treatment_af_dyn = 0,
      treatment_af = 0,
      treatment_ct = 0,
      pdb_share = 0
    )
  ) %>%
  select(c(g, g_af, g_ct, time), everything()) %>%
  arrange(pi_id, time)

# bfill pdb_share, resolution, R_free
sb_data_qtly <- sb_data_qtly %>%
  select(-institution_country_code) %>%
  group_by(pi_id) %>%
  fill(pdb_share, resolution, R_free) %>%
  ungroup()

# take log(1+y) for num_publications, ca_count, patent_count, grant_count
sb_data_qtly <- sb_data_qtly %>%
  mutate(
    num_publications = log1p(num_publications),
    ca_count = log1p(ca_count),
    patent_count = log1p(patent_count),
    grant_count = log1p(grant_count)
  )

# create category columns for whether seed is in af, ct, or other
sb_data_qtly$seed_factor <- as.factor(sb_data_qtly$seed)

# perform year-standardised patent_count
sb_data_qtly <- sb_data_qtly %>%
  group_by(time) %>%
  mutate(
    num_publications_std = scale(num_publications, center = TRUE, scale = TRUE),
    patent_count_std = scale(patent_count, center = TRUE, scale = TRUE),
    grant_count_std = scale(grant_count, center = TRUE, scale = TRUE),
    ca_count_std = scale(ca_count, center = TRUE, scale = TRUE),
    cited_by_count_std = scale(cited_by_count, center = TRUE, scale = TRUE)
  ) %>%
  ungroup()

################################################################################
################################################################################
################################# SUBSET #######################################
################################################################################
################################################################################

### preparing pdb subset
sb_data_qtly_2018_2020 <- sb_data_qtly %>%
  filter(time %in% 1:13)

# calculate the average pdb_share for each pi_id
avg_pdb_share <- sb_data_qtly_2018_2020 %>%
  group_by(pi_id) %>%
  summarise(avg_pdb_share = mean(pdb_share, na.rm = TRUE))

# calculate the 75th percentile of the averages
percentile_75_pdb  <- quantile(avg_pdb_share$avg_pdb_share, 0.75, na.rm = TRUE)

### preparing field_biochemistry_genetics_and_molecular_biology subset
sb_data_qtly_2018_2020 <- sb_data_qtly %>%
  filter(time %in% 1:13)

# calculate the average field_biochemistry_genetics_and_molecular_biology
avg_field_biochemistry_genetics_and_molecular_biology <- sb_data_qtly_2018_2020 %>% # nolint
  group_by(pi_id) %>%
  summarise(
    avg_field_biochemistry_genetics_and_molecular_biology = mean(
      field_biochemistry_genetics_and_molecular_biology,
      na.rm = TRUE
    )
  )

# calculate the 75th percentile of the averages
percentile_75_bio <- quantile(
  avg_field_biochemistry_genetics_and_molecular_biology$avg_field_biochemistry_genetics_and_molecular_biology, # nolint
  0.75,
  na.rm = TRUE
)

### preparing the protein_share subset
sb_data_qtly_2018_2020 <- sb_data_qtly %>%
  filter(time %in% 1:13)

# calculate the average protein_share for each pi_id
avg_protein_share <- sb_data_qtly_2018_2020 %>%
  group_by(pi_id) %>%
  summarise(avg_protein_share = mean(protein_share, na.rm = TRUE))

# calculate the 75th percentile of the averages
percentile_75_prot  <- quantile(avg_protein_share$avg_protein_share, 0.75, na.rm = TRUE) # nolint

### preparing the experimental_share subset
sb_data_qtly_2018_2020 <- sb_data_qtly %>%
  filter(time %in% 1:13)

# calculate the average experimental_share for each pi_id
avg_experimental_share <- sb_data_qtly_2018_2020 %>%
  group_by(pi_id) %>%
  summarise(avg_experimental_share = mean(experimental_share, na.rm = TRUE))

# calculate the 75th percentile of the averages
percentile_75_exp <- quantile(avg_experimental_share$avg_experimental_share, 0.75, na.rm = TRUE) # nolint
percentile_75_prot <- quantile(avg_protein_share$avg_protein_share, 0.75, na.rm = TRUE) # nolint


### SUBSETS ###
sub_samples <- list(
  "all" = sb_data_qtly,
  "af_other" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | seed == "other"),
  "af_ct" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | str_detect(seed, "ct")),
  "af_other_w_pdb" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | seed == "other") %>%
    filter(!is.na(pdb_share)),
  "af_ct_w_high_pdb" = sb_data_qtly %>%
    filter(str_detect(seed, "af") | seed == "ct") %>%
    semi_join(avg_pdb_share %>%
      filter(avg_pdb_share > percentile_75_pdb), by = "pi_id"), # nolint
  "pi_id_high_pdb" = sb_data_qtly %>%
    semi_join(avg_pdb_share %>%
      filter(avg_pdb_share > percentile_75_pdb), by = "pi_id"), # nolint
  "pi_id_high_field" = sb_data_qtly %>%
    semi_join(avg_field_biochemistry_genetics_and_molecular_biology %>%
      filter(
        avg_field_biochemistry_genetics_and_molecular_biology >
          percentile_75_bio
      ), by = "pi_id"),
  "pi_id_high_protein" = sb_data_qtly %>%
    semi_join(avg_protein_share %>%
      filter(avg_protein_share > percentile_75_prot), by = "pi_id"), # nolint
  "pi_id_high_experimental" = sb_data_qtly %>%
    semi_join(avg_experimental_share %>%
      filter(avg_experimental_share > percentile_75_exp), by = "pi_id") # nolint
)

# get na per column in sub_samples[["af_ct"]], export as excel
sub_samples[["af_ct"]] %>%
  summarise_all(~sum(is.na(.))) %>%
  write_csv("data/05_model_output/tables/na_af_ct.csv")

################################################################################
################################################################################
################################ CEM ###########################################
################################################################################
################################################################################

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
  "mesh_C",
  "country",
  "institution_type",
  "institution_works_count"
)

# drop rows with missing "country"
sb_data_qtly <- sb_data_qtly %>%
  filter(!is.na(country))

sb_data_qtly_t0_cem <- sb_data_qtly %>%
  filter(time == 0) %>%
  mutate(
    treatment_af_ct = ifelse(!is.na(g_af) & is.na(g_ct), 1, 0)
  ) %>%
  filter(complete.cases(
    field_biochemistry_genetics_and_molecular_biology, mesh_C
  ))


# Convert the cols vector to a string with + between each variable
cols_str <- paste0(cols, collapse = " + ")

# Perform CEM
match_out_af <- matchit(
  as.formula(paste0("treatment_af_ct ~ ", cols_str)),
  data = sb_data_qtly_t0_cem, method = "cem", k2k = TRUE
)

# groupby seed and count
match.data(match_out_af) %>%
  group_by(seed) %>%
  count()

cem_data <- match.data(match_out_af)

# get the sb_data_qtly data corresponding to the pi_id in cem_data
sb_data_qtly_cem <- sb_data_qtly %>%
  semi_join(cem_data, by = "pi_id")

# add to sublists
sub_samples[["af_ct_cem"]] <- sb_data_qtly_cem
################################################################################
################################################################################
############################## REGRESSIONS #####################################
################################################################################
################################################################################

field_cols <- grep("^field_", names(sb_data_qtly), value = TRUE)
mesh_cols <- grep("^mesh_", names(sb_data_qtly), value = TRUE)
institution_cols <- grep("^institution_", names(sb_data_qtly), value = TRUE)

covs <- list()
covs[["base0"]] <- c("covid_share_2020")
covs[["base1"]] <- c(covs$base0, "pdb_share", "protein_share", "experimental_share") # nolint
covs[["base2"]] <- c(covs$base1, field_cols)
covs[["base3"]] <- c(covs$base2, mesh_cols, institution_cols)
covs[["interact"]] <- c(
  "protein_concept:time_qtly",
  "experimental:time_qtly",
  "strong:time_qtly",
  paste0(field_cols, ":time_qtly"),
  paste0(mesh_cols, ":time_qtly")
)

fes <- list()
fes[["fe1"]] <- c("pi_id")
fes[["fe2"]] <- c(fes$fe1, "time_qtly", "country")

### Extensive ###

# List of dependent variables
dep_vars <- c(
  "num_publications", "tcc", "ct0", "ct1",
  "ca_count", "patent_count", "patent_citation", "cited_by_count_std"
)

# List of covariate sets
cov_sets <- c("base3")

# List of fixed effects
fe_list <- c("fe2")

# List of treatment variables
treat_vars <- c(
  "treatment_af_dyn",
  "treatment_af_dyn + treatment_af_dyn:strong:protein_share + treatment_af_dyn:strong + treatment_af_dyn:protein_share", # nolint
  "treatment_af_dyn + protein_share + experimental_share + strong + treatment_af_dyn:strong + treatment_af_dyn:protein_share + treatment_af_dyn:experimental_share" # nolint
)

form_list <- list()
# Iterate over dependent variables
for (dep_var in dep_vars) { # nolint
  # Iterate over covariate sets
  for (cov_set in cov_sets) {
    # Iterate over fixed effects
    for (fe in fe_list) {
      # Iterate over treatment variables
      for (treat_var in treat_vars) {
        # Check if covs[[cov_set]] is empty
        if (length(covs[[cov_set]]) == 0) {
          # Create formula without '+' before '|'
          form_list[[
            paste0(
              dep_var, "_", cov_set, "_", fe, "_", gsub(" ", "_", treat_var)
            )
          ]] <- as.formula(
            paste0(
              dep_var, " ~ ", treat_var, " |",
              paste0(fes[[fe]], collapse = " + ")
            )
          )
        } else {
          # Create formula with '+' before '|'
          form_list[[
            paste0(
              dep_var, "_", cov_set, "_", fe, "_", gsub(" ", "_", treat_var)
            )
          ]] <- as.formula(
            paste0(
              dep_var, " ~ ", treat_var, " +",
              paste0(covs[[cov_set]], collapse = " + "),
              "|", paste0(fes[[fe]], collapse = " + ")
            )
          )
        }
      }
    }
  }
}

results <- list()
# For each subset, compute feols
for (sub in names(sub_samples)) {
  # For each formula, compute feols
  for (form in names(form_list)) {
    results[[paste0(sub, "_", form)]] <- feols(
      form_list[[form]],
      data = sub_samples[[sub]],
      cluster = "pi_id"
    )
  }
}

variable_interest <- c(
  "treatment_af_dyn", "strong", "protein_share",
  "treatment_af_dyn:strong", "treatment_af_dyn:protein_share",
  "treatment_af_dyn:strong:protein_share"
)
variable_extended_interest <- c(
  "treatment_af_dyn", "strong", "protein_share",
  "experimental_share", "treatment_af_dyn:strong",
  "treatment_af_dyn:protein_share", "treatment_af_dyn:experimental_share"
)

### Clinical Article Citations
etable(results[c(
  "all_ca_count_base3_fe2_treatment_af_dyn",
  "all_ca_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_ca_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_ca_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_interest, file = paste0(tables, "01_applied/01_ca_translational.tex")) # nolint

### Time to Clinical Citation
etable(results[c(
  "all_tcc_base3_fe2_treatment_af_dyn",
  "all_tcc_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_tcc_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_interest, file = paste0(tables, "01_applied/02_tcc_translational.tex")) # nolint

### Patents
etable(results[c(
  "all_patent_count_base3_fe2_treatment_af_dyn",
  "all_patent_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_patent_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_patent_count_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_interest, file = paste0(tables, "01_applied/03_patent_count_translational.tex")) # nolint

etable(results[c(
  "all_patent_citation_base3_fe2_treatment_af_dyn",
  "all_patent_citation_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_patent_citation_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_patent_citation_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_interest, file = paste0(tables, "01_applied/04_patent_citation_translational.tex")) # nolint

### Productivity
etable(results[c(
  "all_ct0_base3_fe2_treatment_af_dyn",
  "all_ct0_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_ct0_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_ct0_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_extended_interest, file = paste0(tables, "01_applied/05_ct0_productivity.tex")) # nolint


etable(results[c(
  "all_ct1_base3_fe2_treatment_af_dyn",
  "all_ct1_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_ct1_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_ct1_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_extended_interest, file = paste0(tables, "01_applied/06_ct1_productivity.tex")) # nolint

etable(results[c(
  "all_num_publications_base3_fe2_treatment_af_dyn",
  "all_num_publications_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_num_publications_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_num_publications_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_extended_interest, file = paste0(tables, "01_applied/07_num_publications_productivity.tex")) # nolint


etable(results[c(
  "all_cited_by_count_std_base3_fe2_treatment_af_dyn",
  "all_cited_by_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_cited_by_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_cited_by_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_extended_interest, file = paste0(tables, "01_applied/08_cited_by_count_std_productivity.tex")) # nolint

### Grants
etable(results[c(
  "all_grant_count_std_base3_fe2_treatment_af_dyn",
  "all_grant_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_grant_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share", # nolint
  "af_ct_w_high_pdb_grant_count_std_base3_fe2_treatment_af_dyn_+_treatment_af_dyn:strong:pdb_share_+_treatment_af_dyn:strong_+_treatment_af_dyn:pdb_share" # nolint
)], keep = variable_interest, file = paste0(tables, "01_applied/09_grant_count_std_productivity.tex")) # nolint


###
basic_str <- as.formula(
  paste0(
    "ca_count ~ treatment_af_dyn + treatment_af_dyn:strong:pdb_share + treatment_af_dyn:strong + treatment_af_dyn:pdb_share +",
    paste0(covs$base3, collapse = " + "),
    "|",
    paste0(fes$fe2, collapse = " + ")
  )
)

reg <- feols(
  basic_str,
  data = sub_samples[["af_ct"]],
  cluster = "pi_id"
)

# print full row 82 of sub_files[["af_ct"]]
print(sub_samples[["af_ct"]][26693:26694, ], width=Inf)

summary(reg)
treat_vars <- c(
  "treatment_af_dyn",
  "treatment_af_dyn + treatment_af_dyn:strong:protein_share + treatment_af_dyn:strong + treatment_af_dyn:protein_share", # nolint
  "treatment_af_dyn + protein_share + experimental_share + strong + treatment_af_dyn:strong + treatment_af_dyn:protein_share + treatment_af_dyn:experimental_share" # nolint
)

rsb_data_qtly <- read_parquet(
  "data/04_outputs/applied_labs/staggered/outputs_quarterly.parquet"
)

# print rows in rsb_data_qtly where pi_id is A5005264477, show column time and strong
print(rsb_data_qtly %>% filter(pi_id == "A5005264477")  %>% select(time, strong), n = 30)

# get size of rsb_data_qtly
dim(rsb_data_qtly)
