# ------------------------------------------------------------------------------
# HELPER FUNCTION TO CHECK FOR PERFECT COLLINEARITY
# ------------------------------------------------------------------------------
# Function to check and fix variables that might cause perfect collinearity
fix_perfect_collinearity <- function(data, fe_vars, dep_var) {
  result_data <- data

  # Get all factor variables that might cause collinearity
  factor_vars <- names(data)[sapply(data, is.factor)]

  # For each fixed effect variable
  for (fe_var in fe_vars) {
    # Check each factor variable
    for (var in factor_vars) {
      # Skip if the variable is one of the fixed effects
      if (var %in% fe_vars) next

      # Group by the fixed effect and check if any factor level appears often
      # Only consider observations where dep_var is not all 0 within the group
      var_counts <- data %>% # nolint
        group_by(!!sym(fe_var)) %>% # nolint
        filter(sum(!!sym(dep_var)) > 0) %>% #
        group_by(!!sym(fe_var), !!sym(var)) %>% # nolint
        summarise(count = n(), .groups = "drop") %>% # nolint
        group_by(!!sym(var)) %>%
        summarise(fe_count = n(), .groups = "drop")

      # If any level appears in only a handful groups, recode it
      if (any(var_counts$fe_count < 1)) {
        message(
          "Fixing potential collinearity in variable: ", var,
          " for fixed effect: ", fe_var
        )

        # Create a new column with recoded values
        problematic_levels <- var_counts %>% # nolint
          filter(fe_count < 1) %>% # nolint
          pull(!!sym(var)) # nolint

        # Recode problematic levels to a common "other" category
        if (length(problematic_levels) > 0) {
          result_data[[var]] <- as.character(result_data[[var]])
          result_data[[var]][result_data[[var]] %in% problematic_levels] <- "other_combined" # nolint
          result_data[[var]] <- as.factor(result_data[[var]])
        }
      }
    }
  }

  return(result_data)
}

# ------------------------------------------------------------------------------