# =============================================================================
# EVENT STUDY COMBINED PLOTS SCRIPT
# =============================================================================
#
# This script combines individual event study plots into publication-ready
# multi-panel figures for the AlphaFold impact analysis. It produces:
#
# 1. Combined event study plots arranging individual PNG plots into 2x4 grids
# 2. Multi-panel figures for different research groups: Laboratories and
#    Established Researchers
# 3. Selected outcome variables: Field-Weighted Citation Impact, Publications,
#    TM-Score, Fold Identity, PDB Submissions, Disease-relevant Structures,
#    Clinical Citations, and Patent Citations
# 4. Publication-ready formatting with proper titles, alignment, and spacing
# 5. High-resolution output (300 DPI) suitable for academic publications
# 6. Handles missing plots gracefully with placeholder text
#
# Loads existing event study plots from did_plots/af/ directories and combines
# them using cowplot and magick packages for image manipulation
# Output: Combined PNG plots saved to data/05_model_output/combined_did_plots/
# =============================================================================

# Clean workspace and set options
rm(list = ls())
options(max.print = 1000, width = 250)

# Load required packages
packages <- c(
  "tidyverse", "did", "ggplot2",
  "scales", "showtext", "cowplot", "magick"
)
invisible(lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

# Set up fonts
font_add_google("Mulish", "muli")
showtext_auto()

# Set working directory and paths
setwd("~/projects/alphafold-impact/")

labs_data_path <- "data/05_model_output/labs/quarterly/did_plots/af"
researchers_data_path <- "data/05_model_output/authors/nonecr/quarterly/did_plots/af" # nolint

# Create output directory
output_dir <- "data/05_model_output/combined_did_plots/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

selected_dep_vars <- c(
  "ln1p_fwci", "num_publications",
  "ln1p_max_tmscore", "ln1p_max_fident",
  "num_pdb_submissions", "num_diseases",
  "ca_count", "patent_count"
)

# Variable labels for plot titles
var_labels <- list(
  "ln1p_fwci" = "Field-Weighted Citation Impact (ln)",
  "num_publications" = "Number of Publications",
  "ln1p_max_tmscore" = "Max TM-Score (ln)",
  "ln1p_max_fident" = "Max Fold Identity (ln)",
  "num_pdb_submissions" = "PDB Submissions",
  "num_diseases" = "Disease-relevant Structures",
  "ca_count" = "Clinical Article Citations",
  "patent_count" = "Patent-Paper Citations"
)

# Function to load and prepare PNG plots
load_png_plot <- function(file_path, title = "") {
  if (file.exists(file_path)) {
    img <- magick::image_read(file_path)
    p <- ggdraw() + draw_image(img) # nolint
    if (title != "") {
      p <- p + labs(title = title) # nolint
    }
    return(p)
  } else {
    # Create empty plot if file doesn't exist
    return(ggplot() + # nolint
      geom_text(aes(x = 0.5, y = 0.5, label = paste("Missing:", basename(file_path))), # nolint
        size = 4, family = "muli"
      ) +
      xlim(0, 1) + # nolint
      ylim(0, 1) + # nolint
      theme_void() + # nolint
      labs(title = title) # nolint
    ) # nolint
  }
}

# # Create combined plots for Labs
# cat("Creating combined plot for Labs...\n")
# labs_plots <- list()

# for (i in seq_along(selected_dep_vars)) {
#   var <- selected_dep_vars[i]
#   file_path <- file.path(labs_data_path, paste0("event_study_", var, ".png"))
#   title <- var_labels[[var]]
#   labs_plots[[i]] <- load_png_plot(file_path, title)
# }

# # Arrange labs plots in 2x4 grid
# labs_combined <- plot_grid(
#   plotlist = labs_plots,
#   ncol = 2, nrow = 4,
#   align = "hv",
#   axis = "tblr"
# )

# # Add main title for labs
# labs_final <- plot_grid(
#   ggdraw() +
#     draw_label("Laboratories",
#       fontface = "bold", size = 64
#     ),
#   labs_combined,
#   ncol = 1,
#   rel_heights = c(0.05, 0.95)
# )

# # Save labs combined plot
# ggsave(
#   filename = file.path(output_dir, "event_study_combined_labs.png"),
#   plot = labs_final,
#   width = 16, height = 20, dpi = 300
# )

# Create combined plots for Researchers
cat("Creating combined plot for Researchers...\n")
researchers_plots <- list()

for (i in seq_along(selected_dep_vars)) {
  var <- selected_dep_vars[i]
  file_path <- file.path(
    researchers_data_path, paste0("event_study_", var, ".png")
  )
  title <- var_labels[[var]]
  researchers_plots[[i]] <- load_png_plot(file_path, title)
}

# Arrange researchers plots in 2x4 grid
researchers_combined <- plot_grid(
  plotlist = researchers_plots,
  ncol = 2, nrow = 4,
  align = "hv",
  axis = "tblr"
)

# Add main title for researchers
researchers_final <- plot_grid(
  ggdraw() +
    draw_label("Established Researchers",
      fontface = "bold", size = 64
    ),
  researchers_combined,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)

# Save researchers combined plot
ggsave(
  filename = file.path(output_dir, "event_study_combined_researchers.png"),
  plot = researchers_final,
  width = 16, height = 20, dpi = 300
)

cat("\n", rep("=", 60), "\n")
cat("EVENT STUDY COMBINED PLOTS COMPLETE\n")
cat(rep("=", 60), "\n")
cat("Combined plots saved to:", output_dir, "\n")
cat("- Labs: event_study_combined_labs.png\n")
cat("- Researchers: event_study_combined_researchers.png\n")
