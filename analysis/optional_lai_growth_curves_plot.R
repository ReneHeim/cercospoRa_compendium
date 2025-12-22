###############################################################################
# Project:     cercospoRa
# Script:      optional_lai_growth_curves_plot.R
# Author:      Nathan Okole
# Affiliation: Institute for Sugar Beet Research
# Date:        2025-11-22
#
# Description:
#   This script uses the LAI maps obtained from the radiative transfer model
#   (PROSAIL) and creates epidemic onset maps. 
#
# Reproducibility:
#   Session information and dependencies are captured at the end of the script
#   using library(usethis) and sessioninfo::session_info()
#
# Usage:
#   Run the entire script or step-by-step using RStudio sections.
#
# Notes:
#   - Data files are not included in the repository due to size constraints.
#
###############################################################################

# ---- Dependency management --------------------------------------------------
required_packages <- c("cercospoRa", "terra", "here", "sessioninfo", "usethis")

# Install any missing packages
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

lapply(required_packages, library, character.only = TRUE)


# Required packages
library(cercospoRa)
library(terra)
library(here)

base_dir    <- here("output", "LAI_maps")
destination <- here("output", "LAI_progression_curve")

# Create destination folder if not present
if (!dir.exists(destination)) {
  dir.create(destination, recursive = TRUE)
}

plot_r <- function(base_dir, platform) {
  # Read available LAI maps (use here-based paths)
  platform_dir <- file.path(base_dir, platform)
  files_full   <- list.files(platform_dir, pattern = "tif", full.names = TRUE)
  files_dates  <- list.files(platform_dir, pattern = "tif", full.names = FALSE)
  
  epidemic_onset_param <-
    read_sb_growth_parameter(
      files_full,
      img_dates = as.POSIXct(as.Date(files_dates, format = "%Y_%m_%d"), tz = "UTC"),
      10
    )
  
  param_rxt <- calc_r_x0(epidemic_onset_param, min_r = 0.02, max_r = 0.05, k = 6)
  r  <- param_rxt$r
  x0 <- param_rxt$x0
  t0 <- param_rxt$t0
  
  table <- data.frame(r = matrix(r), x = matrix(x0))
  table <- na.omit(table)
  
  # Output PNG path using here()
  png_path <- file.path(destination, paste0(platform, '_growth_curve.png'))
  png(png_path, units = "in", width = 6, height = 5, res = 600)
  
  plot_sigma <- function(t, r_i, x_i, ti_0) {
    t    <- as.numeric(t)
    ti_0 <- as.numeric(as.Date(ti_0))
    N <- 6 / (1 + ((6 - x_i) / x_i) * exp(-r_i * (t - ti_0)))
    return(N)
  }
  
  # Plot the first curve
  r_1 <- table$r[1]
  x_1 <- table$x[1]
  t_i <- seq(as.Date("2022-04-01"), as.Date("2022-10-15"), 1)
  N_i <- plot_sigma(t_i, r_1, x_1, t0)
  par(mgp = c(2.5, 0.5, 0))
  plot(t_i, N_i, type = 'l',
       xlab = 'Date',
       ylab = 'Leaf Area Index (m²/m²)',
       ylim = c(0, 6),
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
       cex.axis = 1.25,
       cex.lab  = 1.5,
       main     = paste("LAI Progression Curve (", platform, ")", sep = ""))
  
  # Add the remaining curves
  for (i in 2:nrow(table)) {
    r_i <- table$r[i]
    x_i <- table$x[i]
    N_i <- plot_sigma(t_i, r_i, x_i, t0)
    lines(t_i, N_i, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
  }
  
  dev.off()
}

platforms <- c("S2", "S2_superresolution", "UAV")

for (platform in platforms) {
  plot_r(base_dir, platform)
}

library(usethis)
suppressWarnings(
  sessioninfo::session_info(
    to_file = here("output", "LAI_progression_curve", 
                   "LAI_progression_curve.log")
  )
)
