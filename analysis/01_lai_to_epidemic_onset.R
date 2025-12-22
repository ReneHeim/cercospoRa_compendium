###############################################################################
# Project:     cercospoRa
# Script:      01_lai_to_epidemic_onset.R
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

# Change to your root directory
home <- here::here()

# Load and format weather data
wethr <- read.csv(system.file("extdata", "clean_weather.csv",
                              package = "cercospoRa"))
wethr <- format_weather(wethr, time_zone = "UTC")

generate_EO_maps <- function(platform,
                             start_date = "2022-04-25",
                             end_date   = "2022-09-30",
                             sowing_date = "2022-03-31",
                             cultivar_sus_values = 1:9) {

  print(platform)   # temporary debug
  target_res <- switch(
    platform,
    "S2"                = 10,
    "S2_superresolution" = 10,
    "UAV"               = 10,
    10
  )

  print(target_res) # temporary debug

  # Directory containing LAI maps for this platform
  img_dir <- here("output", "LAI_maps", platform)

  # Read all available LAI maps and fit pixel-wise growth parameters
  epidemic_onset_param <- read_sb_growth_parameter(
    img_files = list.files(img_dir, pattern = "tif$", full.names = TRUE),
    img_dates = as.POSIXct(
      as.Date(
        list.files(img_dir, pattern = "tif$", full.names = FALSE),
        format = "%Y_%m_%d"
      ), tz = "UTC"),
    target_res = target_res
  )

  # Estimate sigmoid parameters (growth rate r and x0) per pixel
  param_rxt <- calc_r_x0(
    epidemic_onset_param,
    min_r = 0.02,
    max_r = 0.05,
    k = 6
  )

  # Determine canopy closure dates (LAI crosses threshold)
  c_closure <- calc_c_closure(
    param_rxt,
    x1 = 1.3,
    k = 6
  )

  # Helper to compute and save epidemic onset map for one cultivar susceptibility
  calc_epidemic_onset_for_sus <- function(cultivar_sus) {

    epidemic_onset_map <- calc_epidemic_onset_from_image(
      start   = as.POSIXct(start_date, tz = "UTC"),
      end     = as.POSIXct(end_date,   tz = "UTC"),
      cc_r    = c_closure,
      weather = wethr,
      cultivar_sus = cultivar_sus
    )

    sowing_num <- as.numeric(as.POSIXct(sowing_date, tz = "UTC"))
    epidemic_onset_map <- (epidemic_onset_map - sowing_num) / (3600 * 24)

    out_dir <- here("output", "Predicted_epidemic_onset", platform)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    file_path <- file.path(
      out_dir,
      paste0("Ep_onset_", cultivar_sus, ".tif")
    )

    writeRaster(
      epidemic_onset_map,
      file_path,
      filetype = "GTiff",
      gdal = "COMPRESS=LZW",
      overwrite = TRUE
    )

    terra::plot(epidemic_onset_map, main = paste(platform, "sus =", cultivar_sus),
                xlab = "Easting (m)",
                ylab = "Northing (m)")
  }

  for (cs in cultivar_sus_values) {
    calc_epidemic_onset_for_sus(cs)
  }
}

platforms <- c("S2", "S2_superresolution", "UAV")

for (platform in platforms) {
  generate_EO_maps(platform)
}

library(usethis)
suppressWarnings(
  sessioninfo::session_info(
    to_file = here("output", "Predicted_epidemic_onset", 
                   "01_lai_to_epidemic_onset_session.log")
    )
)
