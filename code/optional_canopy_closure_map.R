###############################################################################
# Project:     cercospoRa
# Script:      optional_canopy_closure_map.R
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

# # ---- Dependency management --------------------------------------------------
# required_packages <- c("cercospoRa", "terra", "here", "sessioninfo", "usethis")
# 
# # Install any missing packages
# installed <- required_packages %in% rownames(installed.packages())
# if (any(!installed)) {
#   install.packages(required_packages[!installed])
# }
# 
# lapply(required_packages, library, character.only = TRUE)

library(cercospoRa)
library(terra)
library(here)

base_dir    <- here("output", "lai_maps")
destination <- here("output", "cc_maps")

# Ensure destination folder exists
if (!dir.exists(destination)) {
  dir.create(destination, recursive = TRUE)
}

generate_CC <- function(base_dir, platform) {
  platform_dir <- file.path(base_dir, platform)
  files_full   <- list.files(platform_dir, pattern = "tif", full.names = TRUE)
  files_dates  <- list.files(platform_dir, pattern = "tif", full.names = FALSE)
  
  epidemic_onset_param <- read_sb_growth_parameter(
    files_full,
    img_dates = as.POSIXct(as.Date(files_dates, format = "%Y%m%d"), tz = "UTC"),
    10
  )
  
  param_rxt <- calc_r_x0(epidemic_onset_param,
                         min_r = 0.02,
                         max_r = 0.05,
                         k = 6)
  
  c_closure <- calc_c_closure(param_rxt, x1 = 1.3, k = 6)
  
  # Turn into days after sowing (where "2022-03-31" is the sowing date)
  c_closure <- c_closure - as.numeric(as.Date("2022-03-31"))
  names(c_closure) <- "canopy_closure"
  
  # Output path
  tif_path <- file.path(destination, paste0(platform, "_cc.tif"))
  writeRaster(c_closure, tif_path, overwrite = TRUE)
}

platforms <- c("s2", "s2_superresolution", "uas")

for (platform in platforms) {
  generate_CC(base_dir, platform)
}

library(usethis)
suppressWarnings(
  sessioninfo::session_info(
    to_file = here("output", "cc_maps", 
                   "cc_maps.log")
  )
)
