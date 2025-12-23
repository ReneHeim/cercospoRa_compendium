###############################################################################
# Project:     cercospoRa
# Script:      00_create_lai_maps.R
# Author:      Nathan Okole
# Affiliation: Institute for Sugar Beet Research
# Date:        2025-11-22
#
# Description:
#   This script generates LAI maps using the PROSAIL radiative transfer model. 
#
# Reproducibility:
#   Session information and dependencies are captured at the end of the script
#   using library(usethis) and sessioninfo::session_info()
#
# Usage:
#   Run the entire script or step-by-step using RStudio sections.
#
# Notes:
#   - Data available: https://github.com/ReneHeim/cercospoRa_compendium/tree/main/data 
#
###############################################################################


# # ---- Dependency management --------------------------------------------------
# cran_packages <- c("terra", "here", "sessioninfo", "usethis", "dplyr", "sf")
# github_packages <- c("prosail")  # from jbferet/prosail
# 
# # 1) Install missing CRAN packages
# missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
# if (length(missing_cran) > 0) {
#   message("Installing missing CRAN packages: ",
#           paste(missing_cran, collapse = ", "))
#   install.packages(missing_cran)
# }
# 
# # 2) Install missing GitHub packages
# if (!requireNamespace("prosail", quietly = TRUE)) {
#   if (!requireNamespace("remotes", quietly = TRUE)) {
#     install.packages("remotes")
#   }
#   remotes::install_github("jbferet/prosail")
# }
# 
# # 3) Load all required packages (CRAN + GitHub)
# required_packages <- c(cran_packages, github_packages)
# invisible(lapply(required_packages, library, character.only = TRUE))

# ---- Custom function --------------------------------------------------------
tail <- function(x) {
  y <- x[length(x)]
  return(y)
}

# ---- Session info for reproducibility ---------------------------------------
sessioninfo::session_info()


print(paste("You are here:", here::here()))
library(sf)
library(terra)
library(prosail)
library(dplyr)

# Generate look up table (LUT)
# We are obtaining the parameters for sugar beet according to Jay, Sylvain, 
# Fabienne Maupas, Ryad Bendoula, and Nathalie Gorretta."Retrieving LAI, 
# chlorophyll and nitrogen contents in sugar beet crops from multi-angular optical 
# remote sensing: Comparison of vegetation indices and PROSAIL inversion for field 
# phenotyping." Field Crops Research 210 (2017): 33-46.


# ---- PROSAIL parameter grid -------------------------------------------------

solar_zenith  <- 10    # Solar zenith angle in degrees (sun elevation geometry)
deltaazimuth  <- 180   # Relative azimuth angle between sun and sensor (deg)
sensor_zenith <- 0     # Sensor zenith angle in degrees (nadir view)

LAI0 <- seq(0.1, 6, 0.01)      # Leaf Area Index values to explore in the LUT
N0   <- seq(1, 3.5, 0.5)       # Leaf structure parameter (N) range
Cab0 <- seq(20, 65, 7)         # Leaf chlorophyll content (µg/cm² or similar) range
Cw0  <- seq(0.03, 0.09, 0.015) # Leaf equivalent water thickness (cm or similar) range

# Generate all combinations of biophysical parameters (full factorial LUT)
InputPROSAIL <- expand.grid(
  N   = N0,    # Leaf structure
  CHL = Cab0,  # Chlorophyll content
  EWT = Cw0,   # Equivalent water thickness
  lai = LAI0   # Leaf Area Index
)

# Add canopy architecture and background parameters used by PROSAIL
InputPROSAIL$TypeLidf <- 1       # Leaf angle distribution type (1 = predefined distribution)
InputPROSAIL$LIDFa    <- -0.35   # LIDF parameter a (spherical distribution)
InputPROSAIL$LIDFb    <- -0.15   # LIDF parameter b (spherical distribution)
InputPROSAIL$q        <- 0.01    # Hot-spot parameter (controls brightness peak in backward direction)
InputPROSAIL$psoil    <- 0.3     # Fraction of dry soil background (0–1)
InputPROSAIL$tto      <- sensor_zenith  # Observation (sensor) zenith angle for each record
InputPROSAIL$tts      <- solar_zenith   # Solar zenith angle for each record
InputPROSAIL$psi      <- deltaazimuth   # Relative azimuth angle for each record

# Run PROSAIL to generate the spectral LUT for all parameter combinations
spect <- Generate_LUT_PROSAIL(
  SAILversion = "4SAIL",                       # Canopy model version
  InputPROSAIL = InputPROSAIL,                # Table of all input parameter sets
  SpecPROSPECT = prospect::SpecPROSPECT_FullRange, # Leaf optical properties (PROSPECT)
  SpecSOIL = SpecSOIL,                        # Soil reflectance spectrum
  SpecATM = SpecATM                           # Atmospheric transmittance/irradiance info
)

# Use only the BRF component from the PROSAIL output
spect <- spect$BRF

# ---- Sensor spectral definition ---------------------------------------------
# Spectral bands (center wavelength + FWHM) for the target sensor
# Here: Sentinel‑2 bands corresponding to MicaSense Altum channels
data_resampling_matrix <- data.frame(
  center = c(492, 560, 664, 704, 833),   # band center wavelengths (nm)
  fwhm   = c(66,  36,  31,  15, 106)    # full width at half maximum (nm)
)

# ---- LAI retrieval function --------------------------------------------------
# Retrieve LAI at ROI locations for a given reflectance image using a PROSAIL LUT.
retrieve_lai <- function(image_date,
                         folder_directory,
                         roi,
                         spect,                     # matrix: PROSAIL spectra (rows = LUT entries)
                         data_resampling_matrix,    # data frame: sensor band centers/FWHM
                         band_index = c(2:5, 8),    # bands to use from the image
                         platform = "s2",           # label for output folder (e.g. "s2" or "uas")
                         harmonized = FALSE,        # is image already harmonized (no offset)?
                         normalized = FALSE,        # is image already scaled to reflectance (0–1)?
                         aggregate = FALSE,         # aggregate (downscale) for UAV images?
                         agg_factor = 50,           # aggregation factor for spatial resolution
                         lai = InputPROSAIL$lai) {  # vector of LAI values for each LUT row
  
  # ---- Resample LUT to sensor bands -----------------------------------------
  # Build a synthetic sensor definition from the resampling matrix
  sensor <- get_spec_sensor(
    SensorName       = "sensor",               # arbitrary sensor name
    SpectralProps    = data_resampling_matrix, # center & FWHM for each band
    Path_SensorResponse = "./"                # path to store/read sensor response
  )
  
  # Retrieve the spectral response functions (SRF) for the synthetic sensor
  SRF_sensor <- GetRadiometry("sensor")
  
  # Convolve high‑resolution PROSPECT/PROSAIL spectra to sensor bands
  spect_resample <- applySensorCharacteristics(
    wvl   = prospect::SpecPROSPECT_FullRange$lambda, # original wavelength grid
    SRF   = SRF_sensor,                               # sensor spectral response
    InRefl = spect                                    # PROSAIL BRF LUT
  )
  
  # Name rows by band identifiers and convert to data.frame (LUT rows x bands)
  rownames(spect_resample) <- SRF_sensor$Spectral_Bands
  spect_resample <- as.data.frame(t(spect_resample))  # transpose so rows = LUT entries
  
  # ---- Read and preprocess image --------------------------------------------
  # Build full path to the image (expects <date>.tif in folder_directory)
  image_path <- file.path(folder_directory, paste0(image_date, ".tif"))
  
  # Read all bands as a SpatRaster
  allbands <- terra::rast(image_path)
  
  # Keep only the bands used for inversion (matching spect_resample)
  allbands <- allbands[[band_index]]  # e.g. c(2:5, 8) for S2
  
  # Crop and mask image to region of interest (reproject ROI if needed)
  cat(paste0("Convert to SpatRaster and crop: ", image_date, "\n"))
  allbands <- terra::crop(
    allbands,
    st_transform(roi, crs(allbands)),  # ensure ROI in same CRS as raster
    mask = TRUE
  )
  
  # Apply radiometric preprocessing:
  # - If not harmonized, remove additive offset (e.g. S2 L2A offset = 1000)
  if (!harmonized) {
    allbands <- allbands - 1000
  }
  
  # - If not normalized, scale to [0, 1] reflectance (e.g. divide by 10000)
  if (!normalized) {
    allbands <- allbands / 10000
  }
  
  # Optionally aggregate (downsample) for UAV data to reduce resolution/variance
  if (aggregate) {
    allbands <- terra::aggregate(allbands, fact = agg_factor, fun = "mean")
  }
  
  # ---- LAI inversion ---------------------------------------------------------
  # Copy LUT and LAI vector (kept for clarity; both passed to wrapper)
  spect_resample_ <- spect_resample   # PROSAIL LUT at sensor bands
  lai_vec         <- lai              # LAI associated with each LUT spectrum
  
  # Apply the inversion function pixel‑wise over the image bands
  # maximize_cos_wrapper() should return estimated LAI given pixel spectrum.
  lai_inversion <- terra::app(
    allbands,
    maximize_cos_wrapper(T_var = spect_resample_, lai = lai_vec),
    cores = 12                          # number of parallel cores
  )
  
  # ---- Output handling -------------------------------------------------------
  # Base output folder for LAI maps
  lai_folder_directory <-  here::here("output", "lai_maps")
  if (!dir.exists(lai_folder_directory)) dir.create(lai_folder_directory,
                                                    recursive = TRUE)
  
  # Subfolder by platform (e.g. "s2", "uas")
  output_foldername <- file.path(lai_folder_directory, platform)
  if (!dir.exists(output_foldername)) dir.create(output_foldername, 
                                                 recursive = TRUE)
  
  # Full output file name: <output/lai_maps/platform/date.tif>
  output_filename <- file.path(output_foldername, paste0(image_date, ".tif"))
  
  # Write the LAI raster to disk
  terra::writeRaster(lai_inversion, filename = output_filename, overwrite = TRUE)
  
  # Log completion
  message(output_filename, " done")
}

# ---- Cosine-similarity LUT inversion ----------------------------------------
# Returns a function that, for each spectrum v (e.g. pixel), finds the LUT
# entry in T_var with maximum cosine similarity and returns its LAI.
maximize_cos_wrapper <- function(T_var, lai) {
  
  # Precompute LUT norms once for efficiency
  T_mat   <- as.matrix(T_var)                     # ensure numeric matrix (LUT_size x n_bands)
  T_norms <- sqrt(rowSums(T_mat^2))               # Euclidean norm of each LUT spectrum
  
  function(v) {
    # v: numeric vector of band reflectances for a single pixel
    
    # Flatten to numeric vector and handle missing data
    vect <- unlist(v)                             # ensure simple numeric vector
    if (any(is.na(vect))) return(NaN)             # return NaN if any band is NA
    
    # Compute cosine similarity between v and each LUT spectrum:
    # cos(theta) = (T · v) / (||T|| * ||v||)
    v_norm <- sqrt(sum(vect^2))                   # norm of pixel spectrum
    if (v_norm == 0) return(NaN)                  # avoid division by zero
    
    cov    <- T_mat %*% vect                      # dot product for each LUT row
    cos    <- as.numeric(cov) / (T_norms * v_norm) # cosine similarity for all LUT entries
    
    # Index of LUT spectrum with maximum cosine similarity
    max_cos <- which.max(cos)
    
    # Return corresponding LAI value from LUT
    return(lai[max_cos])
  }
}

# ---- Start LAI inversion job -----------------------------------------------

# set and read shapefile path used to crop the S2 image
roi <- sf::st_read(here::here("data", "field_perimeter.gpkg"))  # expects a GPKG with field polygon

# # Start the actual job


# Loop through platform
platforms = list(list(platform="s2", harmonized = FALSE, normalized = FALSE, band_index = c(2:5, 8), aggregate = FALSE, agg_factor=NA,
                      data_resampling_matrix = data.frame(wl=c(492, 560, 664, 704, 833),
                                                          fwhm=c(66, 36, 31, 15, 106))),

                 list(platform="s2_superresolution", harmonized = TRUE, normalized = FALSE, band_index = c(1:4, 7), aggregate = FALSE, agg_factor=NA,
                      data_resampling_matrix = data.frame(wl=c(492, 560, 664, 704, 833),
                                                          fwhm=c(66, 36, 31, 15, 106))),

                 list(platform="uas", harmonized = TRUE, normalized = TRUE, band_index = 1:5, aggregate = TRUE, agg_factor=105, #GSD = 1 m
                      data_resampling_matrix = data.frame(wl=c(475, 560, 668, 717, 842),
                                                          fwhm=c(32, 27, 14, 12, 57)))
)

# platforms = list(list(platform="UAV", harmonized = TRUE, normalized = TRUE, band_index = 1:5, aggregate = TRUE, agg_factor=105, #GSD = 1 m
#                       data_resampling_matrix = data.frame(wl=c(475, 560, 668, 717, 842),
#                                                           fwhm=c(32, 27, 14, 12, 57)))
# )

for(platform_ID in platforms){
  #unlist
  platform = platform_ID$platform
  harmonized = platform_ID$harmonized
  band_index = platform_ID$band_index
  data_resampling_matrix = platform_ID$data_resampling_matrix
  normalized = platform_ID$normalized
  aggregate = platform_ID$aggregate
  agg_factor = platform_ID$agg_factor

  # find list of images
  folder_directory <- file.path(here::here("data"), platform)
  image_list <- list.files(folder_directory, recursive = FALSE)
  image_list <- image_list[endsWith(image_list, '.tif')]
  image_date_list <- strsplit(image_list, ".tif")

  # read all dates available for images
  image_date_vector <- c()
  for (i in 1:length(image_date_list)){
    j <- tail(image_date_list[[i]])
    image_date_vector <- c(image_date_vector, j)
  }


  message("Start ", platform_ID$platform," loop")
  for(this_image in image_date_vector){
    lapply(this_image,
           retrieve_lai,
           folder_directory = folder_directory,
           roi,
           spect,
           data_resampling_matrix,
           band_index=band_index,
           platform=platform,
           harmonized=harmonized,
           normalized=normalized,
           aggregate=aggregate,
           agg_factor=agg_factor)
  }
}
