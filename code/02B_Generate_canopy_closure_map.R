library(cercospoRa)
library(terra)

setwd("E:\\Clean_Directory\\Data")

base_dir <- file.path("LAI_maps")
destination <- "CC_maps"

if (!dir.exists(destination)) {
  dir.create(destination)
}

generate_CC <- function(base_dir, platform){
  # function from Cercospora to read all available LAI maps at different dates and stack them together in a list
  epidemic_onset_param <-
    read_sb_growth_parameter(
      list.files(file.path(base_dir, platform),pattern = "tif",
                 full.names = TRUE),
      img_dates = as.POSIXct(as.Date(list.files(file.path(base_dir, platform),pattern = "tif",
                                                full.names = FALSE), 
                                     format = "%Y_%m_%d"),
                             tz = "UTC"),
      10)
  
  # Use the previously generated list to generate a map of parameters for the pixel-wise sigmoids
  param_rxt <- calc_r_x0(epidemic_onset_param,
                         min_r = 0.02,
                         max_r = 0.05,
                         k = 6)
  
  # create CC maps
  c_closure <- calc_c_closure(param_rxt,
                              x1 = 1.3,
                              k=6 )
  
  # Turn into days after sowing
  c_closure <- c_closure - as.numeric(as.Date("2022-03-31"))
  
  # create destination folder
  if(!dir.exists(file.path(destination))) dir.create(file.path(destination), recursive = TRUE)
  
  names(c_closure) <- "Canopy_closure"
  writeRaster(c_closure,
              file.path(destination,
                        paste(platform, '_CC.tif', sep = '')), 
              overwrite=TRUE)
}

# loop through platforms with the function
platforms <- c("S2", "S2_superresolution", "UAV")

for(platform in platforms){
  generate_CC(base_dir, platform)
}
