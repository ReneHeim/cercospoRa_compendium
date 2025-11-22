#library(readr)
library(sf)
#library(raster)
library(terra)
library("hsdar") # For RTM modeling

# setwd("E:\\Clean_Directory")
tail <- function(x){
  y <- x[length(x)]
  return(y)
}
#generate look up table (LUT)
#define parameters for sugar beet according to
#Jay, Sylvain, Fabienne Maupas, Ryad Bendoula, and Nathalie Gorretta.
#"Retrieving LAI, chlorophyll and nitrogen contents in sugar beet crops
#from multi-angular optical remote sensing: Comparison of vegetation indices
#and PROSAIL inversion for field phenotyping." Field Crops Research 210 (2017): 33-46.
tail <- function(x){
  y <- x[length(x)]
  return(y)
}


solar_zenith = 10
deltaazimuth = 180
sensor_zenith = 0
LAI0 = seq(0.1,6,0.01)
N0 = seq(1,3.5,0.5)
Cab0 = seq(20,65,7)
Cw0 = seq(0.03,0.09,0.015)

LAI =  rep(LAI0, times = length(Cab0)*length(Cw0)*length(N0))
Cab = rep(rep(Cab0, each = length(LAI0)), times=length(Cw0)*length(N0))
Cw = rep(rep(Cw0, each = length(Cab0)*length(LAI0)), times = length(N0))
N = rep(N0, each = length(Cab0)*length(LAI0)*length(Cw0))

param = data.frame(LAI = LAI)
spect=PROSAIL(parameterList = param,
              psi=deltaazimuth,
              tts=solar_zenith,
              tto=sensor_zenith,
              N=N,
              Cab=Cab,
              Cw=Cw
)

# spect <- prosail::PRO4SAIL(lai = LAI[1],
#                           psi=deltaazimuth,
#                           tts=solar_zenith,
#                           tto=sensor_zenith,
#                           N=N,
#                           CHL=Cab,
#                           EWT=Cw)




#spectral bands captured by S2 corresponding to the Micasense Altum
data_resampling_matrix <- data.frame(center=c(492, 560, 664, 704, 833),
                                     fwhm=c(66, 36, 31, 15, 106))


# function to retrieve LAI values at specific points for a given UAV image
#  containing reflectance values
retrieve_lai <- function(image_date,
                         folder_directory,
                         roi,
                         spect,
                         data_resampling_matrix,
                         band_index = c(2:5, 8),
                         platform = "S2",
                         harmonized = FALSE,
                         normalized = FALSE,
                         aggregate = FALSE,
                         agg_factor = 50){

  # Resample spectrum according to sensor
  spect_resample <- spectralResampling(spect, data_resampling_matrix)

  # read raster
  image_path <- file.path(folder_directory, paste(image_date,'.tif', sep = ''))
  allbands <- terra::rast(image_path)
  allbands <- allbands[[band_index]] #allbands[[c(2:5, 8)]]

  # Preprocess the rasters
  cat(paste0("convert to SpatRaster: ", image_date, "\n"))
  allbands <- terra::crop(allbands, st_transform(roi, crs(allbands)), mask=TRUE)
  if(!harmonized) allbands <- allbands - 1000 #additive offset
  if(!normalized) allbands <- (allbands)/10000 # harmonized S2
  if(!aggregate) allbands <- terra::aggregate(allbands, fact=agg_factor, fun='mean') # aggregate in case it is UAV


  # Extract spectrum and LAI
  spect_resample_ <- hsdar::spectra(spect_resample)
  lai <- hsdar::SI(spect_resample)$LAI

  # Apply the function
  lai_inversion = terra::app(allbands, maximize_cos_wrapper(T_var=spect_resample_, lai=lai), cores = 8)

  # Create saving path
  lai_folder_directory <- file.path("Data", paste0("LAI_maps"))
  if(!dir.exists(lai_folder_directory)) dir.create(lai_folder_directory)

  output_foldername <- file.path(lai_folder_directory, platform)
  if(!dir.exists(output_foldername)) dir.create(output_foldername)
  output_filename <- file.path(output_foldername, paste(image_date,'.tif', sep = ''))
  terra::writeRaster(lai_inversion, filename = output_filename)
  print(paste(output_filename, 'done',  sep = ' '))
}

# cosine distance
maximize_cos_wrapper <- function(T_var, lai) {
  function(v) {
    vect <- unlist(v)
    if(any(is.na(vect))) return(NaN)
    cov <- as.matrix(T_var) %*% vect
    sqrvar1 <- sqrt(sum(vect^2))
    sqrvar2 <- sqrt(apply(as.matrix(T_var)^2, 1, sum))
    cos <- cov / (sqrvar1 * sqrvar2)
    max_cos <- which.max(cos)
    return(lai[max_cos])
  }
}


# Start the actual job
# set and read shapefile path used to crop the S2 image
roi <- st_read("Data\\ROI\\Field_perimeter.gpkg")


# Loop through platform
platforms = list(list(platform="S2", harmonized = FALSE, normalized = FALSE, band_index = c(2:5, 8), aggregate = FALSE, agg_factor=NA,
                      data_resampling_matrix = data.frame(center=c(492, 560, 664, 704, 833),
                                                          fwhm=c(66, 36, 31, 15, 106))),

                 list(platform="S2_superresolution", harmonized = TRUE, normalized = FALSE, band_index = c(1:4, 7), aggregate = FALSE, agg_factor=NA,
                      data_resampling_matrix = data.frame(center=c(492, 560, 664, 704, 833),
                                                          fwhm=c(66, 36, 31, 15, 106))),

                 list(platform="UAV", harmonized = TRUE, normalized = TRUE, band_index = 1:5, aggregate = TRUE, agg_factor=50, #GSD = 0.5 m
                      data_resampling_matrix = data.frame(center=c(475, 560, 668, 717, 842),
                                                           fwhm=c(32, 27, 14, 12, 57)))
                 )

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
  folder_directory <- file.path("data/", platform)
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
