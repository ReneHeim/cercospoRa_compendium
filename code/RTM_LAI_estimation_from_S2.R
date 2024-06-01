library(readr)
library(sf)
library(raster)
library("hsdar") # For RTM modeling


#function to extract last element of a vector
tail <- function(x){
  y <- x[length(x)]
  return(y)
}


#function to minimize RMSE. Takes an input vector (spectrum) 
#and a list of vectors to match the first vector to
#then selects the vector in the matching list that has the lowest
#RMSE in comparison to the input spectrum
minimize_RMSE <- function(v, T){
  rownum <- dim(T)[1]
  tab <- matrix(rep(v, each = rownum), nrow = rownum)
  diff <- (T-tab)^2
  RMSE <- apply(diff, 1, mean)
  minim_RMSE <- which.min(RMSE)
  return(minim_RMSE)
}

maximize_cos <- function(v, T){
  rownum <- dim(T)[1]
  tab <- matrix(rep(v, each = rownum), nrow = rownum)
  prod <- T*tab
  cov <- apply(prod, 1, sum)
  sqrvar1 <- sqrt(apply(tab^2, 1, sum))
  sqrvar2 <- sqrt(apply(T^2, 1, sum))
  cos <- cov/(sqrvar1*sqrvar2)
  max_cos <- which.max(cos)
  return(max_cos)
}



# function to retrieve LAI values at specific points for a given UAV image containig reflectance values
retrieve_lai <- function(image_date, folder_directory, sf_bonitur){
  #generate look up table (LUT)
  #define parameters for sugar beet according to 
  #Jay, Sylvain, Fabienne Maupas, Ryad Bendoula, and Nathalie Gorretta. 
  #"Retrieving LAI, chlorophyll and nitrogen contents in sugar beet crops 
  #from multi-angular optical remote sensing: Comparison of vegetation indices 
  #and PROSAIL inversion for field phenotyping." Field Crops Research 210 (2017): 33-46.
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
  
  #spectral bands captured by S2 corresponding to the Micasense Altum
  data_resampling_matrix <- data.frame(center=c(492, 560, 664, 704, 833),
                                       fwhm=c(66, 36, 31, 15, 106))
  
  #Resample spectra from the LUT  to spectral bands captured by S2 sensors
  spect_resample <- spectralResampling(spect, data_resampling_matrix)
  
  
  #define dataframe that will contain LAI after extracting spectra extracted from S2 images at each
  #ground sampling point according to the respective date.
  Spectra_extraction <- data.frame()
  
  # read raster  
  image_path <- file.path(folder_directory, paste(image_date,'.tif', sep = ''))
  allbands <- brick(image_path)
  allbands <- allbands[[c(2:5, 8)]]
  
  #Crop S2 image to our region of interest (defined by a polygon shapefile) to reduce computational burden
  sf_this_spectra <- raster::crop(allbands, sf_bonitur)
  sf_this_spectra <- (sf_this_spectra-1000)/10000
  
  #extract LAI at each sampling location and store it
  #lai_RMSE <- c()
  lai_cos <- sf_this_spectra[[1]]#lai_COS <- c()
  for (j in 1:dim(sf_this_spectra)[1]){
    for (i in 1:dim(sf_this_spectra)[2]){
      spect <- sf_this_spectra[j,i]
      #lai_i_RMSE <- SI(spect_resample)$LAI[minimize_RMSE(spect, spectra(spect_resample))] # this in case you want the RMSE distance
      lai_i_COS <- SI(spect_resample)$LAI[maximize_cos(spect, spectra(spect_resample))]
      if(is.na(lai_i_COS[1])) lai_i_COS <- NA
      #lai_RMSE <- c(lai_RMSE, lai_i_RMSE)
      #lai_COS <- c(lai_COS, lai_i_COS)
      lai_cos[j,i]  <- lai_i_COS
    }  
  }
  
  #return(Spectra_extraction)
  output_filename <- file.path(folder_directory, 'lai', paste(image_date,'.tif', sep = ''))
  writeRaster(lai_cos, filename = output_filename, format = "GTiff")
  print(paste(output_filename, 'done',  sep = ' '))
}





# Start the actual job
# set and read shapefile path used to crop the S2 image
sf_path <- r"(D:\cercospoRa\raw\scoring\SE04_field_perimeter.gpkg)"
sf_field_perimeter <- read_sf(sf_path)

# find list of images
folder_directory <- r"(D:\cercospoRa\raw\S2A)"
image_list <- list.files(folder_directory, recursive = FALSE)
image_list <- image_list[endsWith(image_list, '.tif')]
image_date_list <- strsplit(image_list, ".tif")
  
# read all dates available for images
image_date_vector <- c()
for (i in 1:length(image_date_list)){
    j <- tail(image_date_list[[i]])
    image_date_vector <- c(image_date_vector, j)
}


library(parallel) # for parallel computing

# Calculate the number of cores available for parallel computing
no_cores <- detectCores() - 1
no_cores <- min(5, no_cores)

# Initialize cluster for parallel computing
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=c("retrieve_lai", "tail", "maximize_cos", "sf_field_perimeter"))

parLapply(cl, 
          image_date_vector,
          function(x){
            library(readr)
            library(sf)
            library(raster)
            library("hsdar") # For RTM modeling
            folder_directory <- r"(D:\cercospoRa\raw\S2A)"
            retrieve_lai(x, folder_directory, sf_field_perimeter)
            })
stopCluster(cl)
