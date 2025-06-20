library(cercospoRa)
library(terra)
wethr <- read.csv(system.file("extdata", "clean_weather.csv",
                              package = "cercospoRa"))
wethr <- format_weather(wethr,time_zone = "UTC")

#setwd("E:\\Clean_Directory")


generate_EO_maps = function(platform){
  # Directory containing LAI maps
  img_dir <- file.path("data", platform)

  # Function from Cercospora to read all available LAI maps at different dates
  #  and stack them together in a list
  epidemic_onset_param <-
    read_sb_growth_parameter(
      list.files(img_dir,
                 pattern = "\\.tif$",
                 full.names = TRUE),
      img_dates = as.POSIXct(as.Date(list.files(img_dir,pattern = "\\.tif$",
                                                full.names = FALSE),
                                     format = "%Y_%m_%d"),
                             tz = "UTC"),
      10)

  # Use the previously generated list to generate a map of parameters for the pixel-wise sigmoids
  param_rxt <- calc_r_x0(epidemic_onset_param,
                         min_r = 0.02,
                         max_r = 0.05,
                         k = 6)

  # Determine canopy closure dates (date when the sigmoid crosse LAI = 3)
  c_closure <- calc_c_closure(param_rxt,
                              x1 = 1.3,
                              k=6 )

  # Partially fill the funtion to calculate epidemic onset dates from the previosly generated CC dates,
  # keeping cultivar susceptibility as the only variable parameter
  calc_epidemic_onset_from_image_different_sus <- function(cultivar_sus){
    # this takes about 20 sec to run
    epidemic_onset_map <- calc_epidemic_onset_from_image(start = as.POSIXct("2022-04-25",tz = "UTC"),
                                                         end = as.POSIXct("2022-09-30",tz = "UTC"),
                                                         cc_r = c_closure,
                                                         weather = wethr,
                                                         cultivar_sus = cultivar_sus)

    sowing_date <- as.numeric(as.POSIXct("2022-03-31",tz = "UTC"))
    epidemic_onset_map <- (epidemic_onset_map-sowing_date)/(3600*24)

    # Output file name and folder
    file_path <- file.path("Data\\Predicted_epidemic_onset",
                           paste0(platform),
                           paste0("Ep_onset_", cultivar_sus, ".tif"))

    if(!dir.exists(file.path("Data\\Predicted_epidemic_onset",
                             paste0(platform)))){
      dir.create(file.path("Data\\Predicted_epidemic_onset",
                           paste0(platform)),
                 recursive = TRUE)}

    # Save the generated epidemic onset map
    writeRaster(epidemic_onset_map, file_path)

    # plot the generated epidemic onset map
    terra::plot(epidemic_onset_map)
  }

  # use the previous function for different culvitar susceptbility values
  #(choose the right value for your cultivar)
  cultivar_sus <- seq(from = 7, to = 12, by = 1)

  sapply(cultivar_sus, calc_epidemic_onset_from_image_different_sus)
}

# loop through platforms with the function
platforms <- c("s2", "s2s", "uas")

for(platform in platforms){
  generate_EO_maps(platforms)
}
