library(cercospoRa)
library(terra)

setwd("E:\\Clean_Directory\\Data")
  
base_dir <- file.path("LAI_maps") # file.path("LAI_maps")
destination <- "LAI_progression_curve"

if (!dir.exists(destination)) {
  dir.create(destination)
}

plot_r <- function(base_dir, platform){
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
  
  # Extract sigmoid parameter
  r <- param_rxt$r
  x0 <- param_rxt$x0
  t0 <- param_rxt$t0
  
  # put parameters in a table
  table <- data.frame(r=matrix(r), x=matrix(x0))
  table <- na.omit(table)
  
  # create destination folder
  if(!dir.exists(file.path(destination))) dir.create(file.path(destination), recursive = TRUE)
  
  png(file.path(destination, paste(platform, '_growth_curve.png', sep = '')), 
      units="in", 
      width=6, 
      height=5, 
      res=600)
  
  r_1 <- table$r[1]
  x_1 <- table$x[1]
  
  plot_sigma <- function(t,r_i, x_i, ti_0){
    t <- as.numeric(t)
    ti_0 <- as.numeric(as.Date(ti_0))
    N <- 6/(1 + ((6-x_i)/x_i)*exp(-r_i*(t-ti_0)))
    return(N)
  }
  
  t_i <- seq(as.Date("2022-04-01"), as.Date("2022-10-15"), 1)
  N_i <- plot_sigma(t_i, r_1, x_1, t0)
  par(mgp = c(2.5, 0.5, 0))
  plot(t_i, N_i, t='l', xlab = 'Dates', ylab = 'LAI', ylim = c(0,6), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2), cex.axis=1.25, cex.lab = 1.5)
  
  
  for(i in 2:dim(table)[1]){
    r_i <- table$r[i]
    x_i <- table$x[i]
    
    
    t_i <- seq(as.Date("2022-04-01"), as.Date("2022-10-15"), 1)
    N_i <- plot_sigma(t_i, r_i, x_i, t0)
    lines(t_i, N_i, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
  }
  dev.off()
}

# loop through platforms with the function
platforms <- c("S2", "S2_superresolution", "UAV")

for(platform in platforms){
  plot_r(base_dir, platform)
}
