library(terra)
library(here)

# Setup and detect “best” CRS

root_dir <- here("output", "lai_maps")  # directory that contains the 3 subfolders

# list all raster files in all subdirectories
ras_files <- list.files(root_dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# inspect CRSs
crs_list <- sapply(ras_files, function(f) {
  r <- rast(f)
  crs(r)              # WKT string or empty if undefined
})

# tabulate CRSs
table(crs_list)

# choose the most frequent defined CRS as target
defined <- crs_list[crs_list != ""]
target_crs <- names(sort(table(defined), decreasing = TRUE))[1]
target_crs

# Batch reproject and write _repro rasters

for (f in ras_files) {
  r <- rast(f)
  this_crs <- crs(r)
  
  # skip if already in target CRS
  if (this_crs == target_crs || this_crs == "") {
    # you might want to skip undefined CRS or handle them separately
    next
  }
  
  # reproject
  r_repro <- project(r, target_crs)  # may take time for big rasters[web:131][web:128]
  
  # output path: same folder
  out_f <- sub("\\.tif$", ".tif", f)
  
  writeRaster(r_repro, out_f, overwrite = TRUE)
}

# Check that reprojection succeeded

repro_files <- list.files(root_dir, pattern = ".tif$", full.names = TRUE, recursive = TRUE)

check_crs <- sapply(repro_files, function(f) {
  r <- rast(f)
  crs(r) == target_crs
})

all(check_crs)
which(!check_crs)  # any that failed

