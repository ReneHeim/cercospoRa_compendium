###############################################################################
# Project:     cercospoRa
# Script:      00b_lai_map_validation.R
# Author:      Rene Heim
# Affiliation: University of Goettingen
# Date:        2025-12-28
#
# Description:
#   This script validates LAI maps derived from UAS, Sentinel-2, and
#   super-resolved Sentinel-2 (S2S) data against the corresponding
#   multispectral imagery. It
#     - imports LAI and multispectral rasters for multiple dates,
#     - rescales Sentinel-2/S2S DNs to surface reflectance,
#     - harmonizes CRS, extent, and resolution within each LAI–imagery pair,
#     - renames bands to physically meaningful names,
#     - computes the CIred-edge index (CIre),
#     - evaluates structural similarity between LAI and CIre using SSIM/SIP, and
#     - exports aligned rasters (imagery, LAI, CIre, SIP) as GeoTIFFs
#       for visualization and further analysis (e.g., in QGIS).
#
# Reproducibility:
#   - Relative paths are handled with the here package (here::here()).
#   - All processing steps are scripted and can be rerun end-to-end.
#   - Optionally, append sessioninfo::session_info() at the end to
#     record package versions.
#
# Usage:
#   Run the entire script or step-by-step using RStudio sections.
#   Ensure the directory structure matches the cercospoRa compendium:
#   data/uas, data/s2, data/s2_superresolution, and output/lai_maps/*.
#
# Notes:
#   - Input data: https://github.com/ReneHeim/cercospoRa_compendium/tree/main/data
#   - Output GeoTIFFs are written to the 'output' directory with
#     source/date-specific prefixes.
#
###############################################################################


library(terra)
library(SSIMmap)
library(ggplot2)
library(here)

# Change to your root directory
home <- here::here()

# Get UAS + LAI Raster
uas1  <- rast(here("data", "uas", "20220614_uas.tif")) 
uas2 <- rast(here("data", "uas", "20220628_uas.tif"))
lai_uas1 <- rast(here("output", "lai_maps", "uas", "20220614_uas.tif")) 
lai_uas2 <- rast(here("output", "lai_maps", "uas", "20220628_uas.tif")) 

# Get S2 + LAI Raster

s21  <- rast(here("data", "s2", "20220610_s2.tif")) 
s22 <- rast(here("data", "s2", "20220630_s2.tif"))  
lai_s21 <- rast(here("output", "lai_maps", "s2", "20220610_s2.tif")) 
lai_s22 <- rast(here("output", "lai_maps", "s2", "20220630_s2.tif")) 

# S2S + LAI Raster

s2s1  <- rast(here("data", "s2_superresolution", "20220610_s2s.tif")) 
s2s2 <- rast(here("data", "s2_superresolution", "20220630_s2s.tif")) 
lai_s2s1 <- rast(here("output", "lai_maps", "s2_superresolution", "20220610_s2s.tif")) 
lai_s2s2 <- rast(here("output", "lai_maps", "s2_superresolution", "20220630_s2s.tif")) 


# Sentinel‑2 reflectance scaling (DN → reflectance)
s21  <- s21  / 10000
s22  <- s22  / 10000
s2s1 <- s2s1 / 10000
s2s2 <- s2s2 / 10000


# Harmonize function

check_same_grid <- function(r_list) {
  crs_ref <- crs(r_list[[1]])
  ext_ref <- ext(r_list[[1]])
  res_ref <- res(r_list[[1]])
  
  for (i in seq_along(r_list)) {
    if (crs(r_list[[i]]) != crs_ref)
      stop("CRS mismatch at element ", i)
    if (!all.equal(ext(r_list[[i]]), ext_ref))
      stop("Extent mismatch at element ", i)
    if (!all.equal(res(r_list[[i]]), res_ref))
      stop("Resolution mismatch at element ", i)
  }
  
  list(crs = crs_ref, ext = ext_ref, res = res_ref)
}

align_pair <- function(img, lai) {
  if (prod(res(lai)) >= prod(res(img))) {
    template    <- lai
    img_aligned <- project(img, template, method = "bilinear")
    lai_aligned <- lai
  } else {
    template    <- img
    lai_aligned <- project(lai, template, method = "bilinear")
    img_aligned <- img
  }
  
  e <- intersect(ext(img_aligned), ext(lai_aligned))
  img_aligned <- crop(img_aligned, e)
  lai_aligned <- crop(lai_aligned, e)
  
  out <- list(img = img_aligned, lai = lai_aligned)
  
  grid_info <- check_same_grid(out)
  
  cat("\nCommon grid:\n")
  cat("CRS:\n", as.character(grid_info$crs), "\n\n")
  cat("Extent:\n", as.vector(grid_info$ext), "\n\n")
  cat("Resolution:\n", grid_info$res, "\n\n")
  
  out
}


uas1_aligned  <- align_pair(uas1,  lai_uas1)
uas2_aligned  <- align_pair(uas2,  lai_uas2)

s21_aligned   <- align_pair(s21,   lai_s21)
s22_aligned   <- align_pair(s22,   lai_s22)

s2s1_aligned  <- align_pair(s2s1,  lai_s2s1)
s2s2_aligned  <- align_pair(s2s2,  lai_s2s2)




uas_bands <- c("blue", 
               "green", 
               "red", 
               "re", 
               "nir", 
               "tir")

s2_bands <- c(
  "coastal",
  "blue",
  "green",
  "red",
  "re1",
  "re2",
  "re3",
  "nir"
)

s2s_bands <- c(
  "blue",
  "green",
  "red",
  "re1",
  "re2",
  "re3",
  "nir",
  "8a",
  "b11",
  "b12"
)

names(uas1_aligned$img) <- uas_bands
names(uas2_aligned$img) <- uas_bands

names(s21_aligned$img)  <- s2_bands
names(s22_aligned$img)  <- s2_bands
names(s2s1_aligned$img) <- s2s_bands
names(s2s2_aligned$img) <- s2s_bands


# 1. Compute a vegetation index (example: CIred-edge)
# adapt band names or indices (e.g. "red", "nir") to your data

add_CIre <- function(pair, nir_name, re_name) {
  img <- pair$img
  nir <- img[[nir_name]]
  re  <- img[[re_name]]
  
  cire <- (nir / re) - 1
  names(cire) <- "CIre"
  
  pair$CIre <- cire
  pair
}

uas1_aligned <- add_CIre(uas1_aligned, nir_name = "nir",    re_name = "re")
uas2_aligned <- add_CIre(uas2_aligned, nir_name = "nir",    re_name = "re")

s21_aligned  <- add_CIre(s21_aligned,  nir_name = "nir",  re_name = "re1")
s22_aligned  <- add_CIre(s22_aligned,  nir_name = "nir",  re_name = "re1")
s2s1_aligned <- add_CIre(s2s1_aligned, nir_name = "nir",  re_name = "re1")
s2s2_aligned <- add_CIre(s2s2_aligned, nir_name = "nir",  re_name = "re1")


# Helper to add SSIM/SIP to a pair

add_SSIM_SIP <- function(pair, w = 21) {
  
  lai <- pair$lai
  max_w <- min(2 * nrow(lai) - 1, 2 * ncol(lai) - 1)
  if (w > max_w) {
    w <- max_w
    message("Adjusted SSIM window to w = ", w, " for this pair")
  }
  
  lai  <- pair$lai
  cire <- pair$CIre
  
  # 1) Global metrics
  g <- ssim_raster(lai, cire, global = TRUE, w = w)
  # g is of class "SSIM" with components: SSIM, SIM, SIV, SIP 
  
  pair$ssim_global <- list(
    SSIM = g$SSIM,
    SIM  = g$SIM,
    SIV  = g$SIV,
    SIP  = g$SIP
  )
  
  # 2) Local SIP raster (SIM/SIV/SIP per pixel; we extract the SIP layer)
  loc <- ssim_raster(lai, cire, global = FALSE, w = w)
  # loc is a SpatRaster with layers "SSIM", "SIM", "SIV", "SIP" 
  
  pair$SIP_raster <- loc[["SIP"]]
  
  pair
}

uas1_aligned <- add_SSIM_SIP(uas1_aligned, w = 21)
uas2_aligned <- add_SSIM_SIP(uas2_aligned, w = 21)

s21_aligned  <- add_SSIM_SIP(s21_aligned,  w = 3)
s22_aligned  <- add_SSIM_SIP(s22_aligned,  w = 3)

s2s1_aligned <- add_SSIM_SIP(s2s1_aligned, w = 21)
s2s2_aligned <- add_SSIM_SIP(s2s2_aligned, w = 21)


write_pair <- function(pair, prefix, out_dir = "export") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Multispectral image
  writeRaster(
    pair$img,
    filename = file.path(out_dir, paste0(prefix, "_img.tif")),
    filetype = "GTiff",
    overwrite = TRUE
  )
  
  # LAI
  writeRaster(
    pair$lai,
    filename = file.path(out_dir, paste0(prefix, "_lai.tif")),
    filetype = "GTiff",
    overwrite = TRUE
  )
  
  # CIre
  writeRaster(
    pair$CIre,
    filename = file.path(out_dir, paste0(prefix, "_cire.tif")),
    filetype = "GTiff",
    overwrite = TRUE
  )
  
  # SIP raster
  writeRaster(
    pair$SIP_raster,
    filename = file.path(out_dir, paste0(prefix, "_sip.tif")),
    filetype = "GTiff",
    overwrite = TRUE
  )
}

write_pair(uas1_aligned,  "uas_20220614",  out_dir = "output")
write_pair(uas2_aligned,  "uas_20220628",  out_dir = "output")

write_pair(s21_aligned,   "s2_20220610",   out_dir = "output")
write_pair(s22_aligned,   "s2_20220630",   out_dir = "output")

write_pair(s2s1_aligned,  "s2s_20220610",  out_dir = "output")
write_pair(s2s2_aligned,  "s2s_20220630",  out_dir = "output")


# Reproducibility: record session info to file
if (!requireNamespace("sessioninfo", quietly = TRUE)) {
  install.packages("sessioninfo")
}
sessioninfo::session_info() |> 
  capture.output(file = here::here("output", "session_info_00b_lai_map_validation.txt"))

# OPTIONAL

# # 2. Basic variability: mean, sd, CV
# lai_vals  <- values(lai1, na.rm = TRUE)
# cire_vals <- values(cire1, na.rm = TRUE)
# 
# lai_mean <- mean(lai_vals);  lai_sd <- sd(lai_vals)
# lai_cv   <- lai_sd / lai_mean
# 
# cire_mean <- mean(cire_vals); cire_sd <- sd(cire_vals)
# cire_cv   <- cire_sd / cire_mean
# 
# c(LAI_mean = lai_mean, LAI_sd = lai_sd, LAI_CV = lai_cv,
#   CIre_mean = cire_mean, CIre_sd = cire_sd, CIre_CV = cire_cv)
# 
# # 3. Rank relationship (should reflect similar spatial pattern)
# cor_spearman <- cor(lai_vals, cire_vals, method = "spearman",
#                     use = "complete.obs")
# cor_spearman
# 
# # Semivario
# 
# # Make a two‑layer stack
# stk <- c(lai1, cire1)
# names(stk) <- c("LAI", "CIre")
# 
# # Extract coordinates + values as a data frame
# df <- as.data.frame(stk, xy = TRUE, na.rm = TRUE)
# names(df) <- c("x", "y", "LAI", "CIre")
# 
# # Thin for speed, e.g. every 5th point
# #thin_factor <- 5
# #df_thin <- df[seq(1, nrow(df), by = thin_factor), ]
# 
# # STANDARDIZE (z-score)
# 
# df$LAI_z  <- scale(df$LAI)
# df$CIre_z <- scale(df$CIre)
# 
# #colnames(df_thin) <- c("x", "y", "LAI", "CIre", "LAI_norm", "CIre_norm")
# 
# #Treat x and y as coordinates
# 
# coordinates(df) <- ~ x + y
# 
# 
# 
# # choose a reasonable cutoff: e.g. half of max field dimension
# field_width  <- diff(range(df@coords[,1]))
# field_height <- diff(range(df@coords[,2]))
# cutoff_dist  <- 0.5 * max(field_width, field_height)
# 
# vg_lai  <- variogram(LAI_z  ~ 1, df, cutoff = cutoff_dist)
# vg_cire <- variogram(CIre_z ~ 1, df, cutoff = cutoff_dist)
# 
# par(mfrow = c(2,1))
# plot(vg_lai,  main = "Semivariogram – LAI")
# plot(vg_cire, main = "Semivariogram – CIre")
# 
# 
# # model
# 
# model_ini <- vgm(psill = 1.2, model = "Lin", range = 80, nugget = 0.2)
# fit_lai   <- fit.variogram(vg_lai,  model_ini)
# fit_cire  <- fit.variogram(vg_cire, model_ini)
# 
# fit_lai
# fit_cire
# 
# # assuming vg_lai and vg_cire share the same 'dist' bins
# rmse_gamma <- sqrt(mean((vg_lai$gamma - vg_cire$gamma)^2))
# rmse_gamma
# #####################################
