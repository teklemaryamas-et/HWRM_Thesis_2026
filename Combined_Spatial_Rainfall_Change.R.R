# ============================================================
# COMBINE ALL SEASONAL GeoTIFF MAPS IN ONE RUN
# Seasons:
# 1. Annual
# 2. Belg_MAM
# 3. Kiremt_JJAS
# 4. Bega_ONDJF
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c("terra", "ggplot2", "dplyr", "patchwork")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT AND OUTPUT FOLDERS
# ============================================================

tif_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Spatial_Rainfall_Change/GeoTIFF"
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Spatial_Rainfall_Change/Combined_Figures"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ============================================================
# 3. DEFINE ALL SEASONS
# ============================================================

season_list <- c("Annual", "Belg_MAM", "Kiremt_JJAS", "Bega_ONDJF")

# ============================================================
# 4. FUNCTION TO CONVERT RASTER TO DATA FRAME
# ============================================================

raster_to_df <- function(raster_file) {
  r <- rast(raster_file)
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "change"
  return(df)
}

# ============================================================
# 5. FUNCTION TO MAKE ONE MAP
# ============================================================

make_map <- function(df, panel_title, color_limit) {
  ggplot(df, aes(x = x, y = y, fill = change)) +
    geom_raster() +
    coord_equal() +
    scale_fill_gradient2(
      low = "brown3",
      mid = "white",
      high = "darkgreen",
      midpoint = 0,
      limits = c(-color_limit, color_limit),
      name = "Change (%)"
    ) +
    labs(
      title = panel_title,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      panel.grid = element_line(color = "grey90"),
      legend.position = "right"
    )
}

# ============================================================
# 6. LOOP THROUGH ALL SEASONS
# ============================================================

for (season_name in season_list) {
  
  cat("\nProcessing season:", season_name, "\n")
  
  # File names
  file_ssp245_near <- file.path(
    tif_folder,
    paste0("MIROC6_ssp245_Near_future_2030_2050_", season_name, "_percentage_change.tif")
  )
  
  file_ssp245_far <- file.path(
    tif_folder,
    paste0("MIROC6_ssp245_Far_future_2051_2080_", season_name, "_percentage_change.tif")
  )
  
  file_ssp585_near <- file.path(
    tif_folder,
    paste0("MIROC6_ssp585_Near_future_2030_2050_", season_name, "_percentage_change.tif")
  )
  
  file_ssp585_far <- file.path(
    tif_folder,
    paste0("MIROC6_ssp585_Far_future_2051_2080_", season_name, "_percentage_change.tif")
  )
  
  files <- c(file_ssp245_near, file_ssp245_far, file_ssp585_near, file_ssp585_far)
  
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    cat("Skipping season", season_name, "because some files are missing:\n")
    print(missing_files)
    next
  }
  
  # Read rasters
  df1 <- raster_to_df(file_ssp245_near)
  df2 <- raster_to_df(file_ssp245_far)
  df3 <- raster_to_df(file_ssp585_near)
  df4 <- raster_to_df(file_ssp585_far)
  
  # Common color scale for the current season
  all_values <- c(df1$change, df2$change, df3$change, df4$change)
  color_limit <- max(abs(all_values), na.rm = TRUE)
  
  # Make maps
  p1 <- make_map(df1, "SSP2-4.5 | Near future 2030-2050", color_limit)
  p2 <- make_map(df2, "SSP2-4.5 | Far future 2051-2080", color_limit)
  p3 <- make_map(df3, "SSP5-8.5 | Near future 2030-2050", color_limit)
  p4 <- make_map(df4, "SSP5-8.5 | Far future 2051-2080", color_limit)
  
  # Combine maps
  combined_map <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = paste0("Projected Spatial Change in ", season_name, " Rainfall over the Upper Blue Nile Basin"),
      subtitle = "PT-corrected MIROC6 projections relative to the 1994-2014 XGBoost-corrected CHIRPS baseline",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  print(combined_map)
  
  # Save output
  output_png <- file.path(
    output_folder,
    paste0("Combined_", season_name, "_Spatial_Rainfall_Change.png")
  )
  
  ggsave(
    filename = output_png,
    plot = combined_map,
    width = 13,
    height = 10,
    dpi = 300
  )
  
  cat("Saved:", output_png, "\n")
}

cat("\nAll available seasonal figures have been generated successfully.\n")

