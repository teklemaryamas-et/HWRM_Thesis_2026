# ============================================================
# HIGH-QUALITY SPATIAL RAINFALL CHANGE MAPS
# WITH BASIN BOUNDARY + NMA STATION DOTS
# Study area: Upper Blue Nile Basin, Ethiopia
# GCM: MIROC6
# Scenarios: SSP2-4.5 and SSP5-8.5
# Periods: 2030–2050 and 2051–2080
# Baseline: 1994–2014 XGBoost-corrected CHIRPS
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c("terra", "sf", "ggplot2", "dplyr", "patchwork")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT AND OUTPUT FOLDERS
# ============================================================

# Folder containing spatial percentage-change GeoTIFF files
tif_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Spatial_Rainfall_Change/GeoTIFF"

# Output folder
output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Spatial_Rainfall_Change/Combined_Figures_With_Stations"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ============================================================
# 3. BASIN SHAPEFILE PATH
# ============================================================

# Your basin shapefile folder
basin_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/abay"

# Automatically find the .shp file inside the basin folder
basin_shp_files <- list.files(
  path = basin_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(basin_shp_files) == 0) {
  warning("No basin shapefile found in basin_folder. Maps will be created without basin boundary.")
  basin_sf <- NULL
} else {
  basin_shapefile <- basin_shp_files[1]
  cat("Using basin shapefile:\n")
  cat(basin_shapefile, "\n")
  
  basin_sf <- sf::st_read(basin_shapefile, quiet = TRUE)
  basin_sf <- sf::st_transform(basin_sf, 4326)
}

# Optional river shapefile
# If you have river shapefile, paste the path. Otherwise keep NA.
river_shapefile <- NA

river_sf <- NULL

if (!is.na(river_shapefile) && file.exists(river_shapefile)) {
  river_sf <- sf::st_read(river_shapefile, quiet = TRUE)
  river_sf <- sf::st_transform(river_sf, 4326)
}

# ============================================================
# 4. NMA RAINFALL STATIONS USED IN YOUR STUDY
# Station dots will be shown, but station names will not be printed on map.
# ============================================================

stations_df <- data.frame(
  ID = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  NAME = c(
    "Assosa",
    "Bahir Dar",
    "Debre Birhan",
    "Debre Markos",
    "Dessie",
    "Debre Tabor",
    "Finote Selam",
    "Gondar",
    "Nekemte"
  ),
  Lat = c(
    10.046,
    11.595,
    9.670,
    10.326,
    11.118,
    11.867,
    10.682,
    12.521,
    9.083
  ),
  Long = c(
    34.546,
    37.360,
    39.513,
    37.739,
    39.635,
    37.995,
    37.263,
    37.432,
    36.549
  ),
  ELEVATION = c(
    1541,
    1800,
    3206,
    2446,
    2553,
    2612,
    1872,
    1973,
    2119
  )
)

write.csv(
  stations_df,
  file.path(output_folder, "NMA_Rainfall_Stations_Used_in_Study.csv"),
  row.names = FALSE
)

station_names_text <- paste(stations_df$NAME, collapse = ", ")

# ============================================================
# 5. DEFINE SEASONS FOR UPPER BLUE NILE BASIN
# ============================================================

season_info <- data.frame(
  season_code = c("Annual", "Belg_MAM", "Kiremt_JJAS", "Bega_ONDJF"),
  season_name = c(
    "Annual rainfall",
    "Belg rainfall, March-May",
    "Kiremt rainfall, June-September",
    "Bega rainfall, October-February"
  )
)

# ============================================================
# 6. FUNCTION TO READ RASTER AS DATA FRAME
# ============================================================

raster_to_df <- function(raster_file) {
  
  r <- terra::rast(raster_file)
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  
  names(df)[3] <- "change"
  
  return(df)
}

# ============================================================
# 7. FUNCTION TO MAKE ONE MAP PANEL
# ============================================================

make_map <- function(df, scenario_label, period_label, mean_change, color_limit) {
  
  panel_title <- paste0(
    scenario_label,
    " | ",
    period_label,
    "\nMean change = ",
    sprintf("%.2f", mean_change),
    "%"
  )
  
  p <- ggplot() +
    geom_raster(
      data = df,
      aes(x = x, y = y, fill = change)
    ) +
    coord_equal(expand = FALSE) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "#f7f7f7",
      high = "#1a9850",
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
      plot.title = element_text(
        size = 10.5,
        face = "bold",
        hjust = 0.5,
        lineheight = 1.05
      ),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8, color = "black"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Add basin boundary
  if (!is.null(basin_sf)) {
    p <- p +
      geom_sf(
        data = basin_sf,
        fill = NA,
        color = "black",
        linewidth = 0.8,
        inherit.aes = FALSE
      )
  }
  
  # Add river network if available
  if (!is.null(river_sf)) {
    p <- p +
      geom_sf(
        data = river_sf,
        color = "blue",
        linewidth = 0.45,
        inherit.aes = FALSE
      )
  }
  
  # Add station dots only
  p <- p +
    geom_point(
      data = stations_df,
      aes(x = Long, y = Lat),
      shape = 21,
      fill = "red",
      color = "black",
      size = 2.6,
      stroke = 0.6,
      inherit.aes = FALSE
    )
  
  return(p)
}

# ============================================================
# 8. EMPTY OUTPUT TABLES
# ============================================================

mean_summary_all <- data.frame()
summary_notes_all <- data.frame()

# ============================================================
# 9. LOOP THROUGH ALL SEASONS
# ============================================================

for (i in 1:nrow(season_info)) {
  
  season_code <- season_info$season_code[i]
  season_name <- season_info$season_name[i]
  
  cat("\n================================================\n")
  cat("Processing:", season_name, "\n")
  cat("================================================\n")
  
  # ------------------------------------------------------------
  # Input GeoTIFF files
  # ------------------------------------------------------------
  
  file_ssp245_near <- file.path(
    tif_folder,
    paste0("MIROC6_ssp245_Near_future_2030_2050_", season_code, "_percentage_change.tif")
  )
  
  file_ssp245_far <- file.path(
    tif_folder,
    paste0("MIROC6_ssp245_Far_future_2051_2080_", season_code, "_percentage_change.tif")
  )
  
  file_ssp585_near <- file.path(
    tif_folder,
    paste0("MIROC6_ssp585_Near_future_2030_2050_", season_code, "_percentage_change.tif")
  )
  
  file_ssp585_far <- file.path(
    tif_folder,
    paste0("MIROC6_ssp585_Far_future_2051_2080_", season_code, "_percentage_change.tif")
  )
  
  files <- c(
    file_ssp245_near,
    file_ssp245_far,
    file_ssp585_near,
    file_ssp585_far
  )
  
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    cat("Skipping", season_code, "because these files are missing:\n")
    print(missing_files)
    next
  }
  
  # ------------------------------------------------------------
  # Read raster values
  # ------------------------------------------------------------
  
  df1 <- raster_to_df(file_ssp245_near)
  df2 <- raster_to_df(file_ssp245_far)
  df3 <- raster_to_df(file_ssp585_near)
  df4 <- raster_to_df(file_ssp585_far)
  
  # ------------------------------------------------------------
  # Calculate statistics
  # ------------------------------------------------------------
  
  stats_fun <- function(df) {
    data.frame(
      mean_change = mean(df$change, na.rm = TRUE),
      min_change = min(df$change, na.rm = TRUE),
      max_change = max(df$change, na.rm = TRUE),
      sd_change = sd(df$change, na.rm = TRUE)
    )
  }
  
  s1 <- stats_fun(df1)
  s2 <- stats_fun(df2)
  s3 <- stats_fun(df3)
  s4 <- stats_fun(df4)
  
  season_summary <- data.frame(
    Season = season_name,
    Season_Code = season_code,
    Scenario = c("SSP2-4.5", "SSP2-4.5", "SSP5-8.5", "SSP5-8.5"),
    Period = c(
      "Near future 2030-2050",
      "Far future 2051-2080",
      "Near future 2030-2050",
      "Far future 2051-2080"
    ),
    Mean_Percentage_Change = c(
      s1$mean_change,
      s2$mean_change,
      s3$mean_change,
      s4$mean_change
    ),
    Minimum_Percentage_Change = c(
      s1$min_change,
      s2$min_change,
      s3$min_change,
      s4$min_change
    ),
    Maximum_Percentage_Change = c(
      s1$max_change,
      s2$max_change,
      s3$max_change,
      s4$max_change
    ),
    SD_Percentage_Change = c(
      s1$sd_change,
      s2$sd_change,
      s3$sd_change,
      s4$sd_change
    )
  )
  
  mean_summary_all <- rbind(mean_summary_all, season_summary)
  
  cat("\nSummary statistics for", season_name, ":\n")
  print(season_summary)
  
  # ------------------------------------------------------------
  # Common color scale for the season
  # ------------------------------------------------------------
  
  all_values <- c(df1$change, df2$change, df3$change, df4$change)
  
  color_limit <- ceiling(max(abs(all_values), na.rm = TRUE) / 5) * 5
  
  if (!is.finite(color_limit) || color_limit == 0) {
    color_limit <- 10
  }
  
  # ------------------------------------------------------------
  # Create map panels
  # ------------------------------------------------------------
  
  p1 <- make_map(
    df = df1,
    scenario_label = "SSP2-4.5",
    period_label = "Near future 2030-2050",
    mean_change = s1$mean_change,
    color_limit = color_limit
  )
  
  p2 <- make_map(
    df = df2,
    scenario_label = "SSP2-4.5",
    period_label = "Far future 2051-2080",
    mean_change = s2$mean_change,
    color_limit = color_limit
  )
  
  p3 <- make_map(
    df = df3,
    scenario_label = "SSP5-8.5",
    period_label = "Near future 2030-2050",
    mean_change = s3$mean_change,
    color_limit = color_limit
  )
  
  p4 <- make_map(
    df = df4,
    scenario_label = "SSP5-8.5",
    period_label = "Far future 2051-2080",
    mean_change = s4$mean_change,
    color_limit = color_limit
  )
  
  # ------------------------------------------------------------
  # Combine panels
  # ------------------------------------------------------------
  
  combined_map <- (p1 + p2) / (p3 + p4) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  combined_map <- combined_map +
    patchwork::plot_annotation(
      title = paste0(
        "Projected Spatial Change in ",
        season_name,
        " over the Upper Blue Nile Basin"
      ),
      subtitle = "PT-corrected MIROC6 projections relative to the 1994-2014 XGBoost-corrected CHIRPS baseline",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
      )
    )
  
  # Do not print the large plot in RStudio
  # print(combined_map)
  
  # ------------------------------------------------------------
  # Save high-quality outputs
  # ------------------------------------------------------------
  
  output_png <- file.path(
    output_folder,
    paste0("Combined_", season_code, "_Spatial_Rainfall_Change_Stations_600dpi.png")
  )
  
  output_tiff <- file.path(
    output_folder,
    paste0("Combined_", season_code, "_Spatial_Rainfall_Change_Stations_600dpi.tiff")
  )
  
  ggsave(
    filename = output_png,
    plot = combined_map,
    width = 14,
    height = 10,
    dpi = 600
  )
  
  ggsave(
    filename = output_tiff,
    plot = combined_map,
    width = 14,
    height = 10,
    dpi = 600,
    compression = "lzw"
  )
  
  cat("Saved PNG:", output_png, "\n")
  cat("Saved TIFF:", output_tiff, "\n")
  
  # ------------------------------------------------------------
  # Automatic summary note
  # ------------------------------------------------------------
  
  highest_change <- season_summary %>%
    filter(Mean_Percentage_Change == max(Mean_Percentage_Change, na.rm = TRUE)) %>%
    slice(1)
  
  lowest_change <- season_summary %>%
    filter(Mean_Percentage_Change == min(Mean_Percentage_Change, na.rm = TRUE)) %>%
    slice(1)
  
  general_mean <- mean(season_summary$Mean_Percentage_Change, na.rm = TRUE)
  
  tendency <- ifelse(
    general_mean >= 0,
    "an increasing tendency",
    "a decreasing tendency"
  )
  
  note_text <- paste0(
    season_name,
    " shows ",
    tendency,
    " across the evaluated scenarios and future periods. The highest basin-mean change occurs under ",
    highest_change$Scenario,
    " during ",
    highest_change$Period,
    " with a mean change of ",
    sprintf("%.2f", highest_change$Mean_Percentage_Change),
    "%. The lowest basin-mean change occurs under ",
    lowest_change$Scenario,
    " during ",
    lowest_change$Period,
    " with a mean change of ",
    sprintf("%.2f", lowest_change$Mean_Percentage_Change),
    "%. The basin boundary helps identify where projected rainfall changes occur inside the Upper Blue Nile Basin, while the NMA station dots represent the observation network used in this study, including ",
    station_names_text,
    "."
  )
  
  summary_notes_all <- rbind(
    summary_notes_all,
    data.frame(
      Season = season_name,
      Season_Code = season_code,
      Summary_Note = note_text
    )
  )
}

# ============================================================
# 10. SAVE SUMMARY TABLES AND NOTES
# ============================================================

summary_csv <- file.path(
  output_folder,
  "Spatial_Rainfall_Change_Statistics_All_Seasons_Stations.csv"
)

write.csv(
  mean_summary_all,
  summary_csv,
  row.names = FALSE
)

notes_csv <- file.path(
  output_folder,
  "Spatial_Rainfall_Change_Summary_Notes_Stations.csv"
)

write.csv(
  summary_notes_all,
  notes_csv,
  row.names = FALSE
)

notes_txt <- file.path(
  output_folder,
  "Spatial_Rainfall_Change_Summary_Notes_Stations.txt"
)

writeLines(
  paste0(
    summary_notes_all$Season,
    "\n",
    summary_notes_all$Summary_Note,
    "\n"
  ),
  con = notes_txt
)

# ============================================================
# 11. FINAL OUTPUT MESSAGE
# ============================================================

cat("\n================================================\n")
cat("FINAL SPATIAL RAINFALL CHANGE STATISTICS\n")
cat("================================================\n")
print(mean_summary_all)

cat("\n================================================\n")
cat("SUMMARY NOTES FOR EACH OUTPUT\n")
cat("================================================\n")
print(summary_notes_all)

cat("\nDONE SUCCESSFULLY!\n")
cat("High-quality maps with basin boundary and station dots saved in:\n")
cat(output_folder, "\n")
cat("\nStation table saved as:\n")
cat(file.path(output_folder, "NMA_Rainfall_Stations_Used_in_Study.csv"), "\n")
cat("\nStatistics table:\n")
cat(summary_csv, "\n")
cat("\nSummary notes:\n")
cat(notes_txt, "\n")

