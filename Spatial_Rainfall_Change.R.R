# ============================================================
# SPATIAL CHANGE IN ANNUAL AND ETHIOPIAN SEASONAL RAINFALL
# Study area: Upper Blue Nile Basin, Ethiopia
# Historical baseline: XGBoost-corrected CHIRPS, 1994–2014
# Future rainfall: Final BCCA-PT downscaled MIROC6
# Future periods: 2030–2050 and 2051–2080
# Scenarios: SSP2-4.5 and SSP5-8.5
#
# Ethiopian rainfall/agricultural seasons:
# Annual = Jan–Dec
# Belg   = Mar–May
# Kiremt = Jun–Aug
# Meher  = Sep–Nov
# Bega   = Dec–Feb
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c(
  "terra", "dplyr", "ggplot2", "tidyr",
  "patchwork", "stringr", "scales"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT FOLDERS
# ============================================================

chirps_folder <- "D:/Desktop/HWRM_Thesis/From Gh/XGBOOST_corrected"

future_downscaled_folder <- "D:/Desktop/HWRM_Thesis/From Gh/BCCA_PT_Downscaled_MIROC6_FINAL/03_Future_Downscaled"

basin_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Gebres/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/abay"

output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/Spatial_Rainfall_Change_FINAL_Ethiopian_Seasons"

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "GeoTIFF"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Figures", "Individual_Maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Figures", "Combined_Seasonal_Maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_folder, "Tables"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. BASIC SETTINGS
# ============================================================

selected_gcm <- "MIROC6"

ssp_scenarios <- c("ssp245", "ssp585")

historical_start <- 1994
historical_end <- 2014

future_periods <- data.frame(
  period = c("Near_future_2030_2050", "Far_future_2051_2080"),
  start_year = c(2030, 2051),
  end_year = c(2050, 2080)
)

seasons <- list(
  Annual = 1:12,
  Belg_MAM = c(3, 4, 5),
  Kiremt_JJA = c(6, 7, 8),
  Meher_SON = c(9, 10, 11),
  Bega_DJF = c(12, 1, 2)
)

season_labels <- c(
  Annual = "Annual rainfall",
  Belg_MAM = "Belg rainfall (March–May)",
  Kiremt_JJA = "Kiremt rainfall (June–August)",
  Meher_SON = "Meher rainfall (September–November)",
  Bega_DJF = "Bega rainfall (December–February)"
)

season_labels_short <- c(
  Annual = "Annual rainfall",
  Belg_MAM = "Belg rainfall (Mar–May)",
  Kiremt_JJA = "Kiremt rainfall (Jun–Aug)",
  Meher_SON = "Meher rainfall (Sep–Nov)",
  Bega_DJF = "Bega rainfall (Dec–Feb)"
)

# KEEPING THE PREVIOUS COLOR COMBINATIONS
season_palettes <- list(
  Annual = c(low = "#8c510a", mid = "#f7f7f7", high = "#1b7837"),
  Belg_MAM = c(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac"),
  Kiremt_JJA = c(low = "#a6611a", mid = "#ffffbf", high = "#018571"),
  Meher_SON = c(low = "#b35806", mid = "#f7f7f7", high = "#542788"),
  Bega_DJF = c(low = "#762a83", mid = "#f7f7f7", high = "#1b7837")
)

# ============================================================
# 4. FUNCTIONS
# ============================================================

extract_year_month <- function(file_path) {
  
  fname <- basename(file_path)
  
  match1 <- regmatches(
    fname,
    regexpr("(19[0-9]{2}|20[0-9]{2})[_\\.-](0[1-9]|1[0-2])", fname)
  )
  
  if (length(match1) > 0 && match1 != "") {
    return(gsub("[_\\.-]", "-", match1))
  }
  
  match2 <- regmatches(
    fname,
    regexpr("(19[0-9]{2}|20[0-9]{2})(0[1-9]|1[0-2])", fname)
  )
  
  if (length(match2) > 0 && match2 != "") {
    yr <- substr(match2, 1, 4)
    mn <- substr(match2, 5, 6)
    return(paste0(yr, "-", mn))
  }
  
  return(NA)
}

make_inventory <- function(folder, dataset_name) {
  
  files <- list.files(
    path = folder,
    pattern = "\\.(tif|tiff)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  files <- files[
    !grepl("BACKUP", files, ignore.case = TRUE) &
      !grepl("BAD", files, ignore.case = TRUE) &
      !grepl("Tables", files, ignore.case = TRUE) &
      !grepl("GeoTIFF", files, ignore.case = TRUE)
  ]
  
  if (length(files) == 0) {
    stop(paste("ERROR: No valid tif files found in:", folder))
  }
  
  inv <- data.frame(
    date_id = sapply(files, extract_year_month),
    file = files,
    dataset = dataset_name
  )
  
  inv <- inv[!is.na(inv$date_id), ]
  inv$year <- as.numeric(substr(inv$date_id, 1, 4))
  inv$month <- as.numeric(substr(inv$date_id, 6, 7))
  inv <- inv[order(inv$date_id), ]
  
  return(inv)
}

safe_percent_change_raster <- function(future_raster, baseline_raster) {
  
  pct <- ((future_raster - baseline_raster) / baseline_raster) * 100
  
  pct <- ifel(
    baseline_raster <= 0 | is.na(baseline_raster),
    NA,
    pct
  )
  
  return(pct)
}

calculate_period_mean_total <- function(stack, dates_table, months_selected,
                                        start_year, end_year, season_name) {
  
  dates_table$layer_index <- seq_len(nrow(dates_table))
  
  use_table <- dates_table[
    dates_table$year >= start_year &
      dates_table$year <= end_year &
      dates_table$month %in% months_selected,
  ]
  
  if (nrow(use_table) == 0) {
    return(NULL)
  }
  
  seasonal_list <- list()
  count <- 1
  
  if (season_name == "Bega_DJF") {
    
    dates_table$season_year <- ifelse(
      dates_table$month == 12,
      dates_table$year + 1,
      dates_table$year
    )
    
    season_years <- sort(unique(dates_table$season_year))
    
    for (syear in season_years) {
      
      idx <- dates_table$layer_index[
        dates_table$season_year == syear &
          dates_table$month %in% c(12, 1, 2)
      ]
      
      if (length(idx) != 3) next
      if (syear < start_year | syear > end_year) next
      
      total_raster <- sum(stack[[idx]], na.rm = TRUE)
      seasonal_list[[count]] <- total_raster
      count <- count + 1
    }
    
  } else {
    
    years <- sort(unique(use_table$year))
    
    for (yr in years) {
      
      idx <- use_table$layer_index[
        use_table$year == yr &
          use_table$month %in% months_selected
      ]
      
      if (length(idx) == 0) next
      
      total_raster <- sum(stack[[idx]], na.rm = TRUE)
      seasonal_list[[count]] <- total_raster
      count <- count + 1
    }
  }
  
  if (length(seasonal_list) == 0) {
    return(NULL)
  }
  
  period_stack <- rast(seasonal_list)
  mean_total <- mean(period_stack, na.rm = TRUE)
  
  return(mean_total)
}

scenario_label <- function(ssp) {
  ifelse(ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")
}

period_label <- function(period) {
  ifelse(
    period == "Near_future_2030_2050",
    "Near future",
    "Far future"
  )
}

period_label_full <- function(period) {
  ifelse(
    period == "Near_future_2030_2050",
    "Near future (2030–2050)",
    "Far future (2051–2080)"
  )
}

save_change_map <- function(r, title_text, subtitle_text, output_png,
                            season_name, limit_abs,
                            basin_df = NULL,
                            stations_df = NULL) {
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  
  if (ncol(df) < 3) {
    stop("Raster dataframe has no value column.")
  }
  
  names(df)[3] <- "change"
  
  pal <- season_palettes[[season_name]]
  
  p <- ggplot(df, aes(x = x, y = y, fill = change)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_gradient2(
      low = pal["low"],
      mid = pal["mid"],
      high = pal["high"],
      midpoint = 0,
      limits = c(-limit_abs, limit_abs),
      oob = scales::squish,
      name = "Change (%)"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10,
        margin = margin(b = 5)
      ),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 8),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.position = "right",
      plot.margin = margin(6, 6, 6, 6)
    )
  
  # Basin outline
  if (!is.null(basin_df) && all(c("x", "y", "geom", "part") %in% names(basin_df))) {
    p <- p +
      geom_path(
        data = basin_df,
        aes(x = x, y = y, group = interaction(geom, part)),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      )
  }
  
  # NMA station dots only
  if (!is.null(stations_df)) {
    p <- p +
      geom_point(
        data = stations_df,
        aes(x = lon, y = lat),
        inherit.aes = FALSE,
        shape = 21,
        size = 2.4,
        stroke = 0.7,
        fill = "yellow",
        color = "black"
      )
  }
  
  ggsave(
    filename = output_png,
    plot = p,
    width = 8.5,
    height = 6.2,
    dpi = 400
  )
  
  return(p)
}

create_season_sentence <- function(summary_table, season_key) {
  
  df <- summary_table %>%
    filter(season == season_key) %>%
    arrange(mean_percentage_change)
  
  if (nrow(df) == 0) return("")
  
  label <- unique(df$season_label)[1]
  
  lowest <- df[1, ]
  highest <- df[nrow(df), ]
  
  direction <- ifelse(
    all(df$mean_percentage_change > 0, na.rm = TRUE),
    "shows an increasing tendency",
    ifelse(
      all(df$mean_percentage_change < 0, na.rm = TRUE),
      "shows a decreasing tendency",
      "shows mixed spatial changes"
    )
  )
  
  paste0(
    label, " ", direction, " across the evaluated scenarios and future periods. ",
    "The highest basin-mean change occurs under ", highest$SSP_label,
    " during ", highest$period_label_full, ", with a mean change of ",
    round(highest$mean_percentage_change, 2), "%. ",
    "The lowest basin-mean change occurs under ", lowest$SSP_label,
    " during ", lowest$period_label_full, ", with a mean change of ",
    round(lowest$mean_percentage_change, 2), "%. ",
    "Spatial values range from ",
    round(min(df$min_percentage_change, na.rm = TRUE), 2), "% to ",
    round(max(df$max_percentage_change, na.rm = TRUE), 2),
    "% across the basin.\n\n"
  )
}

# ============================================================
# 5. READ INVENTORIES
# ============================================================

chirps_inv <- make_inventory(chirps_folder, "XGBoost_CHIRPS")

chirps_inv <- chirps_inv[
  chirps_inv$year >= historical_start &
    chirps_inv$year <= historical_end,
]

cat("Historical CHIRPS files:", nrow(chirps_inv), "\n")

if (nrow(chirps_inv) != 252) {
  warning("Expected 252 CHIRPS monthly files for 1994–2014.")
}

future_inv_all <- make_inventory(
  future_downscaled_folder,
  "MIROC6_BCCA_PT_future"
)

future_inv_all <- future_inv_all[
  grepl(selected_gcm, future_inv_all$file, ignore.case = TRUE) &
    grepl("BCCA_PT", future_inv_all$file, ignore.case = TRUE),
]

cat("Future BCCA-PT MIROC6 files:", nrow(future_inv_all), "\n")
cat("SSP245 files:", sum(grepl("ssp245", future_inv_all$file, ignore.case = TRUE)), "\n")
cat("SSP585 files:", sum(grepl("ssp585", future_inv_all$file, ignore.case = TRUE)), "\n")

if (nrow(future_inv_all) != 1224) {
  warning("Expected 1224 future files: 612 for ssp245 and 612 for ssp585.")
}

write.csv(
  chirps_inv,
  file.path(output_folder, "Tables", "CHIRPS_XGBoost_inventory_1994_2014.csv"),
  row.names = FALSE
)

write.csv(
  future_inv_all,
  file.path(output_folder, "Tables", "Future_BCCA_PT_MIROC6_inventory.csv"),
  row.names = FALSE
)

# ============================================================
# 6. READ HISTORICAL CHIRPS STACK
# ============================================================

cat("\nReading historical XGBoost-corrected CHIRPS stack...\n")

chirps_stack <- rast(chirps_inv$file)
ref_raster <- chirps_stack[[1]]

cat("Reference CHIRPS grid:\n")
print(ref_raster)

# ============================================================
# 7. LOAD BASIN BOUNDARY
# ============================================================

basin_shp <- list.files(
  basin_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)[1]

if (is.na(basin_shp)) {
  warning("No basin shapefile found. Maps will be created without basin outline.")
  basin_vect <- NULL
  basin_df <- NULL
} else {
  basin_vect <- vect(basin_shp)
  basin_vect <- project(basin_vect, crs(ref_raster))
  basin_geom <- as.data.frame(geom(basin_vect))
  
  if (all(c("x", "y") %in% names(basin_geom))) {
    basin_df <- basin_geom
  } else {
    basin_df <- NULL
  }
}

# ============================================================
# 7B. NMA STATION COORDINATES
# Dots only, no labels
# Replace coordinates if you have official NMA coordinates
# ============================================================

stations_df <- data.frame(
  station = c(
    "Assosa",
    "Bahir Dar",
    "Debre Birhan",
    "Debre Markos",
    "Debre Tabor",
    "Dessie",
    "Gondar",
    "Nekemte",
    "Finote Selam"
  ),
  lon = c(
    34.53,
    37.39,
    39.53,
    37.73,
    38.02,
    39.63,
    37.47,
    36.55,
    37.27
  ),
  lat = c(
    10.07,
    11.60,
    9.68,
    10.34,
    11.85,
    11.13,
    12.60,
    9.08,
    10.70
  )
)

write.csv(
  stations_df,
  file.path(output_folder, "Tables", "NMA_station_coordinates_used_for_maps.csv"),
  row.names = FALSE
)

# ============================================================
# 8. HISTORICAL BASELINE RASTERS
# ============================================================

cat("\nCalculating historical baseline rasters...\n")

hist_dates <- chirps_inv[, c("date_id", "year", "month")]

baseline_rasters <- list()

for (season_name in names(seasons)) {
  
  months_selected <- seasons[[season_name]]
  
  baseline_rasters[[season_name]] <- calculate_period_mean_total(
    stack = chirps_stack,
    dates_table = hist_dates,
    months_selected = months_selected,
    start_year = historical_start,
    end_year = historical_end,
    season_name = season_name
  )
  
  if (is.null(baseline_rasters[[season_name]])) {
    warning(paste("Baseline raster is NULL for", season_name))
    next
  }
  
  out_baseline <- file.path(
    output_folder,
    "GeoTIFF",
    paste0("Baseline_", season_name, "_1994_2014.tif")
  )
  
  writeRaster(
    baseline_rasters[[season_name]],
    out_baseline,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
  )
}

# ============================================================
# 9. PROCESS FUTURE SCENARIOS, PERIODS, AND SEASONS
# ============================================================

summary_table <- data.frame()
map_records <- data.frame()

for (ssp in ssp_scenarios) {
  
  ssp_inv <- future_inv_all[
    grepl(ssp, future_inv_all$file, ignore.case = TRUE),
  ]
  
  if (nrow(ssp_inv) == 0) {
    warning(paste("No future files found for", ssp))
    next
  }
  
  ssp_inv <- ssp_inv[order(ssp_inv$date_id), ]
  
  cat("\nProcessing scenario:", ssp, "\n")
  
  future_stack <- rast(ssp_inv$file)
  
  if (!compareGeom(ref_raster, future_stack[[1]], stopOnError = FALSE)) {
    cat("Resampling future", ssp, "to CHIRPS grid...\n")
    future_stack <- resample(future_stack, ref_raster, method = "bilinear")
  }
  
  for (p_i in 1:nrow(future_periods)) {
    
    period_name <- future_periods$period[p_i]
    sy <- future_periods$start_year[p_i]
    ey <- future_periods$end_year[p_i]
    
    cat("Period:", period_name, "\n")
    
    for (season_name in names(seasons)) {
      
      months_selected <- seasons[[season_name]]
      
      future_mean_raster <- calculate_period_mean_total(
        stack = future_stack,
        dates_table = ssp_inv[, c("date_id", "year", "month")],
        months_selected = months_selected,
        start_year = sy,
        end_year = ey,
        season_name = season_name
      )
      
      if (is.null(future_mean_raster)) next
      
      baseline_raster <- baseline_rasters[[season_name]]
      
      pct_change <- safe_percent_change_raster(
        future_raster = future_mean_raster,
        baseline_raster = baseline_raster
      )
      
      out_name <- paste0(
        selected_gcm, "_", ssp, "_", period_name, "_",
        season_name, "_percentage_change"
      )
      
      future_tif <- file.path(
        output_folder,
        "GeoTIFF",
        paste0(out_name, "_future_mean.tif")
      )
      
      pct_tif <- file.path(
        output_folder,
        "GeoTIFF",
        paste0(out_name, ".tif")
      )
      
      writeRaster(
        future_mean_raster,
        future_tif,
        overwrite = TRUE,
        wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
      )
      
      writeRaster(
        pct_change,
        pct_tif,
        overwrite = TRUE,
        wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
      )
      
      pct_values <- values(pct_change, mat = FALSE)
      
      mean_change <- mean(pct_values, na.rm = TRUE)
      min_change <- min(pct_values, na.rm = TRUE)
      max_change <- max(pct_values, na.rm = TRUE)
      sd_change <- sd(pct_values, na.rm = TRUE)
      
      summary_table <- rbind(
        summary_table,
        data.frame(
          GCM = selected_gcm,
          SSP = ssp,
          SSP_label = scenario_label(ssp),
          period = period_name,
          period_label = period_label(period_name),
          period_label_full = period_label_full(period_name),
          season = season_name,
          season_label = season_labels[season_name],
          mean_percentage_change = mean_change,
          min_percentage_change = min_change,
          max_percentage_change = max_change,
          sd_percentage_change = sd_change,
          pct_raster = pct_tif,
          future_mean_raster = future_tif
        )
      )
      
      map_records <- rbind(
        map_records,
        data.frame(
          GCM = selected_gcm,
          SSP = ssp,
          SSP_label = scenario_label(ssp),
          period = period_name,
          period_label = period_label(period_name),
          period_label_full = period_label_full(period_name),
          season = season_name,
          season_label = season_labels[season_name],
          mean_percentage_change = mean_change,
          pct_raster = pct_tif
        )
      )
    }
  }
}

write.csv(
  summary_table,
  file.path(output_folder, "Tables", "spatial_rainfall_change_summary.csv"),
  row.names = FALSE
)

cat("\nSpatial rainfall change summary:\n")
print(summary_table)

# ============================================================
# 10. CREATE INDIVIDUAL AND COMBINED MAPS
# Main title on top + short title/subtitle in each map
# + station dots
# ============================================================

for (season_name in names(seasons)) {
  
  season_maps <- map_records %>%
    filter(season == season_name)
  
  if (nrow(season_maps) == 0) next
  
  all_vals <- c()
  
  for (i in 1:nrow(season_maps)) {
    rr <- rast(season_maps$pct_raster[i])
    all_vals <- c(all_vals, values(rr, mat = FALSE))
  }
  
  all_vals <- all_vals[is.finite(all_vals)]
  
  if (length(all_vals) == 0) next
  
  limit_abs <- ceiling(max(abs(all_vals), na.rm = TRUE) / 5) * 5
  
  if (!is.finite(limit_abs) || limit_abs == 0) {
    limit_abs <- 10
  }
  
  plot_list <- list()
  
  for (i in 1:nrow(season_maps)) {
    
    rr <- rast(season_maps$pct_raster[i])
    
    short_title <- paste0(
      season_maps$SSP_label[i],
      " | ",
      season_maps$period_label[i]
    )
    
    short_subtitle <- paste0(
      "Mean = ",
      round(season_maps$mean_percentage_change[i], 2),
      "%"
    )
    
    png_file <- file.path(
      output_folder,
      "Figures",
      "Individual_Maps",
      paste0(
        selected_gcm, "_",
        season_maps$SSP[i], "_",
        season_maps$period[i], "_",
        season_name,
        "_spatial_change.png"
      )
    )
    
    p <- save_change_map(
      r = rr,
      title_text = short_title,
      subtitle_text = short_subtitle,
      output_png = png_file,
      season_name = season_name,
      limit_abs = limit_abs,
      basin_df = basin_df,
      stations_df = stations_df
    )
    
    plot_list[[i]] <- p
  }
  
  combined_title <- paste0(
    "Projected spatial change in ",
    season_labels_short[season_name],
    " over the Upper Blue Nile Basin"
  )
  
  combined_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = combined_title,
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 18,
          hjust = 0.5,
          lineheight = 0.95,
          margin = margin(b = 10)
        ),
        plot.margin = margin(14, 14, 14, 14)
      )
    )
  
  combined_png <- file.path(
    output_folder,
    "Figures",
    "Combined_Seasonal_Maps",
    paste0("Combined_", season_name, "_Spatial_Rainfall_Change.png")
  )
  
  ggsave(
    filename = combined_png,
    plot = combined_plot,
    width = 15,
    height = 11,
    dpi = 400
  )
  
  cat("Combined map saved:", combined_png, "\n")
}

# ============================================================
# 11. SUMMARY BARPLOT FOR ALL SEASONS
# ============================================================

summary_plot_data <- summary_table %>%
  mutate(
    season_label = factor(
      season_label,
      levels = season_labels
    )
  )

p_summary <- ggplot(
  summary_plot_data,
  aes(x = period_label_full, y = mean_percentage_change, fill = SSP_label)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    width = 0.7
  ) +
  facet_wrap(~ season_label, scales = "free_y", ncol = 2) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  labs(
    title = "Mean spatial rainfall change over the Upper Blue Nile Basin",
    subtitle = "Final BCCA-PT downscaled MIROC6 projections relative to 1994–2014 XGBoost-corrected CHIRPS baseline",
    x = "Future period",
    y = "Mean rainfall change (%)",
    fill = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(14, 14, 14, 14)
  )

ggsave(
  filename = file.path(output_folder, "Figures", "Mean_Spatial_Rainfall_Change_All_Seasons_Barplot.png"),
  plot = p_summary,
  width = 12,
  height = 8,
  dpi = 400
)

# ============================================================
# 12. CREATE NUMERICAL SUMMARY NOTE
# ============================================================

summary_note <- paste0(
  "Spatial rainfall change summary for the Upper Blue Nile Basin\n\n",
  "Spatial changes in annual and Ethiopian seasonal rainfall were assessed by comparing final BCCA-PT downscaled MIROC6 future rainfall projections under SSP2-4.5 and SSP5-8.5 with the 1994–2014 XGBoost-corrected CHIRPS baseline. Percentage change maps were produced for the near-future period, 2030–2050, and the far-future period, 2051–2080. Seasonal changes were summarized using Ethiopia’s rainfall and agricultural seasons: Belg (March–May), Kiremt (June–August), Meher (September–November), and Bega (December–February). Annual rainfall was calculated using all months.\n\n",
  create_season_sentence(summary_table, "Annual"),
  create_season_sentence(summary_table, "Belg_MAM"),
  create_season_sentence(summary_table, "Kiremt_JJA"),
  create_season_sentence(summary_table, "Meher_SON"),
  create_season_sentence(summary_table, "Bega_DJF"),
  "The basin boundary and the nine NMA station dots were used to show the spatial distribution of projected rainfall change inside the Upper Blue Nile Basin."
)

writeLines(
  summary_note,
  file.path(output_folder, "Tables", "spatial_rainfall_change_numerical_summary_note.txt")
)

cat("\nNumerical summary note:\n")
cat(summary_note)

# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\nDONE SUCCESSFULLY!\n")
cat("Spatial rainfall change outputs saved in:\n")
cat(output_folder, "\n\n")

cat("Summary table:\n")
cat(file.path(output_folder, "Tables", "spatial_rainfall_change_summary.csv"), "\n\n")

cat("Numerical summary note:\n")
cat(file.path(output_folder, "Tables", "spatial_rainfall_change_numerical_summary_note.txt"), "\n\n")

cat("Station coordinate table:\n")
cat(file.path(output_folder, "Tables", "NMA_station_coordinates_used_for_maps.csv"), "\n\n")

cat("Combined maps folder:\n")
cat(file.path(output_folder, "Figures", "Combined_Seasonal_Maps"), "\n\n")

cat("Individual maps folder:\n")
cat(file.path(output_folder, "Figures", "Individual_Maps"), "\n")

