# ============================================================
# CPI + JOC CLASSIFICATION USING GRID-LEVEL RIM RESULTS
# Input: rim_grid_results.csv
# Output: CPI ranking, JOC classes, bar plot, summary text
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

packages <- c("dplyr", "tidyr", "ggplot2", "classInt", "viridis")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# 2. INPUT AND OUTPUT PATHS
# ============================================================

input_csv <- "D:/Desktop/HWRM_Thesis/From Gh/RIM_Grid_Results/rim_grid_results.csv"

output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/CPI_JOC_Result"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

if (!file.exists(input_csv)) {
  stop("ERROR: rim_grid_results.csv was not found. Check input_csv path.")
}

# ============================================================
# 3. READ RIM GRID RESULTS
# ============================================================

rim_data <- read.csv(input_csv)

cat("Input data preview:\n")
print(head(rim_data))

cat("\nColumn names:\n")
print(names(rim_data))

# ============================================================
# 4. CHECK REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "grid_id",
  "GCM",
  "lmg",
  "first",
  "last",
  "pratt",
  "betasq",
  "car",
  "genizi"
)

missing_cols <- setdiff(required_cols, names(rim_data))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "ERROR: Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# 5. DEFINE RIM METRICS
# ============================================================

rim_metrics <- c(
  "lmg",
  "first",
  "last",
  "pratt",
  "betasq",
  "car",
  "genizi"
)

# All RIMs are treated as benefit criteria:
# higher RIM value = closer to ideal performance
benefit_criteria <- rim_metrics

# ============================================================
# 6. CONVERT TO LONG FORMAT
# ============================================================

rim_long <- rim_data %>%
  pivot_longer(
    cols = all_of(rim_metrics),
    names_to = "RIM",
    values_to = "RIM_value"
  )

# ============================================================
# 7. NORMALIZE RIM VALUES AT EACH GRID POINT
# ============================================================
# Benefit criterion:
# normalized = (x - min) / (max - min)
#
# After normalization:
# 1 = ideal / best
# 0 = worst
# ============================================================

rim_norm <- rim_long %>%
  group_by(grid_id, RIM) %>%
  mutate(
    min_value = min(RIM_value, na.rm = TRUE),
    max_value = max(RIM_value, na.rm = TRUE),
    
    normalized_value = ifelse(
      max_value == min_value,
      1,
      (RIM_value - min_value) / (max_value - min_value)
    )
  ) %>%
  ungroup()

# Remove impossible or missing normalized values
rim_norm <- rim_norm %>%
  filter(!is.na(normalized_value))

# ============================================================
# 8. ASSIGN EQUAL WEIGHTS TO RIM METRICS
# ============================================================

rim_weights <- data.frame(
  RIM = rim_metrics,
  weight = rep(1 / length(rim_metrics), length(rim_metrics))
)

cat("\nRIM weights used:\n")
print(rim_weights)

rim_norm <- rim_norm %>%
  left_join(rim_weights, by = "RIM")

# ============================================================
# 9. CALCULATE CPI AT EACH GRID POINT
# ============================================================
# CPI = [sum(wj * |1 - rij|^p)]^(1/p)
#
# rij = normalized RIM value
# wj  = weight
# p   = distance parameter
#
# Lower CPI = closer to ideal = better GCM
# ============================================================

p_value <- 2

cpi_grid <- rim_norm %>%
  group_by(grid_id, GCM) %>%
  summarise(
    CPI = sum(
      weight * abs(1 - normalized_value)^p_value,
      na.rm = TRUE
    )^(1 / p_value),
    .groups = "drop"
  )

# ============================================================
# 10. RANK GCMs AT EACH GRID POINT
# ============================================================

cpi_grid <- cpi_grid %>%
  group_by(grid_id) %>%
  arrange(CPI, .by_group = TRUE) %>%
  mutate(
    grid_rank = row_number()
  ) %>%
  ungroup()

cat("\nGrid-level CPI preview:\n")
print(head(cpi_grid))

# ============================================================
# 11. AGGREGATE CPI ACROSS ALL GRID POINTS
# ============================================================
# Mean CPI is used as basin-level CPI.
# Lower mean_CPI = better model.
# ============================================================

cpi_final <- cpi_grid %>%
  group_by(GCM) %>%
  summarise(
    mean_CPI = mean(CPI, na.rm = TRUE),
    median_CPI = median(CPI, na.rm = TRUE),
    sd_CPI = sd(CPI, na.rm = TRUE),
    min_CPI = min(CPI, na.rm = TRUE),
    max_CPI = max(CPI, na.rm = TRUE),
    mean_rank = mean(grid_rank, na.rm = TRUE),
    number_of_grids = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_CPI) %>%
  mutate(
    final_rank = row_number()
  )

cat("\nAggregated CPI ranking:\n")
print(cpi_final)

# ============================================================
# 12. JOC / JENKS OPTIMAL CLASSIFICATION
# ============================================================
# Class I = lowest CPI = most suitable
# Class V = highest CPI = least suitable
# ============================================================

number_of_classes <- 5

if (nrow(cpi_final) < number_of_classes) {
  number_of_classes <- nrow(cpi_final)
}

jenks_result <- classIntervals(
  cpi_final$mean_CPI,
  n = number_of_classes,
  style = "jenks"
)

jenks_breaks <- unique(jenks_result$brks)

# Safety check for duplicate breaks
if (length(jenks_breaks) <= 2) {
  number_of_classes <- 3
  
  jenks_result <- classIntervals(
    cpi_final$mean_CPI,
    n = number_of_classes,
    style = "jenks"
  )
  
  jenks_breaks <- unique(jenks_result$brks)
}

cat("\nJOC / Jenks breaks:\n")
print(jenks_breaks)

cpi_final$Class_number <- cut(
  cpi_final$mean_CPI,
  breaks = jenks_breaks,
  include.lowest = TRUE,
  labels = FALSE
)

cpi_final <- cpi_final %>%
  mutate(
    JOC_Class = case_when(
      Class_number == 1 ~ "I",
      Class_number == 2 ~ "II",
      Class_number == 3 ~ "III",
      Class_number == 4 ~ "IV",
      Class_number == 5 ~ "V",
      TRUE ~ "Unclassified"
    ),
    
    Suitability = case_when(
      JOC_Class == "I" ~ "Most suitable",
      JOC_Class == "II" ~ "Suitable",
      JOC_Class == "III" ~ "Moderate",
      JOC_Class == "IV" ~ "Low suitability",
      JOC_Class == "V" ~ "Least suitable",
      TRUE ~ "Unclassified"
    )
  )

cat("\nFinal CPI + JOC classification:\n")
print(cpi_final)

# ============================================================
# 13. SAVE OUTPUT TABLES
# ============================================================

write.csv(
  rim_norm,
  file.path(output_folder, "normalized_RIM_values_by_grid.csv"),
  row.names = FALSE
)

write.csv(
  cpi_grid,
  file.path(output_folder, "CPI_by_grid_and_GCM.csv"),
  row.names = FALSE
)

write.csv(
  cpi_final,
  file.path(output_folder, "Final_CPI_JOC_GCM_ranking.csv"),
  row.names = FALSE
)

# ============================================================
# 14. PLOT CPI BAR CHART LIKE REFERENCE FIGURE
# ============================================================

plot_data <- cpi_final %>%
  arrange(desc(mean_CPI))

plot_data$GCM <- factor(
  plot_data$GCM,
  levels = plot_data$GCM
)

p1 <- ggplot(
  plot_data,
  aes(x = mean_CPI, y = GCM, fill = JOC_Class)
) +
  geom_col(width = 0.75, color = "grey30") +
  
  geom_text(
    aes(label = sprintf("%.5f", mean_CPI)),
    hjust = -0.1,
    size = 3.4
  ) +
  
  scale_fill_manual(
    values = c(
      "I" = "#66c2a5",
      "II" = "#d9ef8b",
      "III" = "#bdbdbd",
      "IV" = "#fc8d59",
      "V" = "#80b1d3",
      "Unclassified" = "grey70"
    ),
    name = "Class"
  ) +
  
  labs(
    title = "GCM CPI Scores by JOC Class",
    subtitle = "Evaluation of GCMs using Compromise Programming Index",
    x = "Compromise Programming Index (CPI)",
    y = "General Circulation Models (GCMs)"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  expand_limits(
    x = max(plot_data$mean_CPI, na.rm = TRUE) * 1.15
  )

print(p1)

ggsave(
  filename = file.path(output_folder, "GCM_CPI_scores_by_JOC_class.png"),
  plot = p1,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 15. PLOT GRID-LEVEL CPI HEATMAP
# ============================================================

p2 <- ggplot(
  cpi_grid,
  aes(x = GCM, y = factor(grid_id), fill = CPI)
) +
  geom_tile(color = "white") +
  
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "CPI"
  ) +
  
  labs(
    title = "Grid-Level CPI Heatmap",
    subtitle = "Lower CPI indicates better model performance at each grid point",
    x = "GCM",
    y = "Grid ID"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(face = "bold")
  )

print(p2)

ggsave(
  filename = file.path(output_folder, "Grid_level_CPI_heatmap.png"),
  plot = p2,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================
# 16. CREATE THESIS-STYLE SUMMARY TEXT
# ============================================================

best_model <- cpi_final %>%
  filter(final_rank == 1)

class_I_models <- cpi_final %>%
  filter(JOC_Class == "I") %>%
  arrange(mean_CPI)

class_II_models <- cpi_final %>%
  filter(JOC_Class == "II") %>%
  arrange(mean_CPI)

summary_text <- paste0(
  "The RIM results across each grid point were ranked and aggregated using the ",
  "Compromise Programming Index (CPI), which evaluates the overall proximity of each ",
  "GCM to the ideal performance. The CPI values were subsequently categorized using ",
  "Jenks Optimal Classification (JOC) to identify the most suitable GCMs for future ",
  "projection purposes. The results indicate that ",
  best_model$GCM[1],
  " achieved the lowest CPI value of ",
  sprintf("%.5f", best_model$mean_CPI[1]),
  ", indicating the best overall performance across the evaluated grid cells."
)

writeLines(
  summary_text,
  file.path(output_folder, "CPI_JOC_summary_text.txt")
)

cat("\nSummary text:\n")
cat(summary_text, "\n")

# ============================================================
# 17. FINAL MESSAGE
# ============================================================

cat("\nDONE SUCCESSFULLY!\n")
cat("CPI and JOC results saved in:\n")
cat(output_folder, "\n")
