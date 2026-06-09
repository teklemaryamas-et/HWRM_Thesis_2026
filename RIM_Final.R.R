# =========================================================
# COMPLETE FINAL WORKABLE CODE
# GCM RELATIVE IMPORTANCE + MCDA + PARETO ANALYSIS
# FULLY FIXED VERSION
# =========================================================

# =========================================================
# 1. INSTALL PACKAGES (RUN ONCE)
# =========================================================
install.packages(c(
  "terra",
  "relaimpo",
  "ggplot2",
  "dplyr",
  "reshape2",
  "viridis",
  "patchwork",
  "pheatmap"
))

# =========================================================
# 2. LOAD LIBRARIES
# =========================================================
library(terra)
library(relaimpo)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(patchwork)
library(pheatmap)

# =========================================================
# 3. PATHS
# =========================================================

# OBSERVED DATA FOLDER
obs_path <- "D:/Desktop/HWRM_Thesis/From Gh/TM proposal/data/CHIRPS_Monthly_Abay/XGBOOST_corrected/resampled_0_25deg"

# GCM ROOT FOLDER
gcm_root <- "D:/Desktop/HWRM_Thesis/From Gh/NEX-GDDP-CMIP6"

# =========================================================
# 4. READ OBSERVED DATA FUNCTION
# =========================================================
read_obs_series <- function(folder){
  
  files <- list.files(
    folder,
    pattern="\\.tif$",
    recursive=TRUE,
    full.names=TRUE
  )
  
  cat("\nOBSERVED FILES FOUND:", length(files), "\n")
  
  if(length(files) == 0){
    stop("No observed tif files found.")
  }
  
  # Empty dataframe
  df <- data.frame(
    Date=as.Date(character()),
    OBS=numeric(),
    stringsAsFactors=FALSE
  )
  
  # Loop through files
  for(f in files){
    
    fname <- basename(f)
    
    # Extract date
    nums <- regmatches(
      fname,
      gregexpr("[0-9]{4}_?[0-9]{2}", fname)
    )[[1]]
    
    if(length(nums) == 0) next
    
    date_str <- gsub("_","", nums[1])
    
    # Convert to Date
    date <- tryCatch(
      as.Date(
        paste0(date_str,"01"),
        format="%Y%m%d"
      ),
      error=function(e) NA
    )
    
    if(is.na(date)) next
    
    # Read raster
    r <- tryCatch(
      rast(f),
      error=function(e) NULL
    )
    
    if(is.null(r)) next
    
    # Mean value
    val <- tryCatch(
      global(r,"mean",na.rm=TRUE)[1,1],
      error=function(e) NA
    )
    
    if(is.na(val)) next
    
    temp <- data.frame(
      Date=as.Date(date),
      OBS=as.numeric(val),
      stringsAsFactors=FALSE
    )
    
    df <- rbind(df,temp)
  }
  
  if(nrow(df) == 0){
    stop("No valid observed data extracted.")
  }
  
  # Fix types
  df$Date <- as.Date(df$Date)
  df$OBS <- as.numeric(df$OBS)
  
  # Remove NA
  df <- df[!is.na(df$Date), , drop=FALSE]
  
  # Filter years
  df <- df[
    df$Date >= as.Date("1990-01-01") &
      df$Date <= as.Date("2014-12-01"),
    ,
    drop=FALSE
  ]
  
  # Sort
  df <- df[order(df$Date), , drop=FALSE]
  
  rownames(df) <- NULL
  
  return(df)
}

# =========================================================
# 5. LOAD OBSERVED DATA
# =========================================================
cat("\n=================================\n")
cat("LOADING OBSERVED DATA\n")
cat("=================================\n")

obs <- read_obs_series(obs_path)

cat("\nOBSERVED DATA DIMENSIONS:\n")
print(dim(obs))

print(head(obs))

# =========================================================
# 6. LOAD GCM FILES
# =========================================================
cat("\n=================================\n")
cat("LOADING GCM FILES\n")
cat("=================================\n")

gcm_files <- list.files(
  gcm_root,
  pattern="\\.tif$",
  recursive=TRUE,
  full.names=TRUE
)

cat("\nTOTAL GCM FILES FOUND:", length(gcm_files), "\n")

if(length(gcm_files) == 0){
  
  stop("NO GCM TIFF FILES FOUND.")
}

# =========================================================
# 7. READ GCM DATA
# =========================================================
gcm_df <- data.frame(
  Date=as.Date(character()),
  Model=character(),
  Value=numeric(),
  stringsAsFactors=FALSE
)

for(f in gcm_files){
  
  fname <- basename(f)
  
  cat("Processing:", fname, "\n")
  
  # Extract date
  nums <- regmatches(
    fname,
    gregexpr("[0-9]{4}_?[0-9]{2}", fname)
  )[[1]]
  
  if(length(nums) == 0) next
  
  date_str <- gsub("_","", nums[1])
  
  date <- tryCatch(
    as.Date(
      paste0(date_str,"01"),
      format="%Y%m%d"
    ),
    error=function(e) NA
  )
  
  if(is.na(date)) next
  
  # Filter years
  if(date < as.Date("1990-01-01") |
     date > as.Date("2014-12-01")){
    next
  }
  
  # Extract model name
  parts <- unlist(strsplit(fname,"_"))
  
  model <- parts[1]
  
  model <- gsub("\\.tif","",model)
  
  # Read raster
  r <- tryCatch(
    rast(f),
    error=function(e) NULL
  )
  
  if(is.null(r)) next
  
  # Extract mean
  val <- tryCatch(
    global(r,"mean",na.rm=TRUE)[1,1],
    error=function(e) NA
  )
  
  if(is.na(val)) next
  
  temp <- data.frame(
    Date=as.Date(date),
    Model=as.character(model),
    Value=as.numeric(val),
    stringsAsFactors=FALSE
  )
  
  gcm_df <- rbind(gcm_df,temp)
}

# =========================================================
# 8. CHECK GCM DATA
# =========================================================
if(nrow(gcm_df) == 0){
  
  stop("NO VALID GCM DATA FOUND.")
}

cat("\nGCM DATA DIMENSIONS:\n")
print(dim(gcm_df))

print(head(gcm_df))

# =========================================================
# 9. CONVERT TO WIDE FORMAT
# =========================================================
gcm_wide <- dcast(
  gcm_df,
  Date ~ Model,
  value.var="Value"
)

gcm_wide$Date <- as.Date(gcm_wide$Date)

cat("\nGCM WIDE DIMENSIONS:\n")
print(dim(gcm_wide))

print(head(gcm_wide))

# =========================================================
# 10. MERGE OBSERVED + GCM
# =========================================================
cat("\n=================================\n")
cat("MERGING DATA\n")
cat("=================================\n")

df <- merge(
  obs,
  gcm_wide,
  by="Date",
  all=FALSE
)

df <- na.omit(df)

cat("\nFINAL DATA DIMENSIONS:\n")
print(dim(df))

print(names(df))

# =========================================================
# 11. CHECK PREDICTORS
# =========================================================
predictors <- setdiff(
  names(df),
  c("Date","OBS")
)

if(length(predictors) == 0){
  
  stop("NO GCM PREDICTORS FOUND.")
}

cat("\nORIGINAL PREDICTORS:\n")
print(predictors)

# =========================================================
# 12. FIX INVALID NAMES
# =========================================================
clean_names <- make.names(predictors)

colnames(df)[match(predictors,names(df))] <- clean_names

predictors <- clean_names

cat("\nCLEANED PREDICTOR NAMES:\n")
print(predictors)

# =========================================================
# 13. BUILD MODEL DATA
# =========================================================
model_df <- df[, c("OBS", predictors), drop=FALSE]

model_df <- na.omit(model_df)

# =========================================================
# 14. LINEAR MODEL
# =========================================================
cat("\n=================================\n")
cat("LINEAR MODEL\n")
cat("=================================\n")

formula_text <- paste(
  "OBS ~",
  paste(predictors, collapse=" + ")
)

cat("\nMODEL FORMULA:\n")
cat(formula_text,"\n")

formula_obj <- as.formula(formula_text)

lm_model <- lm(
  formula_obj,
  data=model_df
)

print(summary(lm_model))

# =========================================================
# 15. CHECK MODEL
# =========================================================
if(length(coef(lm_model)) <= 1){
  
  stop("MODEL CONTAINS ONLY INTERCEPT.")
}

# =========================================================
# 16. RELATIVE IMPORTANCE ANALYSIS
# =========================================================
cat("\n=================================\n")
cat("RELATIVE IMPORTANCE ANALYSIS\n")
cat("=================================\n")

rim <- calc.relimp(
  lm_model,
  type=c(
    "lmg",
    "last",
    "first",
    "betasq",
    "pratt",
    "genizi",
    "car"
  ),
  rela=TRUE
)

# =========================================================
# 17. CREATE RESULT TABLE
# =========================================================
rim_table <- data.frame(
  GCM=names(rim$lmg),
  LMG=as.numeric(rim$lmg),
  LAST=as.numeric(rim$last),
  FIRST=as.numeric(rim$first),
  GENIZI=as.numeric(rim$genizi),
  CAR=as.numeric(rim$car),
  BETASQ=as.numeric(rim$betasq),
  PRATT=as.numeric(rim$pratt)
)

cat("\nRIM RESULTS:\n")
print(rim_table)

# =========================================================
# 18. MCDA ANALYSIS
# =========================================================
pos <- c(
  "LMG",
  "FIRST",
  "GENIZI",
  "CAR",
  "BETASQ",
  "PRATT"
)

neg <- c("LAST")

rim_mat <- rim_table[, c(pos, neg), drop=FALSE]

# Normalize
rim_norm <- scale(
  rim_mat,
  center=FALSE,
  scale=sqrt(colSums(rim_mat^2))
)

# Ideal solution
ideal <- c(
  apply(rim_norm[, pos, drop=FALSE],2,max),
  apply(rim_norm[, neg, drop=FALSE],2,min)
)

# Distance
dist <- apply(
  rim_norm,
  1,
  function(x){
    sqrt(sum((x - ideal)^2))
  }
)

rim_table$MCDA_score <- 1/(1 + dist)

# =========================================================
# 19. PARETO ANALYSIS
# =========================================================
X <- rim_table[, c(pos, neg), drop=FALSE]

is_dominated <- function(i,X){
  
  any(apply(X,1,function(x){
    
    all(x >= X[i,]) &
      any(x > X[i,])
    
  }))
}

rim_table$Pareto <- sapply(
  1:nrow(X),
  function(i) !is_dominated(i,X)
)

# =========================================================
# 20. TOP 3 MODELS
# =========================================================
top3 <- rim_table %>%
  filter(Pareto == TRUE) %>%
  arrange(desc(MCDA_score)) %>%
  head(3)

cat("\n=================================\n")
cat("TOP 3 MODELS\n")
cat("=================================\n")

print(top3)

# =========================================================
# 21. SAVE RESULTS
# =========================================================
write.csv(
  rim_table,
  "FULL_RESULTS.csv",
  row.names=FALSE
)

write.csv(
  top3,
  "TOP3_MODELS.csv",
  row.names=FALSE
)

# =========================================================
# 22. HEATMAP DATA
# =========================================================
heatmap_data <- melt(
  rim_table[, c(
    "GCM",
    "LMG",
    "FIRST",
    "GENIZI",
    "CAR",
    "LAST",
    "BETASQ",
    "PRATT"
  )],
  id.vars="GCM"
)

# =========================================================
# 23. HEATMAP PLOT
# =========================================================
p1 <- ggplot(
  heatmap_data,
  aes(GCM, variable, fill=value)
) +
  geom_tile() +
  geom_text(
    aes(label=round(value,2)),
    size=3
  ) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    )
  )

# =========================================================
# 24. PARETO PLOT
# =========================================================
p2 <- ggplot(
  rim_table,
  aes(LMG, CAR, color=Pareto)
) +
  geom_point(size=4) +
  theme_minimal()

# =========================================================
# 25. MCDA BARPLOT
# =========================================================
p3 <- ggplot(
  rim_table,
  aes(
    reorder(GCM,MCDA_score),
    MCDA_score
  )
) +
  geom_col(fill="steelblue") +
  coord_flip() +
  theme_minimal()

# =========================================================
# 26. FINAL FIGURE
# =========================================================
final_plot <- (p1 / p2) | p3

print(final_plot)

ggsave(
  "JOURNAL_FIGURE.png",
  final_plot,
  width=14,
  height=8,
  dpi=600
)

# =========================================================
# 27. CLUSTER HEATMAP
# =========================================================
mat <- as.matrix(
  rim_table[, c(
    "LMG",
    "FIRST",
    "GENIZI",
    "CAR",
    "LAST",
    "BETASQ",
    "PRATT"
  )]
)

rownames(mat) <- rim_table$GCM

pheatmap(
  mat,
  clustering_method="complete",
  fontsize_row=8,
  fontsize_col=10,
  main="Clustered RIM Heatmap"
)

# =========================================================
# 28. FINAL RESULTS
# =========================================================
cat("\n=================================\n")
cat("FINAL TOP MODELS\n")
cat("=================================\n")

for(i in 1:nrow(top3)){
  
  cat(
    paste0(
      i,
      ". ",
      top3$GCM[i],
      " | MCDA = ",
      round(top3$MCDA_score[i],3),
      "\n"
    )
  )
}

cat("\nANALYSIS COMPLETED SUCCESSFULLY.\n")
# =========================================================
# JOURNAL-READY FIGURE (MCDA + PARETO + HEATMAP + TOP 5)
# =========================================================

library(ggplot2)
library(reshape2)
library(patchwork)
library(viridis)
library(dplyr)

rim_table <- as.data.frame(rim_table)

# =========================================================
# 1. METRIC COLUMNS (FIXED CASE FROM YOUR DATA)
# =========================================================

metric_cols <- c("lmg","first","genizi","car","betasq","pratt")
metric_cols <- metric_cols[metric_cols %in% names(rim_table)]

# =========================================================
# 2. MCDA SCORE
# =========================================================

rim_table$MCDA_score <- rowMeans(
  rim_table[, metric_cols, drop = FALSE],
  na.rm = TRUE
)

# =========================================================
# 3. TOP 5 MODELS
# =========================================================

top5_models <- rim_table$GCM[
  order(rim_table$MCDA_score, decreasing = TRUE)
][1:5]

rim_table$Top5 <- rim_table$GCM %in% top5_models
rim_table$Group <- ifelse(rim_table$Top5, "Top 5", "Others")

# =========================================================
# 4. GCM ORDER
# =========================================================

gcm_order <- rim_table$GCM[
  order(rim_table$MCDA_score, decreasing = TRUE)
]
gcm_order <- unique(as.character(gcm_order))

rim_table$GCM <- factor(rim_table$GCM, levels = gcm_order)

# =========================================================
# 5. HEATMAP DATA (CLEAN)
# =========================================================

heatmap_data <- melt(
  rim_table[, c("GCM", metric_cols)],
  id.vars = "GCM"
)

colnames(heatmap_data) <- c("GCM", "Metric", "Value")

heatmap_data$GCM <- factor(heatmap_data$GCM, levels = gcm_order)
heatmap_data$Metric <- factor(heatmap_data$Metric, levels = metric_cols)

# =========================================================
# 6. TOP 5 BOX LABELS
# =========================================================

top5_labels <- subset(rim_table, Top5)

# =========================================================
# 7. HEATMAP (NATURE STYLE)
# =========================================================

p1 <- ggplot(heatmap_data, aes(x = GCM, y = Metric, fill = Value)) +
  
  geom_tile(color = "white", linewidth = 0.6) +
  
  geom_text(aes(label = sprintf("%.2f", Value)),
            size = 3, fontface = "bold") +
  
  scale_fill_viridis_c(option = "C") +
  
  labs(title = "A) MCDA Component Heatmap",
       x = NULL, y = NULL) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# =========================================================
# 8. PARETO PLOT (TOP 5 HIGHLIGHT + BOX LABEL)
# =========================================================

p2 <- ggplot(rim_table, aes(x = lmg, y = car)) +
  
  geom_point(aes(color = Group), size = 4) +
  
  scale_color_manual(values = c("Top 5" = "red", "Others" = "grey70")) +
  
  geom_label(
    data = top5_labels,
    aes(label = GCM),
    fill = "white",
    color = "black",
    fontface = "bold",
    size = 3,
    linewidth = 0.3
  ) +
  
  labs(title = "B) Pareto Optimality (Top 5 Highlighted)",
       x = "LMG", y = "CAR") +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

# =========================================================
# 9. MCDA BAR (TOP 5 BOXED)
# =========================================================

p3 <- ggplot(rim_table, aes(x = GCM, y = MCDA_score)) +
  
  geom_col(aes(fill = Group), width = 0.7) +
  
  scale_fill_manual(values = c("Top 5" = "red", "Others" = "steelblue")) +
  
  geom_label(
    data = top5_labels,
    aes(label = round(MCDA_score, 2)),
    fill = "white",
    fontface = "bold",
    size = 3,
    linewidth = 0.3
  ) +
  
  coord_flip() +
  
  labs(title = "C) MCDA Ranking of GCMs",
       x = NULL, y = "MCDA Score") +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# =========================================================
# 10. FINAL COMPOSITE FIGURE (NATURE STYLE LAYOUT)
# =========================================================

final_plot <- (p1 / p2) | p3

print(final_plot)

# =========================================================
# 11. EXPORT (600 DPI PUBLICATION QUALITY)
# =========================================================

ggsave("Figure_MCDA_Pareto_Heatmap.png",
       final_plot,
       width = 18,
       height = 10,
       dpi = 600,
       bg = "white")

ggsave("Figure_MCDA_Pareto_Heatmap.pdf",
       final_plot,
       width = 18,
       height = 10)
print(p1)
ggsave("heatmap.png", p1, width = 10, height = 6, dpi = 600)
print(p2)
ggsave("pareto.png", p2, width = 8, height = 6, dpi = 600)
print(p3)
ggsave("mcda.png", p3, width = 8, height = 6, dpi = 600)
# ============================================================
# RIM HEATMAP INCLUDING MCDA SCORE ROW
# UBNB, Ethiopia
# ============================================================

# Install and load packages
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require(viridis)) {
  install.packages("viridis")
  library(viridis)
}

# ============================================================
# 1. INPUT YOUR RIM + MCDA VALUES
# ============================================================

heatmap_data <- data.frame(
  
  Variable = c(
    rep("MCDA", 12),
    rep("pratt", 12),
    rep("betasq", 12),
    rep("last", 12),
    rep("car", 12),
    rep("genizi", 12),
    rep("first", 12),
    rep("lmg", 12)
  ),
  
  GCM = rep(
    c(
      "ACCESS-CM2",
      "CMCC-ESM2",
      "EC-Earth3",
      "EC-Earth3-Veg-LR",
      "FGOALS-g3",
      "INM-CM4-8",
      "MIROC6",
      "MPI-ESM1-2-HR",
      "MPI-ESM1-2-LR",
      "NorESM2-LM",
      "NorESM2-MM",
      "TaiESM1"
    ),
    times = 8
  ),
  
  Value = c(
    # MCDA score
    0.13, 0.15, 0.16, 0.17, 0.16, 0.14, 0.18, 0.14, 0.14, 0.12, 0.14, 0.19,
    
    # pratt
    0.06, 0.08, 0.06, 0.12, 0.09, 0.06, 0.14, 0.02, 0.04, 0.03, 0.04, 0.18,
    
    # betasq
    0.00, 0.01, 0.01, 0.02, 0.01, 0.00, 0.03, 0.00, 0.00, 0.00, 0.00, 0.04,
    
    # last
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01,
    
    # car
    0.06, 0.08, 0.08, 0.09, 0.08, 0.07, 0.10, 0.06, 0.06, 0.06, 0.07, 0.11,
    
    # genizi
    0.07, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07, 0.07, 0.09,
    
    # first
    0.66, 0.74, 0.79, 0.77, 0.75, 0.72, 0.83, 0.73, 0.71, 0.63, 0.70, 0.82,
    
    # lmg
    0.06, 0.08, 0.08, 0.08, 0.08, 0.07, 0.09, 0.07, 0.07, 0.06, 0.07, 0.09
  )
)

# ============================================================
# 2. ORDER GCMS BY MCDA SCORE
# Highest MCDA score appears first
# ============================================================

gcm_order <- heatmap_data %>%
  filter(Variable == "MCDA") %>%
  arrange(desc(Value)) %>%
  pull(GCM)

heatmap_data$GCM <- factor(
  heatmap_data$GCM,
  levels = gcm_order
)

heatmap_data$Variable <- factor(
  heatmap_data$Variable,
  levels = c(
    "MCDA",
    "pratt",
    "betasq",
    "last",
    "car",
    "genizi",
    "first",
    "lmg"
  )
)

# ============================================================
# 3. OUTPUT FOLDER
# ============================================================

output_folder <- "D:/Desktop/HWRM_Thesis/From Gh/MCDA_Heatmap_Result"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

write.csv(
  heatmap_data,
  file.path(output_folder, "RIM_heatmap_with_MCDA_values.csv"),
  row.names = FALSE
)

# ============================================================
# 4. DRAW HEATMAP
# ============================================================

p <- ggplot(
  heatmap_data,
  aes(x = GCM, y = Variable, fill = Value)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    color = "black",
    size = 3.2
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    name = "Value",
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +
  
  labs(
    title = "Relative Importance  Heatmap ",
    x = "GCMs",
    y = "RIMs"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      size = 11
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

print(p)

# ============================================================
# 5. SAVE HEATMAP
# ============================================================

ggsave(
  filename = file.path(output_folder, "RIM_heatmap_with_MCDA_score.png"),
  plot = p,
  width = 13,
  height = 7,
  dpi = 300
)

cat("\nDONE SUCCESSFULLY!\n")
cat("Heatmap saved in:\n")
cat(output_folder, "\n")
