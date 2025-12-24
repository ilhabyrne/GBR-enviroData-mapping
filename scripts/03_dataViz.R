# =============================================================================
# Environmental data exploration
# Correlation analysis, spatial patterns, and visualizations
# =============================================================================

library(tidyverse)
library(corrplot)
library(Hmisc)
library(ggrepel)
library(viridis)
library(patchwork)
library(RColorBrewer)


# Load data
setwd("~/Documents/GitHub/GBR-enviroData-mapping/03_viz")
dat <- read.csv("~/Documents/GitHub/GBR-enviroData-mapping/02_eReefs/synthesis/RRAP_GBR1_metadata_targetSpp_2025-12-09_clean.csv")

# -----------------------------------------------------------------------------
# 1. DEFINE VARIABLE GROUPS
# -----------------------------------------------------------------------------

# Temperature variables
temp_vars <- c("mean_temp", "sd_temp", "skew_temp", "kurtosis_temp", 
               "MaxMonthMean", "MinMonthMean", "mean_daily_temp_range", 
               "mean_annual_temp_range", "temporal_autocorr_daily", 
               "temporal_autocorr_annual", "slope_annual_temp_change",
               "mean_duration_DHD", "mean_duration_DCD")

# Wind variables
wind_vars <- c("mean_wspeed_u", "mean_wspeed_v", "mean_wind_speed", 
               "sd_wind_speed", "cv_wind_speed", "min_wind_speed", 
               "max_wind_speed", "wind_speed_range",
               "mean_wind_direction")

# Current variables
current_vars <- c("mean_u", "mean_v", "mean_current_speed", "sd_current_speed",
                  "cv_current_speed", "min_current_speed", "max_current_speed",
                  "current_speed_range", "mean_current_direction")

# Get all numeric variables present in the data
all_numeric <- names(dat)[sapply(dat, is.numeric)]
# Filter to only include environmental variables (exclude lat/lon/year etc)
env_vars <- all_numeric[all_numeric %in% c(temp_vars, wind_vars, current_vars)]

# -----------------------------------------------------------------------------
# 2. CORRELATION ANALYSIS
# -----------------------------------------------------------------------------

# Get unique locations to avoid pseudoreplication in correlation analysis
dat_unique <- dat %>%
  distinct(decimalLatitude, .keep_all = TRUE) %>%
  dplyr::select(all_of(c("locality", "decimalLatitude", env_vars)))

# --- 2a. Full correlation matrix ---
cor_matrix_full <- cor(dat_unique[, env_vars], use = "pairwise.complete.obs")

# Plot full correlation matrix
png("correlation_matrix_full.png", width = 14, height = 14, units = "in", res = 150)
corrplot(cor_matrix_full, method = "color", type = "upper", 
         order = "hclust", tl.cex = 1.5, tl.col = "black",
         col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
         mar = c(0, 0, 2, 0))
dev.off()

# --- 2b. Temperature-only correlation matrix ---
temp_vars_present <- temp_vars[temp_vars %in% names(dat_unique)]
cor_matrix_temp <- cor(dat_unique[, temp_vars_present], use = "pairwise.complete.obs")

png("correlation_matrix_temperature.png", width = 10, height = 10, units = "in", res = 150)
corrplot(cor_matrix_temp, method = "color", type = "upper", 
         order = "hclust", tl.cex = 1, tl.col = "black",
         col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
         mar = c(0, 0, 2, 0))
dev.off()

# --- 2c. Wind-only correlation matrix ---
wind_vars_present <- wind_vars[wind_vars %in% names(dat_unique)]
cor_matrix_wind <- cor(dat_unique[, wind_vars_present], use = "pairwise.complete.obs")

png("correlation_matrix_wind.png", width = 10, height = 10, units = "in", res = 150)
corrplot(cor_matrix_wind, method = "color", type = "upper", 
         order = "hclust", tl.cex = 1, tl.col = "black",
         col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
         mar = c(0, 0, 2, 0))
dev.off()

# --- 2d. Current-only correlation matrix ---
current_vars_present <- current_vars[current_vars %in% names(dat_unique)]
cor_matrix_current <- cor(dat_unique[, current_vars_present], use = "pairwise.complete.obs")

png("correlation_matrix_current.png", width = 10, height = 10, units = "in", res = 150)
corrplot(cor_matrix_current, method = "color", type = "upper", 
         order = "hclust", tl.cex = 1, tl.col = "black",
         col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
         mar = c(0, 0, 2, 0))
dev.off()

# --- 2e. Identify highly correlated pairs ---
find_high_correlations <- function(cor_mat, threshold = 0.7) {
  high_cor <- which(abs(cor_mat) > threshold & abs(cor_mat) < 1, arr.ind = TRUE)
  # Only keep upper triangle to avoid duplicates
  high_cor <- high_cor[high_cor[,1] < high_cor[,2], , drop = FALSE]
  
  if(nrow(high_cor) == 0) {
    return(data.frame(var1 = character(), var2 = character(), correlation = numeric()))
  }
  
  results <- data.frame(
    var1 = rownames(cor_mat)[high_cor[,1]],
    var2 = colnames(cor_mat)[high_cor[,2]],
    correlation = cor_mat[high_cor]
  ) %>%
    arrange(desc(abs(correlation)))
  
  return(results)
}

high_cor_pairs <- find_high_correlations(cor_matrix_full, threshold = 0.7)
print(high_cor_pairs, n = 50)

# Save to file
write.csv(high_cor_pairs, "high_correlation_pairs.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. SPATIAL PATTERNS - TEMPERATURE ACROSS LATITUDE
# -----------------------------------------------------------------------------

# Summarize by site (unique lat/lon)
site_summary <- dat %>%
  group_by(locality, decimalLatitude) %>%
  summarise(
    n_samples = n(),
    mean_temp = first(mean_temp),
    sd_temp = first(sd_temp),
    MaxMonthMean = first(MaxMonthMean),
    MinMonthMean = first(MinMonthMean),
    mean_annual_temp_range = first(mean_annual_temp_range),
    .groups = "drop"
  )

# --- 3a. Mean temperature vs latitude (red = hot, blue = cold) ---
p_temp_lat <- ggplot(site_summary, aes(x = decimalLatitude, y = mean_temp)) +
  geom_point(aes(size = n_samples, color = mean_temp), alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "gray30", fill = "gray70", alpha = 0.2) +
  scale_x_reverse() +
  scale_color_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"),
    name = "Mean Temp (°C)"
  ) +
  labs(x = "Latitude", y = "Mean Temperature (°C)", 
       size = "N samples") +
  theme_bw() +
  theme(legend.position = "right")
p_temp_lat

# --- 3b. Multiple temperature metrics (red = hot, blue = cold) ---
temp_long <- site_summary %>%
  dplyr::select(locality, decimalLatitude, mean_temp, MaxMonthMean, MinMonthMean) %>%
  pivot_longer(cols = c(mean_temp, MaxMonthMean, MinMonthMean),
               names_to = "metric", values_to = "temperature") %>%
  mutate(metric = factor(metric, 
                         levels = c("MaxMonthMean", "mean_temp", "MinMonthMean"),
                         labels = c("Max Month Mean", "Mean Temp", "Min Month Mean")))

p_temp_metrics <- ggplot(temp_long, aes(x = decimalLatitude, y = temperature)) +
  geom_point(aes(color = temperature), alpha = 0.8, size = 3) +
  geom_smooth(aes(linetype = metric), method = "lm", se = FALSE, 
              linewidth = 0.8, color = "gray30") +
  scale_x_reverse() +
  scale_color_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"),
    name = "Temperature (°C)"
  ) +
  #scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = "Latitude", y = "Temperature (°C)", 
       shape = "Metric", linetype = "Metric") +
  theme_bw() +
  theme(legend.position = "right")
p_temp_metrics

ggsave(plot = p_temp_metrics, "mean_temp_metrics.png", dpi = 300, width = 6, height = 4)

# --- 3c. Temperature by site (red = hot, blue = cold) ---
site_summary <- dat %>%
  group_by(locality) %>%
  summarise(
    decimalLatitude = mean(decimalLatitude, na.rm = TRUE),
    n_samples = n(),
    mean_temp = mean(mean_temp, na.rm = TRUE),
    sd_temp = mean(sd_temp, na.rm = TRUE),
    .groups = "drop"
  )

site_order <- site_summary %>%
  arrange(desc(decimalLatitude)) %>%
  pull(locality)

site_summary$locality_ordered <- factor(site_summary$locality, levels = site_order)

site_subset <- site_summary %>%
  filter(!locality %in% c("Dungeness", "Tydeman", "MacGillivray",
                          "Myrmidon", "Magnetic_Island", "North_Keppel",
                          "Miall", "Middle", "Halfway", "Great_Keppel", "Lizard", "Moore_Reef"))

p_temp_site <- ggplot(site_subset, aes(x = locality_ordered, y = mean_temp)) +
  geom_point(aes(size = n_samples, color = mean_temp)) +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp), 
                width = 0.2, alpha = 0.5) +
  scale_color_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"),
    name = "Mean Temp (°C)") +
  labs(x = "Reef (ordered N to S)", y = "Mean temperature (°C)", 
       subtitle = "",
       size = "N samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
p_temp_site

ggsave(plot = p_temp_site, "mean_SD_temp.png", dpi = 300, width = 6, height = 4)

# -----------------------------------------------------------------------------
# 4. CREATE COLOUR PALETTE
# -----------------------------------------------------------------------------

colR <- colorRampPalette(brewer.pal(11, "RdYlBu"))(29)

# Override problem colors for better distinction
colR[15] <- "#A1D99B"
colR[16] <- "#74C476"
colR[17] <- "#41AB5D"
colR[18] <- "#6BAED6"
colR[19] <- "#4292C6"
colR[27] <- "#4292C6"

# Named vector matching locality column names
locality_colors <- c(
  "Dungeness"         = colR[1],
  "Masig"             = colR[2],
  "Aukane"            = colR[3],
  "Hicks"             = colR[4],
  "No_Name"           = colR[5],
  "Lizard"            = colR[6],
  "Martin"            = colR[6],
  "North_Direction"   = colR[7],
  "Mackay"            = colR[8],
  "St_Crispin"        = colR[9],
  "Moore_Reef"        = colR[10],
  "Fitzroy_Island"    = colR[11],
  "Kelso"             = colR[12],
  "Pelorus"           = colR[13],
  "Orpheus"           = colR[14],
  "Chicken"           = colR[15],
  "Davies"            = colR[16],
  "Little_Broadhurst" = colR[17],
  "East_Cay"          = colR[19],
  "Reef21-550"        = colR[20],
  "North_Keppel"      = colR[21],
  "Miall"             = colR[22],
  "Middle"            = colR[23],
  "Great_Keppel"      = colR[24],
  "Halfway"           = colR[25],
  "Heron"             = colR[26],
  "Sykes"             = colR[27],
  "Fitzroy_Reef"      = colR[28],
  "Lady_Musgrave"     = colR[29]
)

# =============================================================================
# 5. PREPARE DATA - UNIQUE SITES
# =============================================================================

# Get unique locations to avoid pseudoreplication
dat_unique <- dat %>%
  distinct(locality, .keep_all = TRUE) %>%
  dplyr::select(locality, decimalLatitude, all_of(env_vars))

cat("Unique localities:", nrow(dat_unique), "\n")

# =============================================================================
# 5. DIAGNOSE MISSING DATA
# =============================================================================

cat("\n=== Missing Data Diagnosis ===\n")

# Check for NAs in each variable
na_counts <- sapply(dat_unique[, env_vars], function(x) sum(is.na(x)))
if (any(na_counts > 0)) {
  cat("Variables with missing data:\n")
  print(na_counts[na_counts > 0])
}

# Check for infinite values
inf_counts <- sapply(dat_unique[, env_vars], function(x) sum(is.infinite(x)))
if (any(inf_counts > 0)) {
  cat("Variables with infinite values:\n")
  print(inf_counts[inf_counts > 0])
}

# Summary
cat("\nTotal rows:", nrow(dat_unique), "\n")
cat("Complete cases:", sum(complete.cases(dat_unique[, env_vars])), "\n")

# =============================================================================
# 5. HANDLE MISSING DATA FOR PCA
# =============================================================================

# Option A: Keep only variables with <10% missing
na_prop <- na_counts / nrow(dat_unique)
good_vars <- env_vars[na_prop < 0.1]
cat("\nVariables with <10% missing:", length(good_vars), "of", length(env_vars), "\n")

# Option B: Create complete dataset
dat_unique_complete <- dat_unique %>%
  dplyr::select(locality, decimalLatitude, all_of(good_vars)) %>%
  filter(complete.cases(across(all_of(good_vars))))

cat("Complete cases for PCA:", nrow(dat_unique_complete), "\n")

# =============================================================================
# 5. SITE SUMMARY FOR PLOTS
# =============================================================================

site_summary <- dat %>%
  group_by(locality) %>%
  summarise(
    decimalLatitude = mean(decimalLatitude, na.rm = TRUE),
    n_samples = n(),
    mean_temp = first(mean_temp),
    sd_temp = first(sd_temp),
    MaxMonthMean = first(MaxMonthMean),
    MinMonthMean = first(MinMonthMean),
    .groups = "drop"
  )

cat("Site summary created:", nrow(site_summary), "sites\n")

# =============================================================================
# 5. PCA ANALYSIS
# =============================================================================

cat("\n=== PCA Analysis ===\n")

# Run PCA on complete data
pca_result <- prcomp(dat_unique_complete[, good_vars], center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_result)

cat("PCA Variance Explained:\n")
print(pca_summary$importance[, 1:5])

# Extract scores
pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(locality = dat_unique_complete$locality)

# Extract loadings
pca_loadings <- as.data.frame(pca_result$rotation) %>%
  mutate(variable = rownames(.))

# Scale factor for arrows
scale_factor <- min(
  (max(pca_scores$PC1) - min(pca_scores$PC1)) / (max(abs(pca_loadings$PC1)) * 2),
  (max(pca_scores$PC2) - min(pca_scores$PC2)) / (max(abs(pca_loadings$PC2)) * 2)
)

# Top loadings
top_loadings <- pca_loadings %>%
  mutate(loading_magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(loading_magnitude, n = 15)

# Filter out non-focal sites
pca_scores_focal <- pca_scores %>%
  filter(!locality %in% c("Dungeness", "Tydeman", "MacGillivray",
                          "Myrmidon", "Magnetic_Island", "North_Keppel",
                          "Miall", "Middle", "Halfway", "Great_Keppel",
                          "Lizard", "Moore_Reef"))

# Get latitude order (N to S)
locality_lat_order <- dat %>%
  group_by(locality) %>%
  summarise(lat = mean(decimalLatitude, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(lat)) %>%
  pull(locality)

# Apply factor ordering
pca_scores_focal <- pca_scores_focal %>%
  mutate(locality = factor(locality, levels = locality_lat_order))

# --- 5a. PCA biplot by reef ---
p_pca <- ggplot() +
  geom_point(data = pca_scores_focal, aes(x = PC1, y = PC2, color = locality),
             alpha = 0.8, size = 4) +
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0, xend = PC1 * scale_factor * 0.8, yend = PC2 * scale_factor * 0.8),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray40", alpha = 0.7) +
  geom_text_repel(data = top_loadings,
                  aes(x = PC1 * scale_factor * 0.85, y = PC2 * scale_factor * 0.85, label = variable),
                  size = 2.5, color = "gray30", max.overlaps = 20) +
  scale_color_manual(values = locality_colors, name = "Reef") +
  labs(x = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
       title = "PCA of Environmental Variables",
       subtitle = "Top 15 loadings shown") +
  theme_bw() +
  coord_fixed() +
  guides(color = guide_legend(ncol = 1))

p_pca
ggsave("pca_biplot_focal.png", p_pca, width = 10, height = 8, dpi = 300)

# --- 5b. Scree plot ---
var_explained <- data.frame(
  PC = paste0("PC", 1:length(pca_summary$importance[2,])),
  variance = pca_summary$importance[2,] * 100,
  cumulative = pca_summary$importance[3,] * 100
) %>%
  mutate(PC = factor(PC, levels = PC))

p_scree <- ggplot(var_explained[1:10,], aes(x = PC)) +
  geom_col(aes(y = variance), fill = "#0072B2", alpha = 0.7) +
  geom_line(aes(y = cumulative, group = 1), color = "#D55E00", linewidth = 1) +
  geom_point(aes(y = cumulative), color = "#D55E00", size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
  labs(x = "Principal Component", y = "Variance Explained (%)",
       title = "PCA Scree Plot",
       subtitle = "Bars = individual, Line = cumulative") +
  theme_bw()

ggsave("pca_scree.png", p_scree, width = 8, height = 5, dpi = 150)

# --- 5c. PCA loadings by variable type ---
pca_loadings <- pca_loadings %>%
  mutate(var_type = case_when(
    variable %in% temp_vars ~ "Temperature",
    variable %in% wind_vars ~ "Wind",
    variable %in% current_vars ~ "Current",
    TRUE ~ "Other"
  ))

p_loadings <- ggplot(pca_loadings, aes(x = PC1, y = PC2, color = var_type)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = variable), size = 2, max.overlaps = 20) +
  scale_color_manual(values = c("Temperature" = "#D55E00", "Current" = "#0072B2",
                                "Wind" = "#009E73", "Other" = "gray50")) +
  labs(x = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
       title = "",
       color = "Variable type") +
  theme_bw() +
  coord_fixed()

p_loadings
ggsave("pca_loadings_by_type.png", p_loadings, width = 5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 6. WIND ROSE VISUALIZATION
# -----------------------------------------------------------------------------

if (all(c("mean_wspeed_u", "mean_wspeed_v", "mean_wind_speed") %in% names(dat))) {
  
  # Calculate wind direction from u, v components
  wind_vectors <- dat_unique %>%
    dplyr::select(locality, decimalLatitude, mean_wspeed_u, mean_wspeed_v, mean_wind_speed) %>%
    drop_na() %>%
    mutate(
      wind_direction = (270 - atan2(mean_wspeed_v, mean_wspeed_u) * 180/pi) %% 360
    ) %>%
    filter(!locality %in% c("Dungeness", "Tydeman", "MacGillivray",
                            "Myrmidon", "Magnetic_Island", "North_Keppel",
                            "Miall", "Middle", "Halfway", "Great_Keppel",
                            "Lizard", "Moore_Reef"))
  
  # Order by latitude
  wind_vectors <- wind_vectors %>%
    mutate(locality = factor(locality, levels = locality_lat_order))
  
  # U-V scatter plot with arrows
  p_wind_uv <- ggplot(wind_vectors, aes(x = mean_wspeed_u, y = mean_wspeed_v, color = locality)) +
    geom_segment(aes(x = 0, y = 0, xend = mean_wspeed_u, yend = mean_wspeed_v),
                 arrow = arrow(length = unit(0.2, "cm")), alpha = 0.8) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    scale_color_manual(values = locality_colors, name = "Locality") +
    labs(x = "Mean U Component (E-W)", y = "Mean V Component (N-S)",
         title = "Mean Wind Vectors by Locality",
         subtitle = "Arrows from origin show mean wind direction and magnitude") +
    theme_bw() +
    coord_fixed() +
    guides(color = guide_legend(ncol = 1))
  
  ggsave("wind_uv_vectors.png", p_wind_uv, width = 10, height = 8, dpi = 150)
  
  # Polar plot
  p_wind_polar <- ggplot(wind_vectors, aes(x = wind_direction, y = mean_wind_speed, color = locality)) +
    geom_point(size = 4, alpha = 0.8) +
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(0, 315, 45),
                       labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                       limits = c(0, 360)) +
    scale_color_manual(values = locality_colors, name = "Locality") +
    labs(title = "Mean Wind Vector by Locality",
         subtitle = "Angle = direction, Distance from center = speed",
         y = "Mean Wind Speed") +
    theme_minimal() +
    guides(color = guide_legend(ncol = 1))
  
  ggsave("wind_polar_plot.png", p_wind_polar, width = 10, height = 10, dpi = 150)
}

# -----------------------------------------------------------------------------
# WIND ROSES BY LOCALITY
# -----------------------------------------------------------------------------

if(all(c("mean_wind_direction", "mean_wind_speed") %in% names(dat))) {
  
  wind_data <- dat_unique %>%
    dplyr::select(locality, decimalLatitude, mean_wind_direction, mean_wind_speed) %>%
    drop_na() %>%
    mutate(
      direction_sector = cut(mean_wind_direction, 
                             breaks = seq(-22.5, 337.5, 45),
                             labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
      speed_cat = cut(mean_wind_speed, 
                      breaks = quantile(mean_wind_speed, probs = c(0, 0.25, 0.5, 0.75, 1)),
                      labels = c("Low", "Med-Low", "Med-High", "High"),
                      include.lowest = TRUE)
    ) %>%
    mutate(direction_sector = case_when(
      mean_wind_direction > 337.5 | mean_wind_direction <= 22.5 ~ "N",
      TRUE ~ as.character(direction_sector)
    )) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  # Filter to focal localities
  wind_data_focal <- wind_data %>%
    filter(!locality %in% c("Dungeness", "Tydeman", "MacGillivray",
                            "Myrmidon", "Magnetic_Island", "North_Keppel",
                            "Miall", "Middle", "Halfway", "Great_Keppel", 
                            "Lizard", "Moore_Reef"))
  
  # Order by latitude
  locality_lat_order <- wind_data_focal %>%
    group_by(locality) %>%
    summarise(lat = mean(decimalLatitude, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(lat)) %>%
    pull(locality)
  
  wind_data_focal <- wind_data_focal %>%
    mutate(locality = factor(locality, levels = locality_lat_order))
  
  # Create summary for each locality - need all direction sectors represented
  all_sectors <- expand.grid(
    locality = unique(wind_data_focal$locality),
    direction_sector = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  ) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  wind_summary_by_site <- wind_data_focal %>%
    count(locality, direction_sector) %>%
    right_join(all_sectors, by = c("locality", "direction_sector")) %>%
    mutate(n = replace_na(n, 0)) %>%
    group_by(locality) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Faceted wind rose
  p_wind_rose_facet <- ggplot(wind_summary_by_site, aes(x = direction_sector, y = proportion)) +
    geom_col(fill = "#0072B2", alpha = 0.7, width = 0.8) +
    coord_polar(start = 0) +
    facet_wrap(~ locality, ncol = 5) +
    labs(title = "",
         subtitle = "Ordered N to S",
         y = "Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    )
  
  ggsave("wind_rose_by_locality.png", p_wind_rose_facet, width = 7, height = 8, dpi = 300)
  
  # Alternative: color by locality instead of faceting
  p_wind_locality_color <- ggplot(wind_data_focal, 
                                  aes(x = mean_wind_direction, fill = locality)) +
    geom_histogram(binwidth = 22.5, boundary = 0, alpha = 0.7) +
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(0, 315, 45), 
                       labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                       limits = c(0, 360)) +
    scale_fill_manual(values = locality_colors, name = "Locality") +
    labs(title = "",
         x = "", y = "Count") +
    theme_minimal() +
    guides(fill = guide_legend(ncol = 2))
  
  ggsave("wind_direction_by_locality_stacked.png", p_wind_locality_color, 
         width = 10, height = 10, dpi = 300)
}

# -----------------------------------------------------------------------------
# CURRENT ROSES BY LOCALITY
# -----------------------------------------------------------------------------

if(all(c("mean_current_direction", "mean_current_speed") %in% names(dat))) {
  
  current_data <- dat_unique %>%
    dplyr::select(locality, decimalLatitude, mean_current_direction, mean_current_speed) %>%
    drop_na() %>%
    mutate(
      direction_sector = cut(mean_current_direction, 
                             breaks = seq(-22.5, 337.5, 45),
                             labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
      speed_cat = cut(mean_current_speed, 
                      breaks = quantile(mean_current_speed, probs = c(0, 0.25, 0.5, 0.75, 1)),
                      labels = c("Low", "Med-Low", "Med-High", "High"),
                      include.lowest = TRUE)
    ) %>%
    mutate(direction_sector = case_when(
      mean_current_direction > 337.5 | mean_current_direction <= 22.5 ~ "N",
      TRUE ~ as.character(direction_sector)
    )) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  # Filter to focal localities
  current_data_focal <- current_data %>%
    filter(!locality %in% c("Dungeness", "Tydeman", "MacGillivray",
                            "Myrmidon", "Magnetic_Island", "North_Keppel",
                            "Miall", "Middle", "Halfway", "Great_Keppel", 
                            "Lizard", "Moore_Reef"))
  
  # Order by latitude
  locality_lat_order <- current_data_focal %>%
    group_by(locality) %>%
    summarise(lat = mean(decimalLatitude, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(lat)) %>%
    pull(locality)
  
  current_data_focal <- current_data_focal %>%
    mutate(locality = factor(locality, levels = locality_lat_order))
  
  # Create summary with all sectors represented
  all_sectors_current <- expand.grid(
    locality = unique(current_data_focal$locality),
    direction_sector = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  ) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  current_summary_by_site <- current_data_focal %>%
    count(locality, direction_sector) %>%
    right_join(all_sectors_current, by = c("locality", "direction_sector")) %>%
    mutate(n = replace_na(n, 0)) %>%
    group_by(locality) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Faceted current rose
  p_current_rose_facet <- ggplot(current_summary_by_site, aes(x = direction_sector, y = proportion)) +
    geom_col(fill = "#009E73", alpha = 0.7, width = 0.8) +
    coord_polar(start = 0) +
    facet_wrap(~ locality, ncol = 5) +
    labs(title = "",
         subtitle = "Ordered N to S",
         y = "Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 6),
      axis.title.x = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    )
  
  ggsave("current_rose_by_locality.png", p_current_rose_facet, width = 7, height = 8, dpi = 300)
  
  # Alternative: color by locality
  p_current_locality_color <- ggplot(current_data_focal, 
                                     aes(x = mean_current_direction, fill = locality)) +
    geom_histogram(binwidth = 22.5, boundary = 0, alpha = 0.7) +
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(0, 315, 45), 
                       labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                       limits = c(0, 360)) +
    scale_fill_manual(values = locality_colors, name = "Locality") +
    labs(title = "",
         x = "", y = "Count") +
    theme_minimal() +
    guides(fill = guide_legend(ncol = 2))
  
  ggsave("current_direction_by_locality_stacked.png", p_current_locality_color, 
         width = 10, height = 10, dpi = 3000)
}

# -----------------------------------------------------------------------------
# 6. U-V VECTOR PLOTS (Wind and Current)
# -----------------------------------------------------------------------------

if(all(c("mean_wspeed_u", "mean_wspeed_v") %in% names(dat))) {
  
  # U-V scatter plot for wind
  p_wind_uv <- ggplot(dat_unique, aes(x = mean_wspeed_u, y = mean_wspeed_v)) +
    geom_point(aes(color = decimalLatitude), size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_c(option = "viridis") +
    labs(x = "Mean U (E-W) Component", y = "Mean V (N-S) Component",
         subtitle = "U = eastward (+), V = northward (+)",
         color = "Latitude") +
    theme_bw() +
    coord_fixed()
  
  ggsave("wind_uv_scatter.png", p_wind_uv, width = 8, height = 8, dpi = 150)
}

if(all(c("mean_u", "mean_v") %in% names(dat))) {
  
  # U-V scatter plot for currents
  p_current_uv <- ggplot(dat_unique, aes(x = mean_u, y = mean_v)) +
    geom_point(aes(color = decimalLatitude), size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_c(option = "viridis") +
    labs(x = "Mean U (E-W) Component", y = "Mean V (N-S) Component",
         subtitle = "U = eastward (+), V = northward (+)",
         color = "Latitude") +
    theme_bw() +
    coord_fixed()
  
  ggsave("current_uv_scatter.png", p_current_uv, width = 8, height = 8, dpi = 150)
}

# -----------------------------------------------------------------------------
# 7. Hierarchical clustering
# -----------------------------------------------------------------------------

# Hierarchical clustering of variables
var_dist <- as.dist(1 - abs(cor_matrix_full))
var_clust <- hclust(var_dist, method = "complete")

png("variable_dendrogram.png", width = 14, height = 8, units = "in", res = 150)
plot(var_clust, main = "",
     xlab = "", sub = "Based on correlation distance (1 - |r|)")
abline(h = 0.3, col = "red", lty = 2)  # Cut-off for r > 0.7
dev.off()
