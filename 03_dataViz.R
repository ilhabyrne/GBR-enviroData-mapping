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

# Load data
setwd("~/Documents/GitHub/GBR-enviroData-mapping/03_viz")
dat <- read.csv("~/Documents/GitHub/GBR-enviroData-mapping/02_eReefs/processed/RRAP_enviro_metadata_targetSpp_2025-12-09_clean.csv")

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
               "mean_wind_direction", "wind_directional_constancy", 
               "wind_directional_variability", "mean_wave_energy_proxy", "wind_dispersal_symmetry")

# Current variables
current_vars <- c("mean_u", "mean_v", "mean_current_speed", "sd_current_speed",
                  "cv_current_speed", "min_current_speed", "max_current_speed",
                  "current_speed_range", "mean_current_direction",
                  "current_transport_potential", "current_dispersal_symmetry",
                  "current_retention_potential")

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

ggsave(plot = p_temp_metrics, "mean_temp_metrics.png", dpi = 300, width = 8, height = 6)

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
  labs(x = "Reef (ordered N to S)", y = "Mean Temperature (°C)", 
       subtitle = "Error bars = temporal SD",
       size = "N samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
p_temp_site

ggsave(plot = p_temp_site, "mean_SD_temp.png", dpi = 300, width = 6, height = 4)

# -----------------------------------------------------------------------------
# 4. PCA OF ENVIRONMENTAL VARIABLES
# -----------------------------------------------------------------------------

# Prepare data for PCA (remove NAs and scale)
pca_data <- dat_unique %>%
  dplyr::select(locality, all_of(env_vars)) %>%
  drop_na()

# Run PCA
pca_result <- prcomp(pca_data[, env_vars], center = TRUE, scale. = TRUE)

# Summary
pca_summary <- summary(pca_result)
print(pca_summary$importance[, 1:5])

# Extract scores and loadings
pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(locality = pca_data$locality)

pca_loadings <- as.data.frame(pca_result$rotation) %>%
  mutate(variable = rownames(.))

# --- 4a. PCA biplot (PC1 vs PC2) ---
# Calculate scaling factor for loadings
scale_factor <- min(
  (max(pca_scores$PC1) - min(pca_scores$PC1)) / (max(abs(pca_loadings$PC1)) * 2),
  (max(pca_scores$PC2) - min(pca_scores$PC2)) / (max(abs(pca_loadings$PC2)) * 2)
)

# Select top loadings to display
top_loadings <- pca_loadings %>%
  mutate(loading_magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(loading_magnitude, n = 15)

p_pca <- ggplot() +
  # Sites
  geom_point(data = pca_scores, aes(x = PC1, y = PC2), 
             alpha = 0.7, size = 3, color = "#0072B2") +
  geom_text_repel(data = pca_scores, aes(x = PC1, y = PC2, label = locality),
                  size = 2.5, max.overlaps = 10) +
  # Loadings as arrows
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0, xend = PC1 * scale_factor * 0.8, yend = PC2 * scale_factor * 0.8),
               arrow = arrow(length = unit(0.2, "cm")), color = "#D55E00", alpha = 0.7) +
  geom_text_repel(data = top_loadings,
                  aes(x = PC1 * scale_factor * 0.85, y = PC2 * scale_factor * 0.85, label = variable),
                  size = 2, color = "#D55E00", max.overlaps = 20) +
  labs(x = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
       subtitle = "Top 15 loadings shown") +
  theme_bw() +
  coord_fixed()
p_pca

ggsave("pca_biplot.png", p_pca, width = 10, height = 10, dpi = 150)

# --- 4b. PCA by variable type ---
# Color loadings by variable type
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
  scale_color_manual(values = c("Temperature" = "#D55E00", "Wind" = "#0072B2", 
                                "Current" = "#009E73", "Other" = "gray50")) +
  labs(x = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
       color = "Variable Type") +
  theme_bw() +
  coord_fixed()

p_loadings
ggsave("pca_loadings_by_type.png", p_loadings, width = 5, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 5. WIND ROSE VISUALIZATION
# -----------------------------------------------------------------------------

# Check if wind direction and speed are available
if(all(c("mean_wind_direction", "mean_wind_speed") %in% names(dat))) {
  
  # For wind roses, we need the circular package
  # install.packages("circular") if needed
  
  # Simple wind rose using base R approach
  wind_data <- dat_unique %>%
    select(locality, mean_wind_direction, mean_wind_speed) %>%
    drop_na() %>%
    mutate(
      # Convert direction to cardinal sectors (8 sectors)
      direction_sector = cut(mean_wind_direction, 
                             breaks = seq(-22.5, 337.5, 45),
                             labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
      # Speed categories
      speed_cat = cut(mean_wind_speed, 
                      breaks = quantile(mean_wind_speed, probs = c(0, 0.25, 0.5, 0.75, 1)),
                      labels = c("Low", "Med-Low", "Med-High", "High"),
                      include.lowest = TRUE)
    )
  
  # If direction wraps around 360 to 0, handle the N sector
  wind_data <- wind_data %>%
    mutate(direction_sector = case_when(
      mean_wind_direction > 337.5 | mean_wind_direction <= 22.5 ~ "N",
      TRUE ~ as.character(direction_sector)
    ))
  
  # Summarize by direction
  wind_summary <- wind_data %>%
    count(direction_sector) %>%
    mutate(
      direction_sector = factor(direction_sector, 
                                levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")),
      proportion = n / sum(n)
    )
  
  # Wind rose plot
  p_wind_rose <- ggplot(wind_summary, aes(x = direction_sector, y = proportion)) +
    geom_col(fill = "#0072B2", alpha = 0.7, width = 0.8) +
    coord_polar(start = 0) +
    labs(title = "Wind Direction Distribution",
         subtitle = "Mean wind direction across sites",
         y = "Proportion of sites") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank()
    )
  
  ggsave("wind_rose_simple.png", p_wind_rose, width = 8, height = 8, dpi = 150)
  
  # --- More detailed wind rose with speed ---
  wind_detail <- wind_data %>%
    count(direction_sector, speed_cat) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  p_wind_rose_detail <- ggplot(wind_detail, aes(x = direction_sector, y = n, fill = speed_cat)) +
    geom_col(width = 0.8) +
    coord_polar(start = 0) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    labs(title = "Wind Rose by Direction and Speed",
         subtitle = "Site distribution",
         fill = "Wind Speed",
         y = "Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank()
    )
  
  ggsave("wind_rose_with_speed.png", p_wind_rose_detail, width = 8, height = 8, dpi = 150)
}

# Similar for current direction
if(all(c("mean_current_direction", "mean_current_speed") %in% names(dat))) {
  
  current_data <- dat_unique %>%
    select(locality, mean_current_direction, mean_current_speed) %>%
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
    ))
  
  current_summary <- current_data %>%
    count(direction_sector, speed_cat) %>%
    mutate(direction_sector = factor(direction_sector, 
                                     levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")))
  
  p_current_rose <- ggplot(current_summary, aes(x = direction_sector, y = n, fill = speed_cat)) +
    geom_col(width = 0.8) +
    coord_polar(start = 0) +
    scale_fill_viridis_d(option = "viridis", direction = -1) +
    labs(title = "Current Rose by Direction and Speed",
         subtitle = "Site distribution",
         fill = "Current Speed",
         y = "Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title.x = element_blank()
    )
  
  ggsave("current_rose.png", p_current_rose, width = 8, height = 8, dpi = 150)
}

# -----------------------------------------------------------------------------
# 6. U-V VECTOR PLOTS (Wind and Current)
# -----------------------------------------------------------------------------

if(all(c("mean_wspeed_u", "mean_wspeed_v") %in% names(dat))) {
  
  # U-V scatter plot for wind
  p_wind_uv <- ggplot(dat_unique, aes(x = mean_wspeed_u, y = mean_wspeed_v)) +
    geom_point(aes(color = decimalLatitudeOriginal), size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "Mean U (E-W) Component", y = "Mean V (N-S) Component",
         title = "Wind Vector Components",
         subtitle = "U = eastward (+), V = northward (+)",
         color = "Latitude") +
    theme_bw() +
    coord_fixed()
  
  ggsave("wind_uv_scatter.png", p_wind_uv, width = 8, height = 8, dpi = 150)
}

if(all(c("mean_u", "mean_v") %in% names(dat))) {
  
  # U-V scatter plot for currents
  p_current_uv <- ggplot(dat_unique, aes(x = mean_u, y = mean_v)) +
    geom_point(aes(color = decimalLatitudeOriginal), size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_c(option = "viridis") +
    labs(x = "Mean U (E-W) Component", y = "Mean V (N-S) Component",
         title = "Current Vector Components",
         subtitle = "U = eastward (+), V = northward (+)",
         color = "Latitude") +
    theme_bw() +
    coord_fixed()
  
  ggsave("current_uv_scatter.png", p_current_uv, width = 8, height = 8, dpi = 150)
}

# -----------------------------------------------------------------------------
# 7. SUMMARY STATISTICS AND VARIABLE SELECTION GUIDANCE
# -----------------------------------------------------------------------------

cat("\n=== SUMMARY AND RECOMMENDATIONS ===\n\n")

# Count correlations above different thresholds
for(thresh in c(0.9, 0.8, 0.7)) {
  n_pairs <- sum(abs(cor_matrix_full[upper.tri(cor_matrix_full)]) > thresh)
  cat(sprintf("Variable pairs with |r| > %.1f: %d\n", thresh, n_pairs))
}

cat("\n--- Suggested Variable Groups for Reduced Multicollinearity ---\n")
cat("Consider selecting ONE representative from each highly correlated cluster.\n")
cat("Run hierarchical clustering on the correlation matrix to identify clusters.\n\n")

# Hierarchical clustering of variables
var_dist <- as.dist(1 - abs(cor_matrix_full))
var_clust <- hclust(var_dist, method = "complete")

png("variable_dendrogram.png", width = 14, height = 8, units = "in", res = 150)
plot(var_clust, main = "Hierarchical Clustering of Environmental Variables",
     xlab = "", sub = "Based on correlation distance (1 - |r|)")
abline(h = 0.3, col = "red", lty = 2)  # Cut-off for r > 0.7
dev.off()

cat("\nPlots saved:\n")
cat("  - correlation_matrix_full.png\n")
cat("  - correlation_matrix_temperature.png\n")
cat("  - correlation_matrix_wind.png\n")
cat("  - correlation_matrix_current.png\n")
cat("  - temperature_vs_latitude.png\n")
cat("  - temperature_metrics_vs_latitude.png\n")
cat("  - temperature_by_site.png\n")
cat("  - pca_biplot.png\n")
cat("  - pca_scree.png\n")
cat("  - pca_loadings_by_type.png\n")
cat("  - wind_rose_simple.png\n")
cat("  - wind_rose_with_speed.png\n")
cat("  - current_rose.png\n")
cat("  - wind_uv_scatter.png\n")
cat("  - current_uv_scatter.png\n")
cat("  - variable_dendrogram.png\n")
cat("  - high_correlation_pairs.csv\n")