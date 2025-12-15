# Wind and Wave Energy Metrics - Standalone Functions
# Required packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, readr, lubridate, tidyr)

# =============================================================================
# Data Loading Functions (separate from wind rose functions)
# =============================================================================

#' Load eReefs data for metrics calculation
load_ereefs_for_metrics <- function(file_path) {
  df <- read_csv(file_path, 
                 show_col_types = FALSE,
                 na = c("", "NA", "no data", "No data", "NO DATA"))
  df$Date <- mdy(df$Date)
  
  # Convert numeric columns
  numeric_cols <- c("mean", "median", "p5", "p95", "lowest", "highest")
  existing_numeric_cols <- intersect(numeric_cols, names(df))
  
  for (col in existing_numeric_cols) {
    if (is.character(df[[col]])) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  
  return(df)
}

#' Prepare wind data for metrics (not rose diagrams)
prepare_wind_for_metrics <- function(df, statistic = "mean", sites_filter = NULL) {
  
  if (!is.null(sites_filter)) {
    df <- df %>% filter(Site %in% sites_filter)
  }
  
  wind_data <- df %>%
    filter(Variable %in% c("wspeed_u", "wspeed_v")) %>%
    select(Date, Variable, Site, Latitude, Longitude, all_of(statistic)) %>%
    rename(value = all_of(statistic)) %>%
    filter(!is.na(value), !is.na(Date))
  
  wind_wide <- wind_data %>%
    pivot_wider(names_from = Variable, values_from = value, names_prefix = "") %>%
    filter(!is.na(wspeed_u) & !is.na(wspeed_v)) %>%
    mutate(
      wspeed_u = as.numeric(as.character(wspeed_u)),
      wspeed_v = as.numeric(as.character(wspeed_v))
    ) %>%
    filter(!is.na(wspeed_u) & !is.na(wspeed_v)) %>%
    mutate(
      ws = sqrt(wspeed_u^2 + wspeed_v^2),
      wd = (270 - atan2(wspeed_v, wspeed_u) * 180/pi) %% 360,
      year = year(Date),
      month = month(Date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Summer",
        month %in% c(3, 4, 5) ~ "Autumn", 
        month %in% c(6, 7, 8) ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
      )
    ) %>%
    filter(ws > 0.1)
  
  return(wind_wide)
}

# =============================================================================
# Wind Pattern Metrics
# =============================================================================

#' Calculate comprehensive wind statistics
#' 
#' @param wind_data Prepared wind data
#' @return Data frame with wind pattern metrics
calculate_wind_metrics <- function(wind_data) {
  
  wind_stats <- wind_data %>%
    group_by(Site, Latitude, Longitude) %>%
    summarise(
      # Temporal coverage
      n_observations = n(),
      n_years = n_distinct(year),
      date_range = paste(min(Date), "to", max(Date)),
      
      # Basic wind statistics
      mean_wind_speed = mean(ws, na.rm = TRUE),
      median_wind_speed = median(ws, na.rm = TRUE),
      sd_wind_speed = sd(ws, na.rm = TRUE),
      p25_wind_speed = quantile(ws, 0.25, na.rm = TRUE),
      p75_wind_speed = quantile(ws, 0.75, na.rm = TRUE),
      p95_wind_speed = quantile(ws, 0.95, na.rm = TRUE),
      max_wind_speed = max(ws, na.rm = TRUE),
      min_wind_speed = min(ws, na.rm = TRUE),
      
      # Wind variability
      wind_speed_cv = sd(ws, na.rm = TRUE) / mean(ws, na.rm = TRUE),
      wind_speed_range = max(ws, na.rm = TRUE) - min(ws, na.rm = TRUE),
      
      # Directional statistics
      mean_u_component = mean(wspeed_u, na.rm = TRUE),
      mean_v_component = mean(wspeed_v, na.rm = TRUE),
      
      # Vector statistics
      vector_mean_speed = sqrt(mean_u_component^2 + mean_v_component^2),
      vector_mean_direction = (270 - atan2(mean_v_component, mean_u_component) * 180/pi) %% 360,
      
      # Directional constancy (0 = highly variable, 1 = perfectly constant)
      directional_constancy = vector_mean_speed / mean_wind_speed,
      
      # Directional components
      easterly_component = mean(pmax(0, wspeed_u), na.rm = TRUE),
      westerly_component = mean(pmax(0, -wspeed_u), na.rm = TRUE),  
      northerly_component = mean(pmax(0, wspeed_v), na.rm = TRUE),
      southerly_component = mean(pmax(0, -wspeed_v), na.rm = TRUE),
      
      # Trade wind alignment (ESE = 112.5°, typical for GBR)
      trade_wind_component = mean(
        wspeed_u * cos(112.5 * pi/180) + wspeed_v * sin(112.5 * pi/180), 
        na.rm = TRUE
      ),
      
      # Seasonal patterns
      summer_mean_speed = mean(ws[season == "Summer"], na.rm = TRUE),
      autumn_mean_speed = mean(ws[season == "Autumn"], na.rm = TRUE), 
      winter_mean_speed = mean(ws[season == "Winter"], na.rm = TRUE),
      spring_mean_speed = mean(ws[season == "Spring"], na.rm = TRUE),
      
      seasonal_speed_range = max(c(summer_mean_speed, autumn_mean_speed, 
                                   winter_mean_speed, spring_mean_speed), na.rm = TRUE) -
        min(c(summer_mean_speed, autumn_mean_speed,
              winter_mean_speed, spring_mean_speed), na.rm = TRUE),
      
      # Frequency statistics
      percent_calm_periods = sum(ws < 2, na.rm = TRUE) / n() * 100,
      percent_light_winds = sum(ws >= 2 & ws < 5, na.rm = TRUE) / n() * 100,
      percent_moderate_winds = sum(ws >= 5 & ws < 8, na.rm = TRUE) / n() * 100,
      percent_strong_winds = sum(ws >= 8 & ws < 12, na.rm = TRUE) / n() * 100,
      percent_very_strong_winds = sum(ws >= 12, na.rm = TRUE) / n() * 100,
      
      .groups = "drop"
    ) %>%
    
    # Add derived categorical variables
    mutate(
      # Wind exposure category
      wind_exposure_category = cut(mean_wind_speed,
                                   breaks = c(0, 3, 5, 7, 10, Inf),
                                   labels = c("Sheltered", "Moderate", "Exposed", "Very Exposed", "Extreme"),
                                   include.lowest = TRUE),
      
      # Directional consistency category
      directional_consistency_category = cut(directional_constancy,
                                             breaks = c(0, 0.3, 0.6, 0.8, 1),
                                             labels = c("Variable", "Moderate", "Consistent", "Very Consistent"),
                                             include.lowest = TRUE),
      
      # Predominant wind sector (8 directions)
      predominant_wind_sector = case_when(
        vector_mean_direction >= 337.5 | vector_mean_direction < 22.5 ~ "N",
        vector_mean_direction >= 22.5 & vector_mean_direction < 67.5 ~ "NE", 
        vector_mean_direction >= 67.5 & vector_mean_direction < 112.5 ~ "E",
        vector_mean_direction >= 112.5 & vector_mean_direction < 157.5 ~ "SE",
        vector_mean_direction >= 157.5 & vector_mean_direction < 202.5 ~ "S",
        vector_mean_direction >= 202.5 & vector_mean_direction < 247.5 ~ "SW",
        vector_mean_direction >= 247.5 & vector_mean_direction < 292.5 ~ "W",
        vector_mean_direction >= 292.5 & vector_mean_direction < 337.5 ~ "NW"
      ),
      
      # Seasonal variability category
      seasonal_variability = cut(seasonal_speed_range,
                                 breaks = c(0, 1, 2, 3, Inf),
                                 labels = c("Low", "Moderate", "High", "Very High"),
                                 include.lowest = TRUE)
    )
  
  return(wind_stats)
}

# =============================================================================
# Wave Energy Metrics
# =============================================================================

#' Calculate wave energy metrics from wind data
#' 
#' @param wind_data Prepared wind data
#' @return Data frame with wave energy metrics
calculate_wave_energy_metrics <- function(wind_data) {
  
  wave_stats <- wind_data %>%
    group_by(Site, Latitude, Longitude) %>%
    summarise(
      # Basic wave energy proxies (simplified without fetch calculations)
      wave_energy_proxy = mean(ws^3, na.rm = TRUE),  # Wave energy ∝ wind speed³
      wave_power_flux = mean(ws^2, na.rm = TRUE),    # Wave power flux ∝ wind speed²  
      wind_stress = mean(ws^2, na.rm = TRUE) * 1.2,  # Wind stress coefficient
      
      # Statistical wave energy metrics
      median_wave_energy = median(ws^3, na.rm = TRUE),
      p95_wave_energy = quantile(ws^3, 0.95, na.rm = TRUE),
      max_wave_energy = max(ws^3, na.rm = TRUE),
      
      # Directional wave exposure
      easterly_wave_exposure = mean(pmax(0, wspeed_u)^2, na.rm = TRUE),
      westerly_wave_exposure = mean(pmax(0, -wspeed_u)^2, na.rm = TRUE),
      northerly_wave_exposure = mean(pmax(0, wspeed_v)^2, na.rm = TRUE), 
      southerly_wave_exposure = mean(pmax(0, -wspeed_v)^2, na.rm = TRUE),
      
      # Total directional exposure
      total_directional_exposure = easterly_wave_exposure + westerly_wave_exposure +
        northerly_wave_exposure + southerly_wave_exposure,
      
      # Seasonal wave energy
      summer_wave_energy = mean(ws[season == "Summer"]^3, na.rm = TRUE),
      winter_wave_energy = mean(ws[season == "Winter"]^3, na.rm = TRUE),
      seasonal_wave_contrast = abs(summer_wave_energy - winter_wave_energy),
      
      # Extreme wave events
      extreme_wave_events = sum(ws > quantile(ws, 0.95, na.rm = TRUE), na.rm = TRUE),
      percent_extreme_events = extreme_wave_events / n() * 100,
      
      # Fetch-limited approximation (very simplified)
      # Assumes fetch limitation for wind speeds < 8 m/s
      fetch_limited_energy = mean(ifelse(ws < 8, ws^3 * 0.8, ws^3), na.rm = TRUE),
      fully_developed_energy = mean(ifelse(ws >= 8, ws^3, 0), na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    
    # Add derived metrics
    mutate(
      # Wave energy categories based on percentiles
      wave_energy_percentile = round(rank(wave_energy_proxy) / length(wave_energy_proxy) * 100),
      
      wave_exposure_category = cut(wave_energy_proxy,
                                   breaks = quantile(wave_energy_proxy, c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
                                   labels = c("Very Sheltered", "Sheltered", "Moderate", "Exposed", "Very Exposed"),
                                   include.lowest = TRUE),
      
      # Dominant exposure direction
      dominant_exposure_direction = case_when(
        easterly_wave_exposure > 0.5 * total_directional_exposure ~ "Easterly Dominated",
        southerly_wave_exposure > 0.5 * total_directional_exposure ~ "Southerly Dominated",
        westerly_wave_exposure > 0.5 * total_directional_exposure ~ "Westerly Dominated",
        northerly_wave_exposure > 0.5 * total_directional_exposure ~ "Northerly Dominated",
        TRUE ~ "Mixed Exposure"
      ),
      
      # Energy variability
      wave_energy_variability = cut(seasonal_wave_contrast,
                                    breaks = c(0, 50, 100, 200, Inf),
                                    labels = c("Low", "Moderate", "High", "Very High"),
                                    include.lowest = TRUE)
    )
  
  return(wave_stats)
}

# =============================================================================
# Combined Analysis Functions
# =============================================================================

#' Calculate comprehensive wind and wave metrics
#' 
#' @param file_path Path to eReefs CSV file
#' @param site_ids Optional vector of site IDs to analyze
#' @param statistic Which statistic to use ("mean", "median", "p95", etc.)
#' @return List containing wind_metrics and wave_metrics data frames
analyze_wind_wave_metrics <- function(file_path, site_ids = NULL, statistic = "mean") {
  
  # Load and prepare data
  raw_data <- load_ereefs_for_metrics(file_path)
  wind_data <- prepare_wind_for_metrics(raw_data, statistic = statistic, sites_filter = site_ids)
  
  if (nrow(wind_data) == 0) {
    stop("No valid wind data found for analysis")
  }
  
  # Calculate metrics
  wind_metrics <- calculate_wind_metrics(wind_data)
  wave_metrics <- calculate_wave_energy_metrics(wind_data)
  
  return(list(
    wind_metrics = wind_metrics,
    wave_metrics = wave_metrics,
    n_sites = length(unique(wind_data$Site)),
    date_range = paste(range(wind_data$Date), collapse = " to "),
    statistic_used = statistic
  ))
}

#' Create summary table combining key wind and wave metrics
#' 
#' @param wind_metrics Wind metrics data frame
#' @param wave_metrics Wave metrics data frame
#' @return Combined summary data frame
create_summary_table <- function(wind_metrics, wave_metrics) {
  
  summary <- wind_metrics %>%
    left_join(wave_metrics, by = c("Site", "Latitude", "Longitude")) %>%
    select(
      Site, Latitude, Longitude, n_observations, n_years,
      
      # Key wind metrics
      mean_wind_speed, vector_mean_direction, predominant_wind_sector,
      directional_constancy, trade_wind_component, wind_exposure_category,
      
      # Key wave metrics  
      wave_energy_proxy, wave_exposure_category, dominant_exposure_direction,
      wave_energy_percentile,
      
      # Seasonal patterns
      seasonal_speed_range, seasonal_wave_contrast
    ) %>%
    mutate(
      # Round for presentation
      Latitude = round(abs(Latitude), 3),
      Longitude = round(Longitude, 3),
      mean_wind_speed = round(mean_wind_speed, 1),
      vector_mean_direction = round(vector_mean_direction, 0),
      directional_constancy = round(directional_constancy, 2),
      trade_wind_component = round(trade_wind_component, 1),
      wave_energy_proxy = round(wave_energy_proxy, 0),
      seasonal_speed_range = round(seasonal_speed_range, 1),
      seasonal_wave_contrast = round(seasonal_wave_contrast, 0)
    ) %>%
    arrange(Latitude)
  
  return(summary)
}

# =============================================================================
# Export Functions
# =============================================================================

#' Export wind and wave metrics to CSV files
#' 
#' @param file_path Path to eReefs data
#' @param site_ids Optional vector of sites to analyze  
#' @param statistic Which statistic to use
#' @param output_dir Output directory for CSV files
#' @param prefix Prefix for output file names
export_metrics <- function(file_path, site_ids = NULL, statistic = "mean", 
                           output_dir = ".", prefix = "ereefs_metrics") {
  
  # Calculate metrics
  results <- analyze_wind_wave_metrics(file_path, site_ids, statistic)
  summary_table <- create_summary_table(results$wind_metrics, results$wave_metrics)
  
  # Create output file names
  wind_file <- file.path(output_dir, paste0(prefix, "_wind_", statistic, ".csv"))
  wave_file <- file.path(output_dir, paste0(prefix, "_wave_", statistic, ".csv"))
  summary_file <- file.path(output_dir, paste0(prefix, "_summary_", statistic, ".csv"))
  
  # Export to CSV
  write_csv(results$wind_metrics, wind_file)
  write_csv(results$wave_metrics, wave_file)
  write_csv(summary_table, summary_file)
  
  cat("Exported metrics for", results$n_sites, "sites (", results$statistic_used, "statistic)\n")
  cat("Date range:", results$date_range, "\n")
  cat("Files created:\n")
  cat(" -", wind_file, "\n")
  cat(" -", wave_file, "\n") 
  cat(" -", summary_file, "\n")
  
  return(results)
}

# =============================================================================
# Example Usage
# =============================================================================

# Example function showing how to use the metrics
example_metrics_analysis <- function() {
  cat("Wind and Wave Metrics Analysis\n")
  cat("==============================\n\n")
  
  cat("Basic usage:\n")
  cat("results <- analyze_wind_wave_metrics('your_data.csv', statistic = 'mean')\n")
  cat("wind_stats <- results$wind_metrics\n")
  cat("wave_stats <- results$wave_metrics\n\n")
  
  cat("Export all metrics:\n")
  cat("export_metrics('your_data.csv', \n")
  cat("               site_ids = c('CBHE_BA1S', 'ONLI_FL1S'), \n")
  cat("               statistic = 'p95',\n")
  cat("               prefix = 'my_analysis')\n\n")
  
  cat("Available statistics: mean, median, p5, p95, lowest, highest\n")
}

# Run example if in interactive session
if (interactive()) {
  example_metrics_analysis()
}