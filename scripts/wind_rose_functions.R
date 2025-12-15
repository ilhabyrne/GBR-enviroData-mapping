# Clean eReefs Wind Rose Analysis - No verbose output
# Required packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, readr, lubridate, openair, tidyr)

# =============================================================================
# Diagnostic Functions
# =============================================================================

check_data_format <- function(file_path, n_rows = 10) {
  cat("=== Data Format Check ===\n")
  
  # Read first few lines as text to see raw format
  raw_lines <- readLines(file_path, n = 5)
  cat("First few raw lines:\n")
  for (i in 1:min(length(raw_lines), 3)) {
    cat("Line", i, ":", raw_lines[i], "\n")
  }
  
  # Read with read_csv
  df <- read_csv(file_path, n_max = n_rows, show_col_types = FALSE)
  
  cat("\nColumn names:\n")
  print(names(df))
  
  cat("\nColumn types:\n")
  print(sapply(df, class))
  
  cat("\nSample data:\n")
  print(df)
  
  # Check specifically for date and numeric issues
  if ("Date" %in% names(df)) {
    cat("\nDate column analysis:\n")
    cat("- NA count:", sum(is.na(df$Date)), "out of", nrow(df), "\n")
    cat("- Sample values:", paste(head(df$Date, 3), collapse = ", "), "\n")
  }
  
  numeric_cols <- c("mean", "median", "p5", "p95", "lowest", "highest")
  for (col in intersect(numeric_cols, names(df))) {
    cat("\n", col, "column:\n")
    cat("- Type:", class(df[[col]]), "\n")
    cat("- Sample:", paste(head(df[[col]], 3), collapse = ", "), "\n")
  }
  
  return(df)
}

# =============================================================================
# Data Loading and Preparation Functions  
# =============================================================================

load_ereefs_timeseries <- function(file_path) {
  # Read CSV with proper handling of "no data" entries
  df <- read_csv(file_path, 
                 show_col_types = FALSE,
                 na = c("", "NA", "no data", "No data", "NO DATA"))
  
  # Parse dates - now handling m/d/yyyy format
  df$Date <- mdy(df$Date)  # Changed from ymd_hms to mdy
  
  # Convert numeric columns (they should now be properly numeric since "no data" = NA)
  numeric_cols <- c("mean", "median", "p5", "p95", "lowest", "highest")
  existing_numeric_cols <- intersect(numeric_cols, names(df))
  
  # Double-check conversion for any remaining character issues
  for (col in existing_numeric_cols) {
    if (is.character(df[[col]])) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  
  return(df)
}

prepare_wind_data <- function(df, statistic = "mean", sites_filter = NULL) {
  
  # Filter sites if specified
  if (!is.null(sites_filter)) {
    df <- df %>% filter(Site %in% sites_filter)
  }
  
  # Filter for wind speed variables and selected statistic
  wind_data <- df %>%
    filter(Variable %in% c("wspeed_u", "wspeed_v")) %>%
    select(Date, Variable, Site, Latitude, Longitude, all_of(statistic)) %>%
    rename(value = all_of(statistic)) %>%
    filter(!is.na(value))
  
  # DEBUG: Check what dates we actually have
  cat("DEBUG - Unique dates found:", length(unique(wind_data$Date)), "\n")
  cat("Date range:", paste(range(wind_data$Date, na.rm = TRUE), collapse = " to "), "\n")
  cat("First few dates:", paste(head(unique(wind_data$Date)), collapse = ", "), "\n")
  
  # Pivot to wide format
  wind_wide <- wind_data %>%
    pivot_wider(
      names_from = Variable, 
      values_from = value,
      names_prefix = ""
    ) %>%
    filter(!is.na(wspeed_u) & !is.na(wspeed_v))
  
  cat("DEBUG - After pivot, rows:", nrow(wind_wide), "\n")
  cat("DEBUG - Sites with data:", paste(unique(wind_wide$Site), collapse = ", "), "\n")
  
  # CRITICAL FIX: Ensure numeric conversion
  wind_wide <- wind_wide %>%
    mutate(
      # Convert to numeric (handles character/factor data)
      wspeed_u = as.numeric(as.character(wspeed_u)),
      wspeed_v = as.numeric(as.character(wspeed_v))
    ) %>%
    # Remove any rows that became NA after numeric conversion
    filter(!is.na(wspeed_u) & !is.na(wspeed_v)) %>%
    mutate(
      # Now safe to do math operations
      ws = sqrt(wspeed_u^2 + wspeed_v^2),
      wd = (270 - atan2(wspeed_v, wspeed_u) * 180/pi) %% 360,
      year = year(Date),
      month = month(Date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Summer",
        month %in% c(3, 4, 5) ~ "Autumn", 
        month %in% c(6, 7, 8) ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
      ),
      date = Date
    ) %>%
    filter(ws > 0.1)  # Remove very low speeds
  
  cat("DEBUG - Final data: ", nrow(wind_wide), "rows\n")
  cat("DEBUG - Date range in final data:", paste(range(wind_wide$Date, na.rm = TRUE), collapse = " to "), "\n")
  cat("DEBUG - Years:", paste(unique(wind_wide$year), collapse = ", "), "\n")
  cat("DEBUG - Months:", paste(unique(wind_wide$month), collapse = ", "), "\n")
  
  return(wind_wide)
}

# =============================================================================
# Wind Rose Creation Functions
# =============================================================================

create_site_windrose <- function(wind_data, site_id, title = NULL, 
                                 year_range = NULL, months = NULL, seasons = NULL) {
  
  # Filter data for specific site
  site_data <- wind_data %>% filter(Site == site_id)
  
  if (nrow(site_data) == 0) {
    stop(paste("No data found for site:", site_id))
  }
  
  # Apply filters
  if (!is.null(year_range)) {
    site_data <- site_data %>% filter(year >= year_range[1] & year <= year_range[2])
  }
  if (!is.null(months)) {
    site_data <- site_data %>% filter(month %in% months)
  }
  if (!is.null(seasons)) {
    site_data <- site_data %>% filter(season %in% seasons)
  }
  
  # Create title if not provided
  if (is.null(title)) {
    lat <- round(abs(site_data$Latitude[1]), 2)
    lon <- round(site_data$Longitude[1], 2)
    n_years <- length(unique(site_data$year))
    n_months <- nrow(site_data)
    title <- paste0("Wind Rose: ", site_id, " (", n_months, " months, ", n_years, " years)\n", 
                    lat, "°S, ", lon, "°E")
  }
  
  # Set wind speed breaks
  ws_range <- range(site_data$ws, na.rm = TRUE)
  if (max(ws_range) < 12) {
    breaks <- c(0, 2, 4, 6, 8, 10, Inf)
  } else {
    breaks <- c(0, 3, 6, 9, 12, 15, Inf)
  }
  
  # Create wind rose
  windRose(site_data,
           ws = "ws", wd = "wd", 
           main = title,
           cols = c("#E0F7FA", "#80DEEA", "#4DD0E1", "#26C6DA", "#00BCD4", "#006064"),
           breaks = breaks,
           angle = 15,
           key.position = "right",
           key.header = "Wind Speed\n(m/s)",
           key.footer = paste("n =", nrow(site_data), "months"),
           paddle = FALSE)
}

create_seasonal_windrose <- function(wind_data, site_id) {
  site_data <- wind_data %>% filter(Site == site_id)
  if (nrow(site_data) == 0) stop(paste("No data found for site:", site_id))
  
  lat <- round(abs(site_data$Latitude[1]), 2)
  lon <- round(site_data$Longitude[1], 2)
  
  windRose(site_data,
           ws = "ws", wd = "wd", type = "season",
           main = paste0("Seasonal Wind Roses: ", site_id, "\n(", lat, "°S, ", lon, "°E)"),
           cols = c("#E0F7FA", "#80DEEA", "#4DD0E1", "#26C6DA", "#00BCD4", "#006064"),
           breaks = c(0, 2, 4, 6, 8, 10, Inf),
           layout = c(2, 2), angle = 15)
}

create_multisite_windrose <- function(wind_data, site_ids, ncol = 3) {
  multi_site_data <- wind_data %>% 
    filter(Site %in% site_ids) %>%
    mutate(Site = factor(Site, levels = site_ids))
  
  windRose(multi_site_data,
           ws = "ws", wd = "wd", type = "Site",
           main = "Wind Rose Comparison (Monthly Data)",
           cols = c("#E0F7FA", "#80DEEA", "#4DD0E1", "#26C6DA", "#00BCD4", "#006064"),
           breaks = c(0, 2, 4, 6, 8, 10, Inf),
           layout = c(ncol, ceiling(length(site_ids)/ncol)),
           angle = 15)
}

# =============================================================================
# Statistical Functions
# =============================================================================

calculate_wind_statistics <- function(wind_data) {
  wind_stats <- wind_data %>%
    group_by(Site, Latitude, Longitude) %>%
    summarise(
      n_months = n(),
      n_years = n_distinct(year),
      mean_wind_speed = mean(ws, na.rm = TRUE),
      median_wind_speed = median(ws, na.rm = TRUE),
      max_wind_speed = max(ws, na.rm = TRUE),
      p95_wind_speed = quantile(ws, 0.95, na.rm = TRUE),
      mean_u = mean(wspeed_u, na.rm = TRUE),
      mean_v = mean(wspeed_v, na.rm = TRUE),
      vector_mean_speed = sqrt(mean_u^2 + mean_v^2),
      vector_mean_direction = (270 - atan2(mean_v, mean_u) * 180/pi) %% 360,
      directional_constancy = vector_mean_speed / mean_wind_speed,
      predominant_direction = {
        dir_sectors <- cut(wd, breaks = seq(0, 360, 30), 
                           labels = c("N", "NNE", "NE", "ENE", "E", "ESE", 
                                      "SE", "SSE", "S", "SSW", "SW", "WSW"),
                           include.lowest = TRUE)
        names(sort(table(dir_sectors), decreasing = TRUE))[1]
      },
      monthly_speed_cv = sd(ws, na.rm = TRUE) / mean(ws, na.rm = TRUE),
      percent_low_wind_months = sum(ws < 3, na.rm = TRUE) / n() * 100,
      percent_strong_wind_months = sum(ws > 8, na.rm = TRUE) / n() * 100,
      .groups = "drop"
    ) %>%
    mutate(
      exposure_category = cut(mean_wind_speed,
                              breaks = c(0, 3, 5, 7, 10, Inf),
                              labels = c("Sheltered", "Moderate", "Exposed", "Very Exposed", "Extreme"),
                              include.lowest = TRUE),
      directional_consistency = cut(directional_constancy,
                                    breaks = c(0, 0.3, 0.6, 0.8, 1),
                                    labels = c("Variable", "Moderate", "Consistent", "Very Consistent"),
                                    include.lowest = TRUE)
    )
  
  return(wind_stats)
}

create_wind_summary_table <- function(wind_stats) {
  summary_table <- wind_stats %>%
    select(Site, Latitude, Longitude, n_months, n_years,
           mean_wind_speed, vector_mean_direction, 
           predominant_direction, directional_constancy, exposure_category) %>%
    mutate(
      Latitude = round(abs(Latitude), 2),
      Longitude = round(Longitude, 2),
      mean_wind_speed = round(mean_wind_speed, 1),
      vector_mean_direction = round(vector_mean_direction, 0),
      directional_constancy = round(directional_constancy, 2)
    ) %>%
    arrange(Latitude)
  
  return(summary_table)
}

# =============================================================================
# Main Analysis Functions
# =============================================================================

run_wind_analysis <- function(file_path, site_ids = NULL, statistic = "mean", max_sites = 6) {
  
  # Load and prepare data
  raw_data <- load_ereefs_timeseries(file_path)
  available_sites <- unique(raw_data$Site)
  
  cat("DEBUG - Total raw data rows:", nrow(raw_data), "\n")
  cat("DEBUG - Available sites:", length(available_sites), "\n")
  cat("DEBUG - First few sites:", paste(head(available_sites, 5), collapse = ", "), "\n")
  
  # Select sites
  if (is.null(site_ids)) {
    site_ids <- available_sites[1:min(max_sites, length(available_sites))]
  } else {
    missing_sites <- setdiff(site_ids, available_sites)
    if (length(missing_sites) > 0) {
      cat("WARNING - Sites not found:", paste(missing_sites, collapse = ", "), "\n")
      site_ids <- intersect(site_ids, available_sites)
    }
  }
  
  cat("DEBUG - Selected sites:", paste(site_ids, collapse = ", "), "\n")
  
  # Check if selected sites have wind data
  for (site in site_ids) {
    site_wind_data <- raw_data %>% 
      filter(Site == site, Variable %in% c("wspeed_u", "wspeed_v"))
    cat("DEBUG -", site, "has", nrow(site_wind_data), "wind observations\n")
  }
  
  # Prepare wind data
  wind_data <- prepare_wind_data(raw_data, statistic = statistic, sites_filter = site_ids)
  
  if (nrow(wind_data) == 0) {
    cat("ERROR - No wind data survived processing!\n")
    return(NULL)
  }
  
  # Calculate statistics
  wind_stats <- calculate_wind_statistics(wind_data)
  summary_table <- create_wind_summary_table(wind_stats)
  print(summary_table)
  
  # Create wind roses
  cat("DEBUG - Creating wind roses...\n")
  for (site in site_ids) {
    site_data <- wind_data %>% filter(Site == site)
    if (nrow(site_data) > 0) {
      cat("DEBUG - Creating wind rose for", site, "with", nrow(site_data), "observations\n")
      try(create_site_windrose(wind_data, site), silent = FALSE)
    } else {
      cat("DEBUG - No data for", site, "after processing\n")
    }
  }
  
  # Multi-site comparison
  if (length(site_ids) > 1) {
    try(create_multisite_windrose(wind_data, site_ids[1:min(4, length(site_ids))]), silent = FALSE)
  }
  
  # Seasonal analysis
  first_site_data <- wind_data %>% filter(Site == site_ids[1])
  if (nrow(first_site_data) >= 12) {
    try(create_seasonal_windrose(wind_data, site_ids[1]), silent = FALSE)
  }
  
  return(list(
    wind_data = wind_data,
    wind_stats = wind_stats,
    summary_table = summary_table
  ))
}

quick_wind_analysis <- function(file_path, site_id, statistic = "mean", years = NULL) {
  raw_data <- load_ereefs_timeseries(file_path)
  wind_data <- prepare_wind_data(raw_data, statistic = statistic, sites_filter = site_id)
  create_site_windrose(wind_data, site_id, year_range = years)
  return(wind_data %>% filter(Site == site_id))
}

# =============================================================================
# Usage
# =============================================================================

# Example usage (uncomment to run):
# results <- run_wind_analysis("eReefs_hydro_monthly_wind.csv", 
#                              site_ids = c("ONLI_FL1S", "CBHE_BA1S", "TSMA_FR1S"),
#                              statistic = "mean")

# quick_wind_analysis("eReefs_hydro_monthly_wind.csv", "ONLI_FL1S", statistic = "p95")