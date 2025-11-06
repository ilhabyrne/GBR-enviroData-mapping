# Multi-File eReefs Processing - SIMPLIFIED FOR PRE-AGGREGATED DATA
# Data structure: One value per site per month (already temporally aggregated)
# Variables: u, v, wspeed_u, wspeed_v only (no waves)
# =============================================================================

if (!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, readr, lubridate, tidyr, purrr, stringr)

# =============================================================================
# Data Loading
# =============================================================================

#' Load multiple eReefs files (already monthly aggregated)
load_ereefs_files <- function(file_list) {
  
  cat("Loading", length(file_list), "eReefs files...\n")
  
  all_data <- map2_dfr(file_list, names(file_list), function(file, file_name) {
    
    cat("  Loading:", file_name, "\n")
    
    df <- read_csv(file, 
                   show_col_types = FALSE,
                   na = c("", "NA", "no data", "No data", "NO DATA"))
    
    # Parse date - handle different formats (ISO format: YYYY-MM-DDTHH:MM)
    if ("Aggregated Date/Time" %in% names(df)) {
      df$Date <- ymd_hms(df$`Aggregated Date/Time`)
      if (all(is.na(df$Date))) df$Date <- ymd(df$`Aggregated Date/Time`)
    } else if ("Date" %in% names(df)) {
      df$Date <- ymd_hms(df$Date)
      if (all(is.na(df$Date))) df$Date <- ymd(df$Date)
    }
    
    # Convert to Date class (remove time component)
    df$Date <- as.Date(df$Date)
    
    # Convert numeric columns
    numeric_cols <- c("mean", "median", "p5", "p95", "lowest", "highest")
    existing_numeric_cols <- intersect(numeric_cols, names(df))
    
    for (col in existing_numeric_cols) {
      if (is.character(df[[col]])) {
        df[[col]] <- as.numeric(df[[col]])
      }
    }
    
    df$file_source <- file_name
    
    return(df)
  })
  
  # Standardize column names
  all_data <- all_data %>%
    rename(
      Site = any_of(c("Site Name", "Site", "site")),
      Variable = any_of(c("Variable", "variable"))
    )
  
  cat("\nData loaded:\n")
  cat("  Total records:", nrow(all_data), "\n")
  cat("  Sites:", n_distinct(all_data$Site), "\n")
  cat("  Variables:", paste(unique(all_data$Variable), collapse = ", "), "\n")
  cat("  Date range:", format(min(all_data$Date, na.rm = TRUE)), "to",
      format(max(all_data$Date, na.rm = TRUE)), "\n")
  cat("  Unique months:", n_distinct(all_data$Date), "\n\n")
  
  return(all_data)
}

# =============================================================================
# Add Temporal Variables
# =============================================================================

#' Add year, month, season variables
add_temporal_variables <- function(data) {
  data %>%
    mutate(
      year = year(Date),
      month = month(Date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Summer",
        month %in% c(3, 4, 5) ~ "Autumn",
        month %in% c(6, 7, 8) ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
      ),
      year_season = paste(year, season, sep = "_")
    )
}

# =============================================================================
# Site-Level Statistics from Monthly Data
# =============================================================================

#' Calculate site-level statistics
#' Since data is already monthly, just aggregate across months
#' 
#' @param data Combined monthly data
#' @param statistic Which column to use ("mean", "median", "p5", "p95")
#' @return Site-level statistics
calculate_site_statistics <- function(data, statistic = "mean") {
  
  data <- add_temporal_variables(data)
  
  # Select the statistic to work with
  data <- data %>%
    rename(value = all_of(statistic))
  
  # Calculate site-level statistics from monthly values
  site_stats <- data %>%
    group_by(Site, Variable, Latitude, Longitude) %>%
    summarise(
      # Temporal statistics
      overall_mean = mean(value, na.rm = TRUE),
      overall_median = median(value, na.rm = TRUE),
      overall_sd = sd(value, na.rm = TRUE),
      overall_cv = sd(value, na.rm = TRUE) / abs(mean(value, na.rm = TRUE)),
      
      # Range
      overall_min = min(value, na.rm = TRUE),
      overall_max = max(value, na.rm = TRUE),
      overall_range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
      
      # Percentiles
      overall_p25 = quantile(value, 0.25, na.rm = TRUE),
      overall_p75 = quantile(value, 0.75, na.rm = TRUE),
      overall_p95 = quantile(value, 0.95, na.rm = TRUE),
      
      # Sample info
      n_months = n(),
      n_years = n_distinct(year),
      date_start = min(Date),
      date_end = max(Date),
      
      .groups = "drop"
    ) %>%
    mutate(
      variable_type = case_when(
        Variable %in% c("wspeed_u", "wspeed_v") ~ "wind",
        Variable %in% c("u", "v") ~ "current",
        TRUE ~ "other"
      )
    )
  
  cat("Site statistics calculated:\n")
  cat("  Sites:", n_distinct(site_stats$Site), "\n")
  cat("  Mean months per site:", round(mean(site_stats$n_months)), "\n\n")
  
  return(site_stats)
}

#' Calculate seasonal patterns
calculate_seasonal_patterns <- function(data, statistic = "mean") {
  
  data <- add_temporal_variables(data)
  
  data <- data %>%
    rename(value = all_of(statistic))
  
  seasonal_stats <- data %>%
    group_by(Site, Variable, Latitude, Longitude, season) %>%
    summarise(
      seasonal_mean = mean(value, na.rm = TRUE),
      seasonal_sd = sd(value, na.rm = TRUE),
      seasonal_n = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = season,
      values_from = c(seasonal_mean, seasonal_sd, seasonal_n),
      names_sep = "_"
    ) %>%
    mutate(
      # Seasonal range
      seasonal_range = pmax(
        seasonal_mean_Summer, seasonal_mean_Autumn,
        seasonal_mean_Winter, seasonal_mean_Spring, na.rm = TRUE
      ) - pmin(
        seasonal_mean_Summer, seasonal_mean_Autumn,
        seasonal_mean_Winter, seasonal_mean_Spring, na.rm = TRUE
      ),
      
      # Dominant season
      dominant_season = case_when(
        seasonal_mean_Summer == pmax(seasonal_mean_Summer, seasonal_mean_Autumn,
                                     seasonal_mean_Winter, seasonal_mean_Spring,
                                     na.rm = TRUE) ~ "Summer",
        seasonal_mean_Winter == pmax(seasonal_mean_Summer, seasonal_mean_Autumn,
                                     seasonal_mean_Winter, seasonal_mean_Spring,
                                     na.rm = TRUE) ~ "Winter",
        seasonal_mean_Autumn == pmax(seasonal_mean_Summer, seasonal_mean_Autumn,
                                     seasonal_mean_Winter, seasonal_mean_Spring,
                                     na.rm = TRUE) ~ "Autumn",
        TRUE ~ "Spring"
      )
    )
  
  return(seasonal_stats)
}

# =============================================================================
# Vector Statistics for Currents/Winds
# =============================================================================

#' Calculate vector statistics from u/v components
#' Data is already monthly, so we average across months
calculate_vector_statistics <- function(data, statistic = "mean") {
  
  data <- add_temporal_variables(data)
  
  data <- data %>%
    rename(value = all_of(statistic))
  
  # Process winds (wspeed_u, wspeed_v)
  wind_data <- data %>%
    filter(Variable %in% c("wspeed_u", "wspeed_v")) %>%
    select(Site, Latitude, Longitude, Date, year, month, season, Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    filter(!is.na(wspeed_u) & !is.na(wspeed_v))
  
  if (nrow(wind_data) > 0) {
    wind_vectors <- wind_data %>%
      group_by(Site, Latitude, Longitude) %>%
      summarise(
        # Mean components
        mean_u = mean(wspeed_u, na.rm = TRUE),
        mean_v = mean(wspeed_v, na.rm = TRUE),
        
        # Mean scalar speed (average of magnitudes)
        mean_speed = mean(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # SD of speed (temporal variability)
        sd_speed = sd(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # Vector magnitude (resultant)
        vector_magnitude = sqrt(mean_u^2 + mean_v^2),
        
        # Vector direction (0 = North, clockwise)
        vector_direction = (270 - atan2(mean_v, mean_u) * 180/pi) %% 360,
        
        # Directional constancy (0 = variable, 1 = constant)
        directional_constancy = sqrt(mean_u^2 + mean_v^2) / 
                               mean(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # Predominant direction sector
        predominant_sector = case_when(
          vector_direction >= 337.5 | vector_direction < 22.5 ~ "N",
          vector_direction >= 22.5 & vector_direction < 67.5 ~ "NE",
          vector_direction >= 67.5 & vector_direction < 112.5 ~ "E",
          vector_direction >= 112.5 & vector_direction < 157.5 ~ "SE",
          vector_direction >= 157.5 & vector_direction < 202.5 ~ "S",
          vector_direction >= 202.5 & vector_direction < 247.5 ~ "SW",
          vector_direction >= 247.5 & vector_direction < 292.5 ~ "W",
          TRUE ~ "NW"
        ),
        
        # Trade wind component (SE = 135° for GBR)
        trade_wind_component = mean_u * cos(135 * pi/180) + mean_v * sin(135 * pi/180),
        
        # Alongshore/crossshore components
        alongshore_component = mean_v,
        crossshore_component = mean_u,
        
        # Seasonal wind variation
        summer_mean_speed = mean(sqrt(wspeed_u[season == "Summer"]^2 + 
                                     wspeed_v[season == "Summer"]^2), na.rm = TRUE),
        winter_mean_speed = mean(sqrt(wspeed_u[season == "Winter"]^2 + 
                                     wspeed_v[season == "Winter"]^2), na.rm = TRUE),
        
        # Sample size
        n_months = n(),
        
        .groups = "drop"
      ) %>%
      mutate(
        variable_type = "wind",
        seasonal_wind_diff = abs(summer_mean_speed - winter_mean_speed)
      )
  } else {
    wind_vectors <- NULL
  }
  
  # Process currents (u, v)
  current_data <- data %>%
    filter(Variable %in% c("u", "v")) %>%
    select(Site, Latitude, Longitude, Date, year, month, season, Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    filter(!is.na(u) & !is.na(v))
  
  if (nrow(current_data) > 0) {
    current_vectors <- current_data %>%
      group_by(Site, Latitude, Longitude) %>%
      summarise(
        mean_u = mean(u, na.rm = TRUE),
        mean_v = mean(v, na.rm = TRUE),
        mean_speed = mean(sqrt(u^2 + v^2), na.rm = TRUE),
        sd_speed = sd(sqrt(u^2 + v^2), na.rm = TRUE),
        vector_magnitude = sqrt(mean_u^2 + mean_v^2),
        vector_direction = (270 - atan2(mean_v, mean_u) * 180/pi) %% 360,
        directional_constancy = sqrt(mean_u^2 + mean_v^2) / 
                               mean(sqrt(u^2 + v^2), na.rm = TRUE),
        predominant_sector = case_when(
          vector_direction >= 337.5 | vector_direction < 22.5 ~ "N",
          vector_direction >= 22.5 & vector_direction < 67.5 ~ "NE",
          vector_direction >= 67.5 & vector_direction < 112.5 ~ "E",
          vector_direction >= 112.5 & vector_direction < 157.5 ~ "SE",
          vector_direction >= 157.5 & vector_direction < 202.5 ~ "S",
          vector_direction >= 202.5 & vector_direction < 247.5 ~ "SW",
          vector_direction >= 247.5 & vector_direction < 292.5 ~ "W",
          TRUE ~ "NW"
        ),
        trade_wind_component = mean_u * cos(135 * pi/180) + mean_v * sin(135 * pi/180),
        alongshore_component = mean_v,
        crossshore_component = mean_u,
        summer_mean_speed = mean(sqrt(u[season == "Summer"]^2 + 
                                     v[season == "Summer"]^2), na.rm = TRUE),
        winter_mean_speed = mean(sqrt(u[season == "Winter"]^2 + 
                                     v[season == "Winter"]^2), na.rm = TRUE),
        n_months = n(),
        .groups = "drop"
      ) %>%
      mutate(
        variable_type = "current",
        seasonal_current_diff = abs(summer_mean_speed - winter_mean_speed)
      )
  } else {
    current_vectors <- NULL
  }
  
  # Combine
  all_vectors <- bind_rows(wind_vectors, current_vectors)
  
  cat("Vector statistics calculated:\n")
  if (!is.null(wind_vectors)) {
    cat("  Wind sites:", nrow(wind_vectors), "\n")
  }
  if (!is.null(current_vectors)) {
    cat("  Current sites:", nrow(current_vectors), "\n")
  }
  cat("\n")
  
  return(all_vectors)
}

# =============================================================================
# Master Processing Function
# =============================================================================

#' Process all eReefs files (already monthly aggregated)
#' 
#' @param file_paths Named vector of file paths
#' @param statistic Which column to use ("mean", "median", "p5", "p95")
#' @param coral_sites Optional vector of coral site names
#' @return List with all processed data
process_ereefs_monthly_data <- function(file_paths,
                                        statistic = "mean",
                                        coral_sites = NULL) {
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("PROCESSING PRE-AGGREGATED MONTHLY eREEFS DATA\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Step 1: Load all files
  cat("Step 1: Loading files...\n")
  monthly_data <- load_ereefs_files(file_paths)
  
  # Step 2: Calculate site-level statistics
  cat("Step 2: Calculating site-level statistics...\n")
  site_stats <- calculate_site_statistics(monthly_data, statistic)
  
  # Step 3: Calculate seasonal patterns
  cat("Step 3: Calculating seasonal patterns...\n")
  seasonal_patterns <- calculate_seasonal_patterns(monthly_data, statistic)
  
  # Step 4: Calculate vector statistics
  cat("Step 4: Calculating vector statistics...\n")
  vector_stats <- calculate_vector_statistics(monthly_data, statistic)
  
  # Step 5: Create wide format tables
  cat("Step 5: Creating wide format summaries...\n")
  
  wind_wide <- site_stats %>%
    filter(variable_type == "wind") %>%
    select(Site, Variable, Latitude, Longitude,
           overall_mean, overall_sd, overall_cv, overall_p95, n_months) %>%
    pivot_wider(
      names_from = Variable,
      values_from = c(overall_mean, overall_sd, overall_cv, overall_p95),
      names_sep = "_"
    )
  
  current_wide <- site_stats %>%
    filter(variable_type == "current") %>%
    select(Site, Variable, Latitude, Longitude,
           overall_mean, overall_sd, overall_cv, overall_p95, n_months) %>%
    pivot_wider(
      names_from = Variable,
      values_from = c(overall_mean, overall_sd, overall_cv, overall_p95),
      names_sep = "_"
    )
  
  # Step 6: Create coral site environmental matrix
  if (!is.null(coral_sites)) {
    cat("Step 6: Creating environmental matrix for coral sites...\n")
    
    coral_matrix <- vector_stats %>%
      filter(variable_type == "wind") %>%
      filter(Site %in% coral_sites) %>%
      left_join(
        seasonal_patterns %>%
          filter(Variable == "wspeed_u") %>%
          select(Site, seasonal_range, dominant_season),
        by = "Site"
      ) %>%
      mutate(
        # Composite indices
        wind_exposure_index = mean_speed,
        temporal_variability = sd_speed,
        transport_potential = directional_constancy * mean_speed,
        
        # Categories
        exposure_category = cut(
          mean_speed,
          breaks = quantile(mean_speed, c(0, 0.33, 0.67, 1), na.rm = TRUE),
          labels = c("Sheltered", "Moderate", "Exposed"),
          include.lowest = TRUE
        ),
        
        directional_category = cut(
          directional_constancy,
          breaks = c(0, 0.4, 0.7, 1),
          labels = c("Variable", "Moderate", "Consistent"),
          include.lowest = TRUE
        )
      ) %>%
      rename(coral_site = Site)
  } else {
    coral_matrix <- NULL
  }
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("PROCESSING COMPLETE\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  return(list(
    monthly_data = monthly_data,
    site_statistics = site_stats,
    seasonal_patterns = seasonal_patterns,
    vector_statistics = vector_stats,
    wind_wide = wind_wide,
    current_wide = current_wide,
    coral_environmental_matrix = coral_matrix,
    processing_info = list(
      n_files = length(file_paths),
      statistic_used = statistic,
      n_sites = n_distinct(monthly_data$Site),
      n_months = n_distinct(monthly_data$Date),
      date_range = range(monthly_data$Date, na.rm = TRUE),
      processing_date = Sys.time()
    )
  ))
}

# =============================================================================
# Export Functions
# =============================================================================

#' Export all results
export_results <- function(results, output_dir = ".", prefix = "ereefs") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Monthly data
  write_csv(results$monthly_data,
            file.path(output_dir, paste0(prefix, "_monthly_data.csv")))
  
  # Site statistics
  write_csv(results$site_statistics,
            file.path(output_dir, paste0(prefix, "_site_statistics.csv")))
  
  # Seasonal patterns
  write_csv(results$seasonal_patterns,
            file.path(output_dir, paste0(prefix, "_seasonal_patterns.csv")))
  
  # Vector statistics
  write_csv(results$vector_statistics,
            file.path(output_dir, paste0(prefix, "_vector_statistics.csv")))
  
  # Wide format
  write_csv(results$wind_wide,
            file.path(output_dir, paste0(prefix, "_wind_wide.csv")))
  
  if (nrow(results$current_wide) > 0) {
    write_csv(results$current_wide,
              file.path(output_dir, paste0(prefix, "_current_wide.csv")))
  }
  
  # Coral environmental matrix (KEY OUTPUT)
  if (!is.null(results$coral_environmental_matrix)) {
    write_csv(results$coral_environmental_matrix,
              file.path(output_dir, paste0(prefix, "_environmental_matrix.csv")))
  }
  
  cat("\nResults exported to:", output_dir, "\n")
  cat("Files created:\n")
  cat("  ✓ Monthly data\n")
  cat("  ✓ Site statistics\n")
  cat("  ✓ Seasonal patterns\n")
  cat("  ✓ Vector statistics\n")
  if (!is.null(results$coral_environmental_matrix)) {
    cat("  ✓ Environmental matrix for genomics (KEY FILE)\n")
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Simplified eReefs Processing Functions Loaded\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nFor pre-aggregated monthly data (one value per site per month)\n")
cat("Main function: process_ereefs_monthly_data()\n")
cat("Export function: export_results()\n\n")
