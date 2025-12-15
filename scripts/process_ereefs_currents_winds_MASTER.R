# Aggregate raw eReefs data & extract summary statistics for CURRENTS & WIND
# Uses monthly eReefs GBR1 data
# Raw metrics include u, v, wspeed_u, & wspeed_v
# ============================================================================ #

# Load libraries
library(tidyverse)
library(lubridate)
library(stringdist)

# Setwd
setwd("~/Documents/GitHub/GBR-enviroData-mapping/02_eReefs")

# ============================================================================ #
# STEP 1: MERGE ALL FILES ----
# ============================================================================ #

# Merge multiple CSV files into one data frame
merge_ereefs_files <- function(file_paths) {
  
  cat("Merging", length(file_paths), "files...\n\n")
  
  all_data <- map_dfr(file_paths, function(file) {
    cat("  Reading:", basename(file), "\n")
    
    df <- read_csv(file, 
                   show_col_types = FALSE,
                   na = c("", "NA", "no data", "No data", "NO DATA"))
    
    return(df)
  })
  
  # Standardize column names
  all_data <- all_data %>%
    rename(
      Site = any_of(c("Site Name", "Site", "site")),
      Variable = any_of(c("Variable", "variable")),
      Date = any_of(c("Aggregated Date/Time", "Date", "date"))
    )
  
  # Parse date (handle ISO format with timestamp: 2014-12-01T00:00)
  all_data$Date <- ymd(all_data$Date)  # Changed from mdy() to ymd()
  
  # If that doesn't work, try parsing as datetime then extract date
  if (all(is.na(all_data$Date))) {
    all_data$Date <- as.Date(ymd_hms(all_data$Date))
  }
  
  # Convert numeric columns
  numeric_cols <- c("mean", "median", "p5", "p95", "lowest", "highest")
  for (col in intersect(numeric_cols, names(all_data))) {
    if (is.character(all_data[[col]])) {
      all_data[[col]] <- as.numeric(all_data[[col]])
    }
  }
  
  cat("\n✓ Files merged successfully\n")
  cat("  Total records:", nrow(all_data), "\n")
  cat("  Unique sites:", n_distinct(all_data$Site), "\n")
  cat("  Variables:", paste(unique(all_data$Variable), collapse = ", "), "\n")
  cat("  Date range:", format(min(all_data$Date)), "to", format(max(all_data$Date)), "\n")
  cat("  Unique dates:", n_distinct(all_data$Date), "\n\n")
  
  # Check for duplicates (including lat/long to distinguish different locations)
 
  # TRUE duplicates = same site name, variable, date, AND coords
  dups <- all_data %>%
    group_by(Site, Variable, Date, Latitude, Longitude) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # Sites with SAME NAME but DIFFERENT locations (should be kept separate!)
  same_name_diff_location <- all_data %>%
    group_by(Site, Variable) %>%
    summarise(
      n_unique_locations = n_distinct(paste(Latitude, Longitude)),
      .groups = "drop"
    ) %>%
    filter(n_unique_locations > 1)
  
  if (nrow(same_name_diff_location) > 0) {
    cat("WARNING: Found", nrow(same_name_diff_location), 
        "site-variable combinations with SAME NAME but DIFFERENT locations\n")
    cat("  These will be renamed to preserve as unique sites\n\n")
    
    cat("Examples of sites with same name, different locations:\n")
    print(head(same_name_diff_location, 10))
    cat("\n")
    
    # Create detailed report of sites with different locations
    location_report <- all_data %>%
      group_by(Site) %>%
      summarise(
        n_unique_locations = n_distinct(paste(Latitude, Longitude)),
        locations = paste(unique(paste0("(", round(Latitude, 4), ", ", 
                                        round(Longitude, 4), ")")), 
                          collapse = " | "),
        n_records = n(),
        .groups = "drop"
      ) %>%
      filter(n_unique_locations > 1) %>%
      arrange(desc(n_unique_locations))
    
    # Save detailed report
    write_csv(location_report, "duplicate_locations_report.csv")
    
    # Rename sites with different locations
    # First, add Site_original column for ALL data
    all_data <- all_data %>%
      mutate(Site_original = Site)
    
    # Then rename only the sites with multiple locations
    all_data <- all_data %>%
      group_by(Site_original) %>%
      arrange(Latitude, Longitude) %>%
      mutate(
        n_locations = n_distinct(paste(round(Latitude, 6), round(Longitude, 6)))
      ) %>%
      mutate(
        # Only add suffix if site has multiple locations
        Site = if (n_locations[1] > 1) {
          paste0(Site_original, "_v", dense_rank(paste(round(Latitude, 6), round(Longitude, 6))))
        } else {
          Site_original
        }
      ) %>%
      ungroup() %>%
      select(-n_locations)
    
    cat("Sites with different locations have been renamed\n")
    cat("Example: CBHE_LA2S at two locations → CBHE_LA2S_v1, CBHE_LA2S_v2\n\n")
    
    # Also create a mapping table showing old and new names
    location_mapping <- all_data %>%
      filter(Site_original %in% location_report$Site) %>%
      select(Site_original, Site, Latitude, Longitude, Depth) %>%
      distinct() %>%
      arrange(Site_original, Latitude, Longitude)
    
    write_csv(location_mapping, "site_renaming_mapping.csv")
  }
  
  # NOW check for true duplicates (after renaming)
  dups <- all_data %>%
    group_by(Site, Variable, Date) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(dups) > 0) {
    cat("Found", nrow(dups), "TRUE duplicate records (same site, variable, date, location)\n")
    
    # Show some examples
    cat("\nExample true duplicates:\n")
    example_dups <- dups %>%
      select(Site, Variable, Date, Depth, Latitude, Longitude, mean) %>%
      arrange(Site, Variable, Date) %>%
      head(10)
    print(example_dups)
    
    # Average true duplicates
    all_data <- all_data %>%
      group_by(Site, Variable, Date, Latitude, Longitude, Depth) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE),
                Site_original = first(Site_original),
                .groups = "drop")
  } else {
    cat("No true duplicates found\n\n")
  }
  
  return(all_data)
}

# ============================================================================ #
# THEN USE YOUR EXISTING FUNCTIONS
# ============================================================================ #

# Source functions for extracting metrics
source("../scripts/ereefs_currents_wind_function.R")
source("../scripts/ereefs_currents_wind_function_enhanced.R")

# ============================================================================ #
# CONFIGURATION ----
# ============================================================================ #

# Your file paths - EDIT
ereefs_files <- c(
  "raw/RRAP_eReefs_GBR1_0.5m_currents_wind.csv", 
  "raw/RRAP_eReefs_GBR1_2.35m_currents_wind.csv",
  "raw/RRAP_eReefs_GBR1_5.35m_currents_wind.csv",
  "raw/RRAP_eReefs_GBR1_9m_currents_wind.csv",
  "raw/RRAP_eReefs_GBR1_13m_currents_wind.csv"
)

# Your coral sites - replace with list of specific sites if desired
coral_sites <- NULL

# Output directory
output_dir <- "processed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

# ============================================================================ #
# STEP 2: MERGE FILES ----
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("STEP 1: MERGING FILES\n")
cat(rep("=", 70), "\n\n", sep = "")

merged_data <- merge_ereefs_files(ereefs_files)

# Save merged file for reference
write_csv(merged_data, file.path(output_dir, "RRAP_eReefs_GBR1_currents_wind_merged.csv"))
cat("✓ Merged data saved:", file.path(output_dir, "RRAP_eReefs_GBR1_currents_wind_merged.csv"), "\n\n")

# ============================================================================ #
# STEP 3: PROCESS MERGED DATA ----
# ============================================================================ #

# Add temporal variables
merged_data <- merged_data %>%
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

# Calculate site statistics (using "mean" column)
site_stats <- calculate_site_statistics(merged_data, statistic = "mean")

# Calculate seasonal patterns
seasonal_patterns <- calculate_seasonal_patterns(merged_data, statistic = "mean")

# Calculate vector statistics (with biological metrics!)
vector_stats <- calculate_vector_statistics_biological(merged_data, statistic = "mean")

# Also calculate with p95 for extreme conditions
vector_stats_p95 <- calculate_vector_statistics_biological(merged_data, statistic = "p95")

# ============================================================================ #
# STEP 4: MATCH TO CORAL SITES ----
# ============================================================================ #

# OPTIONAL - if you want to just look at specific sites

# Function to match sites
find_site_matches <- function(coral_sites, all_ereefs_sites) {
  
  matches <- tibble(
    coral_site = coral_sites,
    exact_match = NA_character_,
    fuzzy_match = NA_character_,
    fuzzy_distance = NA_integer_,
    matched_site = NA_character_,
    match_method = NA_character_
  )
  
  for (i in seq_along(coral_sites)) {
    cs <- coral_sites[i]
    
    # Exact match
    if (cs %in% all_ereefs_sites) {
      matches$exact_match[i] <- cs
      matches$matched_site[i] <- cs
      matches$match_method[i] <- "exact"
    } else {
      # Fuzzy match
      distances <- stringdist(cs, all_ereefs_sites, method = "lv")
      best_idx <- which.min(distances)
      best_dist <- distances[best_idx]
      
      if (best_dist <= 5) {
        matches$fuzzy_match[i] <- all_ereefs_sites[best_idx]
        matches$fuzzy_distance[i] <- best_dist
        matches$matched_site[i] <- all_ereefs_sites[best_idx]
        matches$match_method[i] <- "fuzzy"
      } else {
        # Substring match
        substring_matches <- grep(cs, all_ereefs_sites, ignore.case = TRUE, value = TRUE)
        if (length(substring_matches) > 0) {
          matches$matched_site[i] <- substring_matches[1]
          matches$match_method[i] <- "substring"
        } else {
          matches$match_method[i] <- "NO_MATCH"
        }
      }
    }
  }
  
  return(matches)
}

site_matching <- find_site_matches(coral_sites, all_sites$Site)
write_csv(site_matching, file.path(output_dir, "site_matching_table.csv"))
  
# ============================================================================ #
# STEP 5: EXPORT ALL RESULTS ----
# ============================================================================ #
  
  # Monthly data
  write_csv(merged_data, 
            file.path(output_dir, "currents_wind_merged.csv"))
  
  # Site statistics
  write_csv(site_stats,
            file.path(output_dir, "currents_wind_general_stats.csv"))
  
  # Seasonal patterns
  write_csv(seasonal_patterns,
            file.path(output_dir, "currents_wind_seasonal_patterns.csv"))
  
  # Vector statistics (all sites)
  write_csv(vector_stats,
            file.path(output_dir, "currents_wind_vector_stats.csv"))
  
  # Vector statistics (extreme)
  write_csv(vector_stats_p95,
            file.path(output_dir, "currents_wind_vector_stats_p95.csv"))
  
  
# ============================================================================ #
# STEP 6: CREATE COMPREHENSIVE MASTER SHEET ----
# ============================================================================ #
  
  # Merge all statistics into one comprehensive sheet
  
  # Start with vector statistics (has most metrics)
  master_wind <- vector_stats %>%
    filter(variable_type == "wind") %>%
    select(-variable_type)
  
  master_current <- vector_stats %>%
    filter(variable_type == "current") %>%
    select(-variable_type) %>%
    # Rename current columns to distinguish from wind
    rename_with(~paste0("current_", .), 
                .cols = c(mean_u:n_years, -Site, -Latitude, -Longitude))
  
  # Add seasonal patterns for wind
  seasonal_wind <- seasonal_patterns %>%
    filter(Variable == "wspeed_u") %>%
    select(-Variable)
  
  # Add extreme (p95) statistics for wind
  extreme_wind <- vector_stats_p95 %>%
    filter(variable_type == "wind") %>%
    select(Site, 
           extreme_mean_wind_speed = mean_wind_speed,
           extreme_vector_magnitude = vector_magnitude,
           extreme_wind_transport_potential = wind_transport_potential,
           extreme_wave_energy_proxy = mean_wave_energy_proxy) %>%
    distinct()
  
  # Merge everything together
  master_sheet <- master_wind %>%
    left_join(master_current, by = c("Site", "Latitude", "Longitude")) %>%
    left_join(seasonal_wind, by = c("Site", "Latitude", "Longitude")) %>%
    left_join(extreme_wind, by = "Site") %>%
    # Calculate some additional composite metrics
    mutate(
      # Wind-current interaction
      wind_current_ratio = if ("current_mean_current_speed" %in% names(.)) {
        mean_wind_speed / current_mean_current_speed
      } else {
        NA_real_
      },
      
      # Hydrodynamic exposure index
      hydrodynamic_exposure_index = scale(mean_wind_speed)[,1] + 
        scale(coalesce(mean_wave_energy_proxy, 0))[,1],
      
      # Environmental variability index
      environmental_variability_index = scale(cv_wind_speed)[,1] + 
        scale(coalesce(seasonal_wind_range, 0))[,1],
      
      # Connectivity index
      connectivity_index = scale(wind_transport_potential)[,1],
      
      # Retention index  
      retention_index = scale(wind_retention_potential)[,1]
    ) %>%
    # Arrange by latitude
    arrange(Latitude)
  
  # Export master sheet
  write_csv(master_sheet, 
            file.path(output_dir, "currents_wind_all_metrics_MASTER.csv"))
  
# ============================================================================ #
# STEP 8: VISUALIZATIONS ----
# ============================================================================ #
  
# OPTIONAL - corresponds to STEP 4
  
fig_dir <- file.path(output_dir, "figures")
  
  ## Plot 1: Wind speed gradient ----
  p1 <- ggplot(env_matrix,
               aes(x = abs(Latitude), y = mean_wind_speed,
                   color = wind_exposure_category, label = coral_site)) +
    geom_point(size = 4) +
    geom_text(hjust = -0.2, size = 3, show.legend = FALSE) +
    geom_smooth(method = "loess", se = TRUE, show.legend = FALSE) +
    labs(x = "Latitude (°S)", y = "Mean wind speed (m/s)") +
    theme_minimal()
  
  ggsave(file.path(fig_dir, "wind_gradient.png"), p1, 
         width = 10, height = 6, dpi = 300)
  
  ## Plot 2: Transport potential vs retention ----
  p2 <- ggplot(env_matrix,
               aes(x = wind_transport_potential, y = wind_retention_potential,
                   label = coral_site)) +
    geom_point(size = 4, alpha = 0.6) +
    geom_text(hjust = -0.2, size = 3) +
    labs(x = "Transport potential (high = long dispersal)",
         y = "Retention potential (high = local retention)") +
    theme_minimal()
  
  ggsave(file.path(fig_dir, "transport_vs_retention.png"), p2,
         width = 9, height = 7, dpi = 300)
  
  ## Plot 3: Directional patterns ----
  p3 <- ggplot(env_matrix,
               aes(x = abs(Latitude), y = directional_constancy,
                   color = predominant_wind_sector, size = mean_wind_speed)) +
    geom_point(alpha = 0.7) +
    geom_text(aes(label = coral_site), hjust = -0.2, size = 3,
              show.legend = FALSE) +
    theme_minimal()
  
  ggsave(file.path(fig_dir, "directional_patterns.png"), p3,
         width = 11, height = 6, dpi = 300)