# Coral Environmental Data Aggregation for Seriatopora hystrix Population Genomics
# Simplified version for pre-filtered depth data (-9m and -2.35m)

# Load required libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(stringr)

# Function to calculate coefficient of variation
cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# Function to calculate Maximum Monthly Mean (MMM) - critical for coral bleaching thresholds
calculate_mmm <- function(data) {
  data %>%
    mutate(
      year = year(Date),
      month = month(Date)
    ) %>%
    group_by(Site, Variable, Depth, year, month) %>%
    summarise(monthly_mean = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Site, Variable, Depth, month) %>%
    summarise(climatological_monthly_mean = mean(monthly_mean, na.rm = TRUE), .groups = "drop") %>%
    group_by(Site, Variable, Depth) %>%
    summarise(mmm = max(climatological_monthly_mean, na.rm = TRUE), .groups = "drop")
}

# Function to calculate monthly variability (coefficient of variation of monthly means)
calculate_monthly_variability <- function(data) {
  data %>%
    mutate(
      year = year(Date),
      month = month(Date)
    ) %>%
    group_by(Site, Variable, Depth, year, month) %>%
    summarise(monthly_mean = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Site, Variable, Depth, month) %>%
    summarise(climatological_monthly_mean = mean(monthly_mean, na.rm = TRUE), .groups = "drop") %>%
    group_by(Site, Variable, Depth) %>%
    summarise(monthly_variability = cv(climatological_monthly_mean), .groups = "drop")
}

# Function to calculate mean annual values
calculate_mean_annual <- function(data) {
  data %>%
    mutate(year = year(Date)) %>%
    group_by(Site, Variable, Depth, year) %>%
    summarise(annual_mean = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Site, Variable, Depth) %>%
    summarise(mean_annual = mean(annual_mean, na.rm = TRUE), .groups = "drop")
}

# Function to calculate same metrics for non-temperature variables
# (MMM, mean annual, monthly variability - same as temperature)
calculate_same_metrics_for_other_vars <- function(data) {
  # Calculate MMM (Maximum Monthly Mean)
  mmm_result <- calculate_mmm(data)
  
  # Calculate mean annual
  annual_result <- calculate_mean_annual(data)
  
  # Calculate monthly variability
  variability_result <- calculate_monthly_variability(data)
  
  # Combine all metrics
  combined_result <- mmm_result %>%
    left_join(annual_result, by = c("Site", "Variable", "Depth")) %>%
    left_join(variability_result, by = c("Site", "Variable", "Depth"))
  
  return(combined_result)
}

# Function to clean numeric columns (handle missing data placeholders)
clean_numeric <- function(x) {
  # Replace common non-numeric placeholders with NA
  x_clean <- str_replace_all(x, c("^-$" = NA_character_, 
                                  "^--$" = NA_character_,
                                  "^NA$" = NA_character_,
                                  "^null$" = NA_character_,
                                  "^NULL$" = NA_character_,
                                  "^$" = NA_character_))  # empty strings
  as.numeric(x_clean)
}

# Main aggregation function
aggregate_environmental_data <- function(file_path) {
  
  # Read the data
  cat("Reading environmental data...\n")
  env_data <- read_csv(file_path) %>%
    mutate(
      Date = as.Date(Date),
      Depth = as.numeric(Depth),
      mean_value = clean_numeric(mean),  # Clean and rename to avoid conflicts
      median = clean_numeric(median),
      lowest = clean_numeric(lowest),
      highest = clean_numeric(highest),
      Latitude = as.numeric(Latitude),
      Longitude = as.numeric(Longitude)
    ) %>%
    # Remove rows with missing critical data
    filter(!is.na(mean_value) & !is.na(Date) & !is.na(Variable))
  
  cat("Data loaded:", nrow(env_data), "records\n")
  cat("Unique sites:", length(unique(env_data$Site)), "\n")
  cat("Unique variables:", length(unique(env_data$Variable)), "\n")
  cat("Unique depths:", paste(sort(unique(env_data$Depth)), collapse = ", "), "\n")
  cat("Date range:", min(env_data$Date), "to", max(env_data$Date), "\n")
  
  # Get site metadata
  site_metadata <- env_data %>%
    group_by(Site) %>%
    summarise(
      Latitude = first(Latitude),
      Longitude = first(Longitude),
      .groups = "drop"
    )
  
  # Identify temperature and other variables
  temp_vars <- c("temp", "temperature", "sst", "sea_surface_temperature", "water_temperature")
  temperature_variables <- unique(env_data$Variable)[tolower(unique(env_data$Variable)) %in% temp_vars]
  other_variables <- setdiff(unique(env_data$Variable), temperature_variables)
  
  cat("Temperature variables:", paste(temperature_variables, collapse = ", "), "\n")
  cat("Other variables:", paste(other_variables, collapse = ", "), "\n")
  
  # Process temperature variables
  if (length(temperature_variables) > 0) {
    cat("Processing temperature variables...\n")
    temp_data <- env_data %>% filter(Variable %in% temperature_variables)
    
    # Calculate temperature metrics
    temp_mmm <- calculate_mmm(temp_data)
    temp_annual <- calculate_mean_annual(temp_data)
    temp_variability <- calculate_monthly_variability(temp_data)
    
    # Combine temperature metrics
    temp_aggregated <- temp_mmm %>%
      left_join(temp_annual, by = c("Site", "Variable", "Depth")) %>%
      left_join(temp_variability, by = c("Site", "Variable", "Depth"))
    
  } else {
    temp_aggregated <- data.frame()
  }
  
  # Process other variables (using the same metrics as temperature)
  if (length(other_variables) > 0) {
    cat("Processing other variables with same metrics as temperature...\n")
    other_data <- env_data %>% filter(Variable %in% other_variables)
    
    # Calculate same metrics (MMM, mean annual, monthly variability) for other variables
    other_aggregated <- calculate_same_metrics_for_other_vars(other_data)
    
  } else {
    other_aggregated <- data.frame()
  }
  
  # Combine all aggregated data
  all_aggregated <- bind_rows(temp_aggregated, other_aggregated)
  
  # Add site metadata
  final_data <- all_aggregated %>%
    left_join(site_metadata, by = "Site") %>%
    # Reorder columns for clarity
    select(Site, Latitude, Longitude, Variable, Depth, everything()) %>%
    arrange(Site, Variable, Depth)
  
  cat("Aggregation complete!\n")
  cat("Final dataset:", nrow(final_data), "records\n")
  
  return(final_data)
}

# Function to create wide format for population genomics analysis
create_wide_format <- function(aggregated_data) {
  
  # Create variable names that include depth information
  wide_data <- aggregated_data %>%
    mutate(
      depth_suffix = paste0("_", abs(Depth), "m"),
      variable_depth = paste0(Variable, depth_suffix)
    ) %>%
    select(-Depth, -Variable, -depth_suffix)
  
  # Pivot to wide format for each metric
  available_metrics <- c("mmm", "mean_annual", "monthly_variability")
  
  wide_list <- list()
  
  for (metric in available_metrics) {
    if (metric %in% names(wide_data)) {
      wide_metric <- wide_data %>%
        select(Site, Latitude, Longitude, variable_depth, all_of(metric)) %>%
        pivot_wider(
          names_from = variable_depth,
          values_from = all_of(metric),
          names_prefix = paste0(metric, "_")
        )
      wide_list[[metric]] <- wide_metric
    }
  }
  
  # Combine all metrics
  if (length(wide_list) > 0) {
    wide_combined <- wide_list[[1]]
    if (length(wide_list) > 1) {
      for (i in 2:length(wide_list)) {
        wide_combined <- wide_combined %>%
          full_join(wide_list[[i]], by = c("Site", "Latitude", "Longitude"))
      }
    }
  } else {
    wide_combined <- data.frame()
  }
  
  return(wide_combined)
}

# Usage example:
# Replace "your_data.csv" with your actual file path
file_path <- "metadata/GBR1_Hydro_monthly_subset.csv"

# Run the aggregation
aggregated_data <- aggregate_environmental_data(file_path)

# View the results
head(aggregated_data)

# Save the aggregated data
write_csv(aggregated_data, "eReefs_hydro_monthly_aggregated.csv")

# Create wide format for population genomics analysis
wide_data <- create_wide_format(aggregated_data)
write_csv(wide_data, "eReefs_hydro_monthly_aggregated_wide.csv")

# Summary statistics
print(summary(aggregated_data))