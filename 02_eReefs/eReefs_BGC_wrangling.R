# Byrne, 2025
# Environmental data wrangling - Bio-geochemical
# Simplified version for pre-filtered depth data (-9m and -2.35m)

# Load required libraries
library(tidyverse)
library(lubridate)

setwd("~/Documents/GitHub/shysDelineation-WGS/metadata")

# Read the data
data <- read_csv("eReefs_BGC_monthly.csv")

# Date parsing
data <- data %>%
  mutate(
    DateTime = as.POSIXct(`Aggregated Date/Time`, format = "%Y-%m-%dT%H:%M"),
    Year = year(DateTime),
    Month = month(DateTime),
    YearMonth = format(DateTime, "%Y-%m")
  )

# Calculate monthly summaries for each site
monthly_by_site <- data %>%
  group_by(`Site Name`, Latitude, Longitude, Depth, Variable, YearMonth) %>%
  summarise(
    # Monthly values for each site
    monthly_mean = mean(mean, na.rm = TRUE),
    monthly_median = median(median, na.rm = TRUE),
    monthly_max = max(highest, na.rm = TRUE),
    monthly_min = min(lowest, na.rm = TRUE),
    monthly_p5 = mean(p5, na.rm = TRUE),
    monthly_p95 = mean(p95, na.rm = TRUE),
    .groups = "drop"
  )


# Calculate overall statistics for each site across all months
site_summary <- monthly_by_site %>%
  group_by(`Site Name`, Latitude, Longitude, Depth, Variable) %>%
  summarise(
    # Number of months with data
    n_months = n(),
    
    # Mean of monthly maximum values
    mean_monthly_max = mean(monthly_max, na.rm = TRUE),
    
    # Median of monthly maximum values
    median_monthly_max = median(monthly_max, na.rm = TRUE),
    
    # Variability measures across months
    sd_monthly_mean = sd(monthly_mean, na.rm = TRUE),
    cv_monthly_mean = (sd(monthly_mean, na.rm = TRUE) / mean(monthly_mean, na.rm = TRUE)) * 100,
    range_monthly_mean = max(monthly_mean, na.rm = TRUE) - min(monthly_mean, na.rm = TRUE),
    iqr_monthly_mean = IQR(monthly_mean, na.rm = TRUE),
    
    # Overall mean and median
    overall_mean = mean(monthly_mean, na.rm = TRUE),
    overall_median = median(monthly_median, na.rm = TRUE),
    
    # Min and max across all months
    overall_min = min(monthly_min, na.rm = TRUE),
    overall_max = max(monthly_max, na.rm = TRUE),
    
    .groups = "drop"
  )

# View the site summary
print(site_summary)


# Filter based on site name endings and depths
site_summary_filtered <- site_summary %>%
  filter(
    # Keep -3 depth for sites ending in S
    (str_ends(`Site Name`, "S") & Depth == -3) |
      # Keep -8.8 depth for sites ending in D
      (str_ends(`Site Name`, "D") & Depth == -8.8) |
      # Keep any depth that is 99999.9 (special case)
      (Depth == 99999.9)
  )

# Convert to wide format
# Create column names that combine variable and statistic
site_summary_wide <- site_summary_filtered %>%
  select(`Site Name`, Latitude, Longitude, Variable, 
         mean_monthly_max, median_monthly_max, cv_monthly_mean, 
         overall_mean, overall_median, sd_monthly_mean) %>%
  pivot_wider(
    names_from = Variable,
    values_from = c(mean_monthly_max, median_monthly_max, cv_monthly_mean, 
                    overall_mean, overall_median, sd_monthly_mean),
    names_glue = "{Variable}_{.value}"
  )

# Reorder columns to group by variable
# Get variable names
vars <- unique(site_summary_filtered$Variable)

# Create ordered column names
ordered_cols <- c("Site Name", "Latitude", "Longitude")
for(var in vars) {
  ordered_cols <- c(ordered_cols,
                    paste0(var, "_mean_monthly_max"),
                    paste0(var, "_median_monthly_max"),
                    paste0(var, "_overall_mean"),
                    paste0(var, "_overall_median"),
                    paste0(var, "_cv_monthly_mean"),
                    paste0(var, "_sd_monthly_mean"))
}

# Reorder columns if they exist
existing_cols <- intersect(ordered_cols, names(site_summary_wide))
site_summary_wide <- site_summary_wide %>%
  select(all_of(existing_cols))

# Save the wide format summary
write_csv(site_summary_wide, "eReefs_BGC_monthly_aggregated_wide.csv")


