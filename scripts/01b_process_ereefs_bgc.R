# Byrne, 2025
# Environmental data wrangling - Bio-geochemical
# Adapted for depth-specific CSV files

library(tidyverse)
library(lubridate)

setwd("~/Documents/GitHub/shysDelineation-WGS/metadata")

# Read and combine all depth files
# Adjust filenames to match your actual file names
depth_files <- list.files(pattern = "*.csv", full.names = TRUE)

# Or specify them explicitly:
# depth_files <- c("depth1.csv", "depth2.csv", "depth3.csv", "depth4.csv", "depth5.csv")

data <- depth_files %>%
  map_dfr(read_csv)

# Standardize column names (adjust if your columns differ)
data <- data %>%
  rename(
    DateTime = data,
    Variable = variable,
    Depth = depth,
    `Site Name` = `site name`,
    Latitude = lat,
    Longitude = lon
  )

# Date parsing - adjust format string to match your date format
data <- data %>%
  mutate(
    DateTime = as.POSIXct(DateTime, format = "%Y-%m-%dT%H:%M"),  # adjust format as needed
    Year = year(DateTime),
    Month = month(DateTime),
    YearMonth = format(DateTime, "%Y-%m")
  )

# Calculate monthly summaries for each site
monthly_by_site <- data %>%
  group_by(`Site Name`, Latitude, Longitude, Depth, Variable, YearMonth) %>%
  summarise(
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
    n_months = n(),
    mean_monthly_max = mean(monthly_max, na.rm = TRUE),
    median_monthly_max = median(monthly_max, na.rm = TRUE),
    sd_monthly_mean = sd(monthly_mean, na.rm = TRUE),
    cv_monthly_mean = (sd(monthly_mean, na.rm = TRUE) / mean(monthly_mean, na.rm = TRUE)) * 100,
    range_monthly_mean = max(monthly_mean, na.rm = TRUE) - min(monthly_mean, na.rm = TRUE),
    iqr_monthly_mean = IQR(monthly_mean, na.rm = TRUE),
    overall_mean = mean(monthly_mean, na.rm = TRUE),
    overall_median = median(monthly_median, na.rm = TRUE),
    overall_min = min(monthly_min, na.rm = TRUE),
    overall_max = max(monthly_max, na.rm = TRUE),
    .groups = "drop"
  )

# Since sites are unique per depth, no depth filtering needed
# Convert directly to wide format
site_summary_wide <- site_summary %>%
  select(`Site Name`, Latitude, Longitude, Depth, Variable,
         mean_monthly_max, median_monthly_max, cv_monthly_mean,
         overall_mean, overall_median, sd_monthly_mean) %>%
  pivot_wider(
    names_from = Variable,
    values_from = c(mean_monthly_max, median_monthly_max, cv_monthly_mean,
                    overall_mean, overall_median, sd_monthly_mean),
    names_glue = "{Variable}_{.value}"
  )

# Reorder columns to group by variable
vars <- unique(site_summary$Variable)

ordered_cols <- c("Site Name", "Latitude", "Longitude", "Depth")
for(var in vars) {
  ordered_cols <- c(ordered_cols,
                    paste0(var, "_mean_monthly_max"),
                    paste0(var, "_median_monthly_max"),
                    paste0(var, "_overall_mean"),
                    paste0(var, "_overall_median"),
                    paste0(var, "_cv_monthly_mean"),
                    paste0(var, "_sd_monthly_mean"))
}

existing_cols <- intersect(ordered_cols, names(site_summary_wide))
site_summary_wide <- site_summary_wide %>%
  select(all_of(existing_cols))

write_csv(site_summary_wide, "eReefs_BGC_monthly_aggregated_wide.csv")