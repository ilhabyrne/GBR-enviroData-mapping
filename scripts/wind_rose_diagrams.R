# Byrne, 2025
# Create rose diagrams using GBR monthly wind speeds

# Wind direction calculated based on trigonometry:
# wd = (270 - atan2(wspeed_v, wspeed_u) * 180/pi) %% 360


source("wind_rose_functions.R")

# Run functions & compare sites
results <- run_wind_analysis("eReefs_hydro_monthly_wind.csv", 
                             site_ids = c("CBHE_BA1S", "ONLI_FL1S", "TSMA_FR1S"),
                             statistic = "mean")



# Wind speed magnitude: sqrt(u² + v²)
# Wind direction: atan2(v, u)
# Wave energy proxy: wind_speed³ (simplified wave energy relationship)
# Directional exposure: Separate east/west/north/south components

# Load the metrics functions
source("wind_exposure_functions.R")  # Save the artifact as this file

# Analyze all wind and wave metrics
results <- analyze_wind_wave_metrics("eReefs_hydro_monthly_wind.csv", 
                                     statistic = "mean")

# Access the results
wind_stats <- results$wind_metrics
wave_stats <- results$wave_metrics

# View summary
summary_table <- create_summary_table(wind_stats, wave_stats)
print(summary_table)
