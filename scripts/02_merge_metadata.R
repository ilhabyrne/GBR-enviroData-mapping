# Byrne, 2025
# Merge metadata, temp summary stats, and wind/current summary stats

# Load libraries
library(dplyr)

# Set workdir
setwd("~/Documents/GitHub/GBR-enviroData-mapping/02_eReefs/processed/")

# Read the files
small_df <- read.csv("GBR1_dailytemp_2025-09-15_temperature_summary_stats.csv") # summary stats
large_df <- read.csv("RRAP_metadata_coordinates_all_with_ereefs_indices_2025-09-15.csv") # all samples

# Left join - keeps all rows from large_df, matches info from small_df
merged <- left_join(large_df, small_df, 
                    by = c("long_GBR1", "lat_GBR1", "k_value_GBR1"))


# Read the files
wind_df <- read.csv("currents_wind_vector_stats_subset.csv") # summary stats
meta_df <- read.csv("../../01_coords/RRAP_metadata_coordinates_all_2025-09-05_clean.csv") # all samples

# Left join - keeps all rows from large_df, matches info from small_df
merged_wind <- left_join(meta_df, wind_df, 
                    by = c("decimalLatitude", "decimalLongitude"))


# Left join - keeps all rows from large_df, matches info from small_df
merged_all <- left_join(merged_wind, merged, 
                         by = c("individualID"))

write.csv(merged_all, "RRAP_enviro_metadata_all_2025-12-09.csv")
