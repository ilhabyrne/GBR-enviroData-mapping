# =============================================================================
# Byrne, 2025
# Environmental data concatenation
# =============================================================================

library(dplyr)
library(readr)

# Load data
setwd("~/Documents/GitHub/GBR-enviroData-mapping")
env_dat <- read.csv("~/Documents/GitHub/GBR-enviroData-mapping/02_eReefs/synthesis/RRAP_enviro_metadata_targetSpp_2025-12-09.csv")

# Summarize by species and reef, calculating means for all numeric columns
env_by_reef <- env_dat %>%
  group_by(scientificName, locality) %>%
  summarise(
    n_individuals = n(),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Split by species and write separate CSVs
species_list <- unique(env_by_reef$scientificName)

for (sp in species_list) {
  env_by_reef %>%
    filter(scientificName == sp) %>%
    write_csv(paste0("env_reef_means_", sp, ".csv"))
}