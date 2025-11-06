# Enhanced Vector Statistics with Biological Metrics
# For coral population genomics and connectivity analyses
# =============================================================================

#' Calculate comprehensive biologically-relevant metrics for winds and currents
#' 
#' @param data Monthly data with u/v components
#' @param statistic Which column to use ("mean", "median", "p5", "p95")
#' @return Vector statistics with biological metrics
calculate_vector_statistics_biological <- function(data, statistic = "mean") {
  
  data <- add_temporal_variables(data)
  
  data <- data %>%
    rename(value = all_of(statistic))
  
  # =============================================================================
  # WIND METRICS
  # =============================================================================
  
  wind_data <- data %>%
    filter(Variable %in% c("wspeed_u", "wspeed_v")) %>%
    select(Site, Latitude, Longitude, Date, year, month, season, Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    filter(!is.na(wspeed_u) & !is.na(wspeed_v))
  
  if (nrow(wind_data) > 0) {
    wind_vectors <- wind_data %>%
      group_by(Site, Latitude, Longitude) %>%
      summarise(
        # =====================================================================
        # BASIC METRICS
        # =====================================================================
        
        # Mean components (vector means)
        mean_u = mean(wspeed_u, na.rm = TRUE),
        mean_v = mean(wspeed_v, na.rm = TRUE),
        
        # Mean speed (scalar mean - average of magnitudes)
        mean_wind_speed = mean(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # Temporal variability
        sd_wind_speed = sd(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        cv_wind_speed = sd(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE) / 
                       mean(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # Speed range
        min_wind_speed = min(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        max_wind_speed = max(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        wind_speed_range = max_wind_speed - min_wind_speed,
        
        # =====================================================================
        # DIRECTIONAL METRICS
        # =====================================================================
        
        # Vector magnitude (resultant vector)
        vector_magnitude = sqrt(mean_u^2 + mean_v^2),
        
        # Mean wind direction (0° = North, clockwise)
        mean_wind_direction = (270 - atan2(mean_v, mean_u) * 180/pi) %% 360,
        
        # Directional constancy (0 = highly variable, 1 = constant direction)
        # This is R/V̄ where R = resultant length, V̄ = mean speed
        directional_constancy = sqrt(mean_u^2 + mean_v^2) / 
                               mean(sqrt(wspeed_u^2 + wspeed_v^2), na.rm = TRUE),
        
        # Directional variability (inverse of constancy)
        directional_variability = 1 - directional_constancy,
        
        # Predominant wind sector (8 directions)
        predominant_wind_sector = case_when(
          mean_wind_direction >= 337.5 | mean_wind_direction < 22.5 ~ "N",
          mean_wind_direction >= 22.5 & mean_wind_direction < 67.5 ~ "NE",
          mean_wind_direction >= 67.5 & mean_wind_direction < 112.5 ~ "E",
          mean_wind_direction >= 112.5 & mean_wind_direction < 157.5 ~ "SE",
          mean_wind_direction >= 157.5 & mean_wind_direction < 202.5 ~ "S",
          mean_wind_direction >= 202.5 & mean_wind_direction < 247.5 ~ "SW",
          mean_wind_direction >= 247.5 & mean_wind_direction < 292.5 ~ "W",
          TRUE ~ "NW"
        ),
        
        # =====================================================================
        # GBR-SPECIFIC METRICS
        # =====================================================================
        
        # Trade wind component (SE trade winds = 135° for GBR)
        # Positive = aligned with trades, negative = against trades
        trade_wind_component = mean_u * cos(135 * pi/180) + mean_v * sin(135 * pi/180),
        
        # Alongshore component (N-S, assuming reef orientation)
        # Positive = northward, negative = southward
        alongshore_wind_component = mean_v,
        
        # Cross-shore component (E-W, offshore/onshore)
        # Positive = eastward (offshore), negative = westward (onshore)
        crossshore_wind_component = mean_u,
        
        # Wind exposure to dominant quadrants
        easterly_exposure = mean(pmax(0, wspeed_u), na.rm = TRUE),  # Winds from east
        westerly_exposure = mean(pmax(0, -wspeed_u), na.rm = TRUE), # Winds from west
        northerly_exposure = mean(pmax(0, wspeed_v), na.rm = TRUE), # Winds from north
        southerly_exposure = mean(pmax(0, -wspeed_v), na.rm = TRUE), # Winds from south
        
        # =====================================================================
        # SEASONAL METRICS (Biologically Important)
        # =====================================================================
        
        # Seasonal means
        summer_wind_speed = mean(sqrt(wspeed_u[season == "Summer"]^2 + 
                                     wspeed_v[season == "Summer"]^2), na.rm = TRUE),
        autumn_wind_speed = mean(sqrt(wspeed_u[season == "Autumn"]^2 + 
                                     wspeed_v[season == "Autumn"]^2), na.rm = TRUE),
        winter_wind_speed = mean(sqrt(wspeed_u[season == "Winter"]^2 + 
                                     wspeed_v[season == "Winter"]^2), na.rm = TRUE),
        spring_wind_speed = mean(sqrt(wspeed_u[season == "Spring"]^2 + 
                                     wspeed_v[season == "Spring"]^2), na.rm = TRUE),
        
        # Seasonal range (environmental variability)
        seasonal_wind_range = pmax(summer_wind_speed, autumn_wind_speed,
                                  winter_wind_speed, spring_wind_speed, na.rm = TRUE) -
                             pmin(summer_wind_speed, autumn_wind_speed,
                                  winter_wind_speed, spring_wind_speed, na.rm = TRUE),
        
        # Dominant season
        dominant_wind_season = case_when(
          summer_wind_speed == pmax(summer_wind_speed, autumn_wind_speed,
                                   winter_wind_speed, spring_wind_speed, na.rm = TRUE) ~ "Summer",
          winter_wind_speed == pmax(summer_wind_speed, autumn_wind_speed,
                                   winter_wind_speed, spring_wind_speed, na.rm = TRUE) ~ "Winter",
          autumn_wind_speed == pmax(summer_wind_speed, autumn_wind_speed,
                                   winter_wind_speed, spring_wind_speed, na.rm = TRUE) ~ "Autumn",
          TRUE ~ "Spring"
        ),
        
        # Summer-winter difference (key for breeding season vs non-breeding)
        summer_winter_wind_diff = abs(summer_wind_speed - winter_wind_speed),
        
        # =====================================================================
        # BIOLOGICAL PROXIES
        # =====================================================================
        
        # Wave energy proxy (wind-generated waves)
        # Wave energy scales with wind speed cubed
        mean_wave_energy_proxy = mean((sqrt(wspeed_u^2 + wspeed_v^2))^3, na.rm = TRUE),
        
        # Extreme event frequency (>95th percentile)
        extreme_wind_frequency = sum(sqrt(wspeed_u^2 + wspeed_v^2) > 
                                    quantile(sqrt(wspeed_u^2 + wspeed_v^2), 0.95, na.rm = TRUE),
                                    na.rm = TRUE) / n() * 100,
        
        # Calm period frequency (<2 m/s)
        calm_wind_frequency = sum(sqrt(wspeed_u^2 + wspeed_v^2) < 2, na.rm = TRUE) / n() * 100,
        
        # =====================================================================
        # LARVAL TRANSPORT METRICS
        # =====================================================================
        
        # Transport potential (for larval connectivity)
        # High constancy × high speed = strong directional transport
        wind_transport_potential = directional_constancy * mean_wind_speed,
        
        # Dispersal symmetry (is transport balanced in all directions?)
        # Low = asymmetric (directional), high = symmetric (omnidirectional)
        wind_dispersal_symmetry = 1 - abs(mean_u / mean_wind_speed),
        
        # Retention potential (low speed × low constancy = local retention)
        wind_retention_potential = (1 - directional_constancy) * (1 / mean_wind_speed),
        
        # =====================================================================
        # SAMPLE SIZE
        # =====================================================================
        
        n_months = n(),
        n_years = n_distinct(year),
        
        .groups = "drop"
      ) %>%
      mutate(
        # ===================================================================
        # CATEGORICAL VARIABLES (for statistical models)
        # ===================================================================
        
        # Wind exposure category
        wind_exposure_category = cut(
          mean_wind_speed,
          breaks = quantile(mean_wind_speed, c(0, 0.33, 0.67, 1), na.rm = TRUE),
          labels = c("Sheltered", "Moderate", "Exposed"),
          include.lowest = TRUE
        ),
        
        # Directional consistency category
        directional_consistency_category = cut(
          directional_constancy,
          breaks = c(0, 0.4, 0.7, 1),
          labels = c("Variable", "Moderate", "Consistent"),
          include.lowest = TRUE
        ),
        
        # Seasonal variability category
        seasonal_variability_category = cut(
          seasonal_wind_range,
          breaks = c(0, 1, 2, Inf),
          labels = c("Low", "Moderate", "High"),
          include.lowest = TRUE
        ),
        
        # Wind regime type (based on speed and constancy)
        wind_regime_type = case_when(
          mean_wind_speed > 5 & directional_constancy > 0.6 ~ "Steady_Strong",
          mean_wind_speed > 5 & directional_constancy <= 0.6 ~ "Variable_Strong",
          mean_wind_speed <= 5 & directional_constancy > 0.6 ~ "Steady_Weak",
          TRUE ~ "Variable_Weak"
        ),
        
        variable_type = "wind"
      )
  } else {
    wind_vectors <- NULL
  }
  
  # =============================================================================
  # CURRENT METRICS (Same structure as wind)
  # =============================================================================
  
  current_data <- data %>%
    filter(Variable %in% c("u", "v")) %>%
    select(Site, Latitude, Longitude, Date, year, month, season, Variable, value) %>%
    pivot_wider(names_from = Variable, values_from = value) %>%
    filter(!is.na(u) & !is.na(v))
  
  if (nrow(current_data) > 0) {
    current_vectors <- current_data %>%
      group_by(Site, Latitude, Longitude) %>%
      summarise(
        # Basic metrics
        mean_u = mean(u, na.rm = TRUE),
        mean_v = mean(v, na.rm = TRUE),
        mean_current_speed = mean(sqrt(u^2 + v^2), na.rm = TRUE),
        sd_current_speed = sd(sqrt(u^2 + v^2), na.rm = TRUE),
        cv_current_speed = sd(sqrt(u^2 + v^2), na.rm = TRUE) / 
                          mean(sqrt(u^2 + v^2), na.rm = TRUE),
        min_current_speed = min(sqrt(u^2 + v^2), na.rm = TRUE),
        max_current_speed = max(sqrt(u^2 + v^2), na.rm = TRUE),
        current_speed_range = max_current_speed - min_current_speed,
        
        # Directional metrics
        vector_magnitude = sqrt(mean_u^2 + mean_v^2),
        mean_current_direction = (270 - atan2(mean_v, mean_u) * 180/pi) %% 360,
        directional_constancy = sqrt(mean_u^2 + mean_v^2) / 
                               mean(sqrt(u^2 + v^2), na.rm = TRUE),
        directional_variability = 1 - directional_constancy,
        predominant_current_sector = case_when(
          mean_current_direction >= 337.5 | mean_current_direction < 22.5 ~ "N",
          mean_current_direction >= 22.5 & mean_current_direction < 67.5 ~ "NE",
          mean_current_direction >= 67.5 & mean_current_direction < 112.5 ~ "E",
          mean_current_direction >= 112.5 & mean_current_direction < 157.5 ~ "SE",
          mean_current_direction >= 157.5 & mean_current_direction < 202.5 ~ "S",
          mean_current_direction >= 202.5 & mean_current_direction < 247.5 ~ "SW",
          mean_current_direction >= 247.5 & mean_current_direction < 292.5 ~ "W",
          TRUE ~ "NW"
        ),
        
        # GBR-specific metrics
        trade_wind_component = mean_u * cos(135 * pi/180) + mean_v * sin(135 * pi/180),
        alongshore_current_component = mean_v,
        crossshore_current_component = mean_u,
        easterly_exposure = mean(pmax(0, u), na.rm = TRUE),
        westerly_exposure = mean(pmax(0, -u), na.rm = TRUE),
        northerly_exposure = mean(pmax(0, v), na.rm = TRUE),
        southerly_exposure = mean(pmax(0, -v), na.rm = TRUE),
        
        # Seasonal metrics
        summer_current_speed = mean(sqrt(u[season == "Summer"]^2 + 
                                       v[season == "Summer"]^2), na.rm = TRUE),
        autumn_current_speed = mean(sqrt(u[season == "Autumn"]^2 + 
                                       v[season == "Autumn"]^2), na.rm = TRUE),
        winter_current_speed = mean(sqrt(u[season == "Winter"]^2 + 
                                       v[season == "Winter"]^2), na.rm = TRUE),
        spring_current_speed = mean(sqrt(u[season == "Spring"]^2 + 
                                       v[season == "Spring"]^2), na.rm = TRUE),
        seasonal_current_range = pmax(summer_current_speed, autumn_current_speed,
                                     winter_current_speed, spring_current_speed, na.rm = TRUE) -
                                pmin(summer_current_speed, autumn_current_speed,
                                     winter_current_speed, spring_current_speed, na.rm = TRUE),
        dominant_current_season = case_when(
          summer_current_speed == pmax(summer_current_speed, autumn_current_speed,
                                      winter_current_speed, spring_current_speed, na.rm = TRUE) ~ "Summer",
          winter_current_speed == pmax(summer_current_speed, autumn_current_speed,
                                      winter_current_speed, spring_current_speed, na.rm = TRUE) ~ "Winter",
          autumn_current_speed == pmax(summer_current_speed, autumn_current_speed,
                                      winter_current_speed, spring_current_speed, na.rm = TRUE) ~ "Autumn",
          TRUE ~ "Spring"
        ),
        summer_winter_current_diff = abs(summer_current_speed - winter_current_speed),
        
        # Biological proxies
        extreme_current_frequency = sum(sqrt(u^2 + v^2) > 
                                       quantile(sqrt(u^2 + v^2), 0.95, na.rm = TRUE),
                                       na.rm = TRUE) / n() * 100,
        calm_current_frequency = sum(sqrt(u^2 + v^2) < 0.05, na.rm = TRUE) / n() * 100,
        
        # Larval transport metrics
        current_transport_potential = directional_constancy * mean_current_speed,
        current_dispersal_symmetry = 1 - abs(mean_u / mean_current_speed),
        current_retention_potential = (1 - directional_constancy) * (1 / mean_current_speed),
        
        # Sample size
        n_months = n(),
        n_years = n_distinct(year),
        
        .groups = "drop"
      ) %>%
      mutate(
        # Categorical variables
        current_exposure_category = cut(
          mean_current_speed,
          breaks = quantile(mean_current_speed, c(0, 0.33, 0.67, 1), na.rm = TRUE),
          labels = c("Low_Flow", "Moderate_Flow", "High_Flow"),
          include.lowest = TRUE
        ),
        
        directional_consistency_category = cut(
          directional_constancy,
          breaks = c(0, 0.4, 0.7, 1),
          labels = c("Variable", "Moderate", "Consistent"),
          include.lowest = TRUE
        ),
        
        seasonal_variability_category = cut(
          seasonal_current_range,
          breaks = c(0, 0.05, 0.1, Inf),
          labels = c("Low", "Moderate", "High"),
          include.lowest = TRUE
        ),
        
        current_regime_type = case_when(
          mean_current_speed > 0.1 & directional_constancy > 0.6 ~ "Steady_Strong",
          mean_current_speed > 0.1 & directional_constancy <= 0.6 ~ "Variable_Strong",
          mean_current_speed <= 0.1 & directional_constancy > 0.6 ~ "Steady_Weak",
          TRUE ~ "Variable_Weak"
        ),
        
        variable_type = "current"
      )
  } else {
    current_vectors <- NULL
  }
  
  # Combine wind and current metrics
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

cat("\nEnhanced biological metrics function loaded\n")
cat("This will replace the standard calculate_vector_statistics() function\n\n")
