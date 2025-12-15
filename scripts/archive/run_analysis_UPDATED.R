# IMPLEMENTATION SCRIPT - Handles All Sites, Then Filters to Coral Sites
# For Pre-Aggregated Monthly Data
# =============================================================================

library(tidyverse)
library(lubridate)

# Source functions
source("process_ereefs_FINAL.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Your eReefs files (edit these paths)
ereefs_files <- c(
  "file1" = "data/RRAP_eReefs_GBR1_0.5m_currents_wind.csv",  # EDIT
  "file2" = "data/RRAP_eReefs_GBR1_0.5m_currents_wind_2.csv",
  "file3" = "data/RRAP_eReefs_GBR1_0.5m_currents_wind_3.csv",
  "file4" = "data/RRAP_eReefs_GBR1_0.5m_currents_wind_4.csv",
  "file5" = "data/RRAP_eReefs_GBR1_0.5m_currents_wind_5.csv"
)

# Your coral sampling sites (the actual sites where you sampled S. hystrix)
coral_sites <- c("CBHE", "CBLM", "PAPE", "OCLB", "OCCH", 
                 "OCDA", "ONMO", "ONLI", "TSAU", "TSMA")

# Output directory
output_dir <- "results/ereefs_final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

# =============================================================================
# STEP 1: PROCESS ALL SITES (Don't filter yet)
# =============================================================================

cat("\n### PROCESSING ALL SITES WITH MEAN STATISTIC ###\n\n")

# Process ALL sites in the files first
results_mean <- process_ereefs_monthly_data(
  file_paths = ereefs_files,
  statistic = "mean",
  coral_sites = NULL  # NULL = process all sites
)

# Export full results (all sites)
export_results(
  results_mean,
  output_dir = file.path(output_dir, "all_sites_mean"),
  prefix = "all_sites"
)

cat("\n### PROCESSING ALL SITES WITH P95 STATISTIC ###\n\n")

# Process extreme conditions for all sites
results_p95 <- process_ereefs_monthly_data(
  file_paths = ereefs_files,
  statistic = "p95",
  coral_sites = NULL
)

export_results(
  results_p95,
  output_dir = file.path(output_dir, "all_sites_p95"),
  prefix = "all_sites_extreme"
)

# =============================================================================
# STEP 2: INSPECT ALL AVAILABLE SITES
# =============================================================================

cat("\n### AVAILABLE SITES IN DATA ###\n\n")

all_sites <- results_mean$vector_statistics %>%
  filter(variable_type == "wind") %>%
  select(Site, Latitude, Longitude, mean_speed, directional_constancy) %>%
  arrange(Latitude)

cat("Total sites found:", nrow(all_sites), "\n\n")

# Save list of all sites
write_csv(all_sites, file.path(output_dir, "all_available_sites.csv"))

cat("All site names:\n")
print(all_sites$Site)
cat("\n")

# =============================================================================
# STEP 3: CREATE SITE MATCHING HELPER
# =============================================================================

cat("### MATCHING CORAL SITES TO eREEFS SITES ###\n\n")

# Function to help match sites
find_site_matches <- function(coral_sites, all_ereefs_sites, method = "fuzzy") {
  
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
    
    # Check for exact match
    if (cs %in% all_ereefs_sites) {
      matches$exact_match[i] <- cs
      matches$matched_site[i] <- cs
      matches$match_method[i] <- "exact"
    } else {
      # Try fuzzy matching
      library(stringdist)
      distances <- stringdist(cs, all_ereefs_sites, method = "lv")
      best_idx <- which.min(distances)
      best_dist <- distances[best_idx]
      
      if (best_dist <= 5) {  # Allow up to 5 character differences
        matches$fuzzy_match[i] <- all_ereefs_sites[best_idx]
        matches$fuzzy_distance[i] <- best_dist
        matches$matched_site[i] <- all_ereefs_sites[best_idx]
        matches$match_method[i] <- "fuzzy"
      } else {
        # Check if coral site is substring of any eReefs site
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

# Try to match your coral sites
site_matching <- find_site_matches(coral_sites, all_sites$Site)

# Save matching table
write_csv(site_matching, file.path(output_dir, "site_matching_table.csv"))

cat("Site matching results:\n")
print(site_matching)
cat("\n")

# Count matches
n_matched <- sum(site_matching$match_method != "NO_MATCH")
n_unmatched <- sum(site_matching$match_method == "NO_MATCH")

cat("Summary:\n")
cat("  Matched:", n_matched, "/", length(coral_sites), "\n")
cat("  Unmatched:", n_unmatched, "\n")

if (n_unmatched > 0) {
  unmatched_sites <- site_matching$coral_site[site_matching$match_method == "NO_MATCH"]
  cat("\n⚠️  Unmatched sites:", paste(unmatched_sites, collapse = ", "), "\n")
  cat("  → Check 'all_available_sites.csv' for similar names\n")
  cat("  → Manually create matching table if needed\n")
}
cat("\n")

# =============================================================================
# STEP 4: CREATE ENVIRONMENTAL MATRIX FOR MATCHED CORAL SITES
# =============================================================================

cat("### CREATING ENVIRONMENTAL MATRIX FOR CORAL SITES ###\n\n")

# Get successfully matched sites
matched_coral_sites <- site_matching %>%
  filter(match_method != "NO_MATCH") %>%
  select(coral_site, ereefs_site = matched_site)

if (nrow(matched_coral_sites) > 0) {
  
  # Extract environmental data for matched sites
  env_matrix <- results_mean$vector_statistics %>%
    filter(variable_type == "wind") %>%
    filter(Site %in% matched_coral_sites$ereefs_site) %>%
    # Add seasonal data
    left_join(
      results_mean$seasonal_patterns %>%
        filter(Variable == "wspeed_u") %>%
        select(Site, seasonal_range, dominant_season),
      by = "Site"
    ) %>%
    # Add extreme values
    left_join(
      results_p95$vector_statistics %>%
        filter(variable_type == "wind") %>%
        select(Site, extreme_speed = mean_speed, extreme_magnitude = vector_magnitude),
      by = "Site"
    ) %>%
    # Match back to coral site names
    left_join(matched_coral_sites, by = c("Site" = "ereefs_site")) %>%
    select(-Site) %>%
    relocate(coral_site, .before = 1) %>%
    # Calculate composite indices
    mutate(
      # Basic indices
      wind_exposure_index = mean_speed,
      temporal_variability = sd_speed,
      transport_potential = directional_constancy * mean_speed,
      
      # Advanced indices
      hydro_index = scale(mean_speed)[,1],
      disturbance_index = scale(temporal_variability)[,1] + 
                         scale(coalesce(extreme_speed - mean_speed, 0))[,1],
      predictability = 1 / (1 + temporal_variability / mean_speed),
      
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
    )
  
  # Export environmental matrix for genomics
  write_csv(env_matrix, 
            file.path(output_dir, "environmental_matrix_for_RDA.csv"))
  
  cat("✓ Environmental matrix created for", nrow(env_matrix), "sites\n")
  cat("  File:", file.path(output_dir, "environmental_matrix_for_RDA.csv"), "\n\n")
  
} else {
  cat("⚠️  No sites matched - check site names manually\n")
  cat("  See 'site_matching_table.csv' and 'all_available_sites.csv'\n\n")
  env_matrix <- NULL
}

# =============================================================================
# STEP 5: VISUALIZATIONS (if we have matched sites)
# =============================================================================

if (!is.null(env_matrix) && nrow(env_matrix) > 0) {
  
  cat("### CREATING VISUALIZATIONS ###\n\n")
  
  fig_dir <- file.path(output_dir, "figures")
  
  # 1. Wind speed along latitudinal gradient
  p1 <- ggplot(env_matrix,
               aes(x = abs(Latitude), y = mean_speed, 
                   color = exposure_category, label = coral_site)) +
    geom_point(size = 4) +
    geom_text(hjust = -0.2, size = 3, show.legend = FALSE) +
    geom_smooth(method = "loess", se = TRUE, color = "darkblue", 
                alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = c("Sheltered" = "green3",
                                  "Moderate" = "orange", 
                                  "Exposed" = "red3")) +
    labs(
      title = "Wind Speed Exposure Along GBR Latitudinal Gradient",
      subtitle = "S. hystrix sampling sites",
      x = "Latitude (°S)",
      y = "Mean Wind Speed (m/s)",
      color = "Exposure"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(file.path(fig_dir, "wind_latitudinal_gradient.png"), p1,
         width = 10, height = 6, dpi = 300)
  
  # 2. Directional patterns
  p2 <- ggplot(env_matrix,
               aes(x = abs(Latitude), y = directional_constancy,
                   color = predominant_sector, size = mean_speed)) +
    geom_point(alpha = 0.7) +
    geom_text(aes(label = coral_site), hjust = -0.2, vjust = 0,
              size = 3, show.legend = FALSE) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = "Wind Directional Patterns Along GBR",
      subtitle = "Size = wind speed | Color = predominant direction",
      x = "Latitude (°S)",
      y = "Directional Constancy (0-1)",
      color = "Wind\nDirection",
      size = "Wind Speed"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, "directional_patterns.png"), p2,
         width = 11, height = 6, dpi = 300)
  
  # 3. Seasonal variability
  p3 <- ggplot(env_matrix,
               aes(x = reorder(coral_site, -seasonal_range),
                   y = seasonal_range, fill = dominant_season)) +
    geom_col() +
    labs(
      title = "Seasonal Wind Variability at Sampling Sites",
      x = "Site",
      y = "Seasonal Range (m/s)",
      fill = "Strongest\nSeason"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(file.path(fig_dir, "seasonal_variability.png"), p3,
         width = 10, height = 6, dpi = 300)
  
  # 4. Mean vs extreme conditions
  if ("extreme_speed" %in% names(env_matrix)) {
    p4 <- ggplot(env_matrix,
                 aes(x = mean_speed, y = extreme_speed, label = coral_site)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(size = 4, alpha = 0.6, color = "darkred") +
      geom_text(hjust = -0.2, size = 3) +
      labs(
        title = "Mean vs Extreme Wind Conditions",
        subtitle = "Dashed line = 1:1 relationship",
        x = "Mean Wind Speed (m/s)",
        y = "Extreme Wind Speed - P95 (m/s)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, "mean_vs_extreme.png"), p4,
           width = 8, height = 8, dpi = 300)
  }
  
  # 5. Environmental space PCA
  env_for_pca <- env_matrix %>%
    select(mean_speed, temporal_variability, directional_constancy, seasonal_range) %>%
    na.omit()
  
  if (nrow(env_for_pca) >= 3) {
    pca_result <- prcomp(env_for_pca, scale = TRUE, center = TRUE)
    
    pca_data <- as.data.frame(pca_result$x) %>%
      mutate(coral_site = env_matrix$coral_site[complete.cases(env_for_pca)])
    
    var_explained <- summary(pca_result)$importance[2, 1:2] * 100
    
    p5 <- ggplot(pca_data, aes(x = PC1, y = PC2, label = coral_site)) +
      geom_point(size = 4, color = "darkblue", alpha = 0.6) +
      geom_text(hjust = -0.2, size = 3) +
      labs(
        title = "Environmental Space - Principal Components",
        subtitle = "Position of sites in multivariate environmental space",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(file.path(fig_dir, "environmental_pca.png"), p5,
           width = 9, height = 7, dpi = 300)
  }
  
  cat("✓ Visualizations created in:", fig_dir, "\n\n")
}

# =============================================================================
# STEP 6: SUMMARY REPORT
# =============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PROCESSING SUMMARY\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Data Structure:\n")
cat("  Files processed:", results_mean$processing_info$n_files, "\n")
cat("  Total sites in data:", results_mean$processing_info$n_sites, "\n")
cat("  Months in dataset:", results_mean$processing_info$n_months, "\n")
cat("  Date range:", format(results_mean$processing_info$date_range[1]), "to",
    format(results_mean$processing_info$date_range[2]), "\n\n")

cat("Coral Site Matching:\n")
cat("  Coral sites requested:", length(coral_sites), "\n")
cat("  Sites matched:", n_matched, "\n")
if (n_unmatched > 0) {
  cat("  ⚠️  Sites NOT matched:", n_unmatched, "\n")
  cat("     →", paste(site_matching$coral_site[site_matching$match_method == "NO_MATCH"], 
                       collapse = ", "), "\n")
}
cat("\n")

if (!is.null(env_matrix)) {
  cat("Environmental Gradients (matched sites):\n")
  cat("  Wind speed range:",
      round(min(env_matrix$mean_speed, na.rm = TRUE), 2), "-",
      round(max(env_matrix$mean_speed, na.rm = TRUE), 2), "m/s\n")
  cat("  Directional constancy range:",
      round(min(env_matrix$directional_constancy, na.rm = TRUE), 2), "-",
      round(max(env_matrix$directional_constancy, na.rm = TRUE), 2), "\n")
  cat("  Seasonal variability range:",
      round(min(env_matrix$seasonal_range, na.rm = TRUE), 2), "-",
      round(max(env_matrix$seasonal_range, na.rm = TRUE), 2), "m/s\n\n")
}

cat(rep("=", 70), "\n", sep = "")
cat("OUTPUT FILES\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("📁", output_dir, "/\n")
cat("├── all_sites_mean/              # All sites in data (mean)\n")
cat("├── all_sites_p95/               # All sites in data (p95)\n")
cat("├── all_available_sites.csv      # List of ALL sites ⭐\n")
cat("├── site_matching_table.csv      # How coral sites matched ⭐\n")
if (!is.null(env_matrix)) {
  cat("├── environmental_matrix_for_RDA.csv  # Matched coral sites ⭐⭐⭐\n")
}
cat("└── figures/                     # Visualizations\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("NEXT STEPS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("1. Check site matching:\n")
cat("   → Open 'site_matching_table.csv'\n")
cat("   → Verify matches are correct\n")
cat("   → If needed, manually edit matching\n\n")

cat("2. If sites didn't match automatically:\n")
cat("   → Check 'all_available_sites.csv' for actual site names\n")
cat("   → Look for similar names to your coral sites\n")
cat("   → Create manual matching table\n\n")

cat("3. Once matches are confirmed:\n")
cat("   → Use 'environmental_matrix_for_RDA.csv' for genomics\n\n")

if (n_unmatched > 0) {
  cat("⚠️  UNMATCHED SITES DETECTED\n")
  cat("=========================\n\n")
  cat("Some of your coral sites could not be automatically matched.\n")
  cat("Please check:\n")
  cat("  1. Are the site names spelled exactly as in the eReefs data?\n")
  cat("  2. Do the eReefs sites use different naming conventions?\n")
  cat("  3. Check 'all_available_sites.csv' for the actual names\n\n")
  cat("To manually match sites, create a file called 'manual_site_matching.csv':\n")
  cat("  coral_site,ereefs_site\n")
  cat("  CBHE,Aken_1191\n")
  cat("  CBLM,SomeOtherSite\n")
  cat("  ...\n\n")
}

# Save R objects
saveRDS(results_mean, file.path(output_dir, "results_mean_all_sites.rds"))
saveRDS(results_p95, file.path(output_dir, "results_p95_all_sites.rds"))
if (!is.null(env_matrix)) {
  saveRDS(env_matrix, file.path(output_dir, "environmental_matrix.rds"))
}

processing_log <- list(
  timestamp = Sys.time(),
  files = names(ereefs_files),
  n_total_sites = results_mean$processing_info$n_sites,
  n_months = results_mean$processing_info$n_months,
  date_range = results_mean$processing_info$date_range,
  coral_sites_requested = coral_sites,
  coral_sites_matched = n_matched,
  coral_sites_unmatched = n_unmatched,
  site_matching = site_matching
)

saveRDS(processing_log, file.path(output_dir, "processing_log.rds"))

cat("\n✓ Processing complete!\n")
cat("✓ All outputs saved to:", output_dir, "\n")
if (!is.null(env_matrix)) {
  cat("✓ Environmental matrix ready for genomics\n\n")
} else {
  cat("⚠️  Please resolve site matching issues before using in genomics\n\n")
}
