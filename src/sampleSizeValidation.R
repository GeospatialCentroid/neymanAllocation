
simulate_sample_sizes <- function(detailed_df, min_sample = 100, step_size = 100, iterations = 50, tolerance_pct = 0.10) {
  
  message("Starting sample size simulation using Spatially Balanced Systematic Selection...")
  
  mlra_year_groups <- detailed_df %>%
    filter(!is.na(Rep_Year)) %>%
    group_split(MLRA_ID, Rep_Year)
  
  simulation_results <- purrr::map_dfr(mlra_year_groups, function(grp) {
    
    current_mlra <- grp$MLRA_ID[1]
    current_year <- grp$Rep_Year[1]
    total_records <- nrow(grp)
    
    # 1. Establish the "Truth" for this specific MLRA/Year
    true_R <- sum(grp$Area_Value_1_M2, na.rm = TRUE) / sum(grp$Total_Area_M2, na.rm = TRUE)
    true_mean_pct <- true_R * 100
    tolerance_val_pct <- true_mean_pct * tolerance_pct
    
    sizes_to_test <- unique(c(
      total_records, 
      seq(from = total_records - (total_records %% step_size), to = min_sample, by = -step_size)
    ))
    
    purrr::map_dfr(sizes_to_test, function(target_size) {
      
      iters_to_run <- ifelse(target_size == total_records, 1, iterations)
      
      purrr::map_dfr(1:iters_to_run, function(i) {
        
        # --- SPATIALLY BALANCED SYSTEMATIC SELECTION ---
        # Sort the group by the hierarchical ID to maintain spatial order
        grp_sorted <- grp %>% arrange(AOI_Site_ID)
        
        # Calculate the fractional skip interval to ensure we get exactly 'target_size'
        skip_interval <- total_records / target_size
        
        # Pick a random starting point within the first interval
        start_index <- runif(1, min = 1, max = skip_interval)
        
        # Generate the sequence of evenly spaced indices
        target_indices <- round(seq(from = start_index, by = skip_interval, length.out = target_size))
        target_indices <- pmin(target_indices, total_records) # Safeguard against rounding overflow
        
        # Slice the dataframe using these evenly spaced indices
        sampled_grp <- grp_sorted %>% slice(target_indices)
        # -----------------------------------------------
        
        n <- target_size
        x <- sampled_grp$Total_Area_M2
        y <- sampled_grp$Area_Value_1_M2
        
        R_hat <- sum(y) / sum(x)
        est_tof_pct <- R_hat * 100
        
        fpc <- 1 - (n / total_records)
        x_bar <- mean(x)
        residuals <- y - (R_hat * x)
        
        # Calculate Variance and Standard Error 
        s2_res <- ifelse(n > 1, sum(residuals^2) / (n - 1), 0)
        se <- ifelse(n > 1, sqrt(fpc / n) * (1 / x_bar) * sqrt(s2_res), 0)
        
        # Convert SE to percentage scale for CI calculation
        se_pct <- se * 100
        ci_lower <- est_tof_pct - (1.96 * se_pct)
        ci_upper <- est_tof_pct + (1.96 * se_pct)
        
        # Evaluate Conditions
        is_accurate <- abs(est_tof_pct - true_mean_pct) <= tolerance_val_pct
        is_covered <- (ci_lower <= true_mean_pct) & (ci_upper >= true_mean_pct)
        
        data.frame(
          MLRA_ID = current_mlra,
          Rep_Year = current_year,
          Sample_Size = target_size,
          Iteration = i,
          True_TOF_Percent = true_mean_pct,
          Estimated_TOF_Percent = est_tof_pct,
          CI_95_Lower = ci_lower,
          CI_95_Upper = ci_upper,
          Is_Accurate = is_accurate,
          Is_Covered = is_covered,
          stringsAsFactors = FALSE
        )
      })
    })
  })
  
  message("Simulation complete.")
  return(simulation_results)
}

generate_tof_map <- function(spatial_sf, tof_df, target_year, target_size = "Full") {
  
  # 1. Parse and Aggregate the Data
  if (target_size == "Full") {
    map_data <- tof_df %>%
      filter(Rep_Year == target_year) %>%
      select(MLRA_ID, Plot_Value = Weighted_Avg_TOF_Percent)
  } else {
    map_data <- tof_df %>%
      filter(Rep_Year == target_year, Sample_Size == as.numeric(target_size)) %>%
      group_by(MLRA_ID) %>%
      summarize(Plot_Value = mean(Estimated_TOF_Percent, na.rm = TRUE), .groups = "drop")
  }
  
  if (nrow(map_data) == 0) {
    stop(sprintf("No data found for Year: %s, Sample Size: %s", target_year, target_size))
  }
  
  # 2. INNER JOIN to keep ONLY the MLRAs that have data values
  spatial_joined <- spatial_sf %>%
    mutate(MLRA_ID = as.character(MLRA_ID)) %>%
    inner_join(
      map_data %>% mutate(MLRA_ID = as.character(MLRA_ID)),
      by = "MLRA_ID"
    )
  
  # 3. Create the Unique Value Scale (Exact % TOF)
  spatial_joined <- spatial_joined %>%
    # Sort by value first so the legend is ordered logically
    arrange(Plot_Value) %>% 
    mutate(
      # Format the exact number to 2 decimal places and append the % sign
      Plot_Class = factor(sprintf("%.2f%%", Plot_Value), levels = unique(sprintf("%.2f%%", Plot_Value)))
    )
  
  # 4. Fetch and Crop State Boundaries
  options(tigris_use_cache = TRUE) 
  usa_states <- tigris::states(cb = TRUE, class = "sf") %>%
    st_transform(st_crs(spatial_joined)) 
  
  active_bbox <- st_bbox(spatial_joined)
  states_cropped <- st_crop(usa_states, active_bbox)
  
  # 5. Build the Map
  title_text <- sprintf("Weighted Average TOF - %s", target_year)
  subtitle_text <- ifelse(target_size == "Full",
                          "Full Target Dataset (~1400 features/MLRA)",
                          sprintf("Simulated Sample Size: %s features/MLRA", target_size))
  
  p <- ggplot() +
    # Layer 1: State boundaries
    geom_sf(data = states_cropped, fill = "grey95", color = "grey60", linewidth = 0.5) +
    
    # Layer 2: Active MLRAs (Mapped to exact % TOF values)
    geom_sf(data = spatial_joined, aes(fill = Plot_Class), color = "white", linewidth = 0.3) +
    
    # Discrete scale for the unique factor levels
    scale_fill_viridis_d(
      option = "mako", 
      name = "% TOF", # Title updated for clarity
      na.value = "grey90",
      direction = -1 # Reverses the color ramp so higher % is lighter/brighter
    ) +
    theme_void() + 
    labs(
      title = title_text,
      subtitle = subtitle_text,
      caption = "Source: Aggregated MLRA Modeling Output"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5),
      plot.caption = element_text(color = "grey50", size = 9),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
  
  return(p)
}


