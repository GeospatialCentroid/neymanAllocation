library(dplyr)
library(readr)

#' Determine Optimal Land Cover Threshold for Sampling Budgets
#'
#' @param req_df Dataframe containing the sample requirements (e.g., sampleRequirementsNebraska.csv)
#' @param stats_df Dataframe containing the land cover statistics (e.g., ALL_MLRA_UNET_summary_stats.csv)
#' @param target_var The land cover column to test (e.g., "Forest", "Herbaceous")
#' @param direction The logical direction for the 'High' bucket ("<=" or ">=")
#' @param high_budget The static budget assigned to MLRAs exceeding the threshold (default: 1400)
#' @return A list containing the optimal threshold, the dynamic low budget, and the minimum RMSE score.
determine_optimal_threshold <- function(req_df, stats_df, target_var, direction = "<=", high_budget = 1400) {
  
  # 1. Prepare Requirements Data (Filter for Systematic & 95% Confidence)
  sys_req <- req_df %>%
    dplyr::filter(Method == "Systematic", Milestone == "3. 95% Acc + 95% CI") %>%
    dplyr::group_by(MLRA) %>%
    # Use max required N across years for safety
    dplyr::summarize(Required_N = max(Required_N, na.rm = TRUE), .groups = "drop")
  
  # 2. Prepare Stats Data (Average land cover across years per MLRA)
  stats_mean <- stats_df %>%
    dplyr::group_by(MLRA_ID) %>%
    dplyr::summarize(target_val = mean(.data[[target_var]], na.rm = TRUE), .groups = "drop")
  
  # 3. Merge Requirements and Statistics
  df <- sys_req %>%
    dplyr::inner_join(stats_mean, by = c("MLRA" = "MLRA_ID"))
  
  # 4. Exhaustive Search Setup
  best_rmse <- Inf
  best_thresh <- NA
  best_low_budget <- NA
  
  # Generate 1000 potential thresholds between the min and max observed values
  thresholds <- seq(min(df$target_val, na.rm = TRUE), 
                    max(df$target_val, na.rm = TRUE), 
                    length.out = 1000)
  
  for (thresh in thresholds) {
    
    # Apply threshold direction logic
    if (direction == ">=") {
      high_mask <- df$target_val >= thresh
    } else {
      high_mask <- df$target_val <= thresh
    }
    
    low_mask <- !high_mask
    
    # Ensure enough data remains to reliably calculate Standard Deviation
    if (sum(low_mask) < 2 || sum(high_mask) < 1) next
    
    # 5. Dynamic Lower Budget Rule (Mean + Standard Deviation of actual requirements)
    actual_low <- df$Required_N[low_mask]
    lower_budget <- mean(actual_low) + sd(actual_low)
    
    # Predict assignments
    assigned <- ifelse(high_mask, high_budget, lower_budget)
    actual <- df$Required_N
    
    # 6. Score via Root Mean Squared Error (RMSE)
    rmse <- sqrt(mean((assigned - actual)^2))
    
    if (rmse < best_rmse) {
      best_rmse <- rmse
      best_thresh <- thresh
      best_low_budget <- lower_budget
    }
  }
  
  return(list(
    Target_Variable = target_var,
    Direction = direction,
    Optimal_Threshold = best_thresh,
    Dynamic_Low_Budget = best_low_budget,
    RMSE = best_rmse
  ))
}

# Load your datasets
req_data <- read_csv("temp/sampleRequirementsNebraska.csv")
stats_data <- read_csv("~/trueNAS/work/neymanSampling/data/derived/summaries_unet")

# Run the optimization for Forest (looking for values LESS THAN the threshold)
forest_opt <- determine_optimal_threshold(
  req_data, 
  stats_data, 
  target_var = "Forest", 
  direction = "<="
)

cat(sprintf("\nOptimal %s Threshold: %s %.3f%%\n-> Dynamic Low Budget: %.1f samples\n-> RMSE: %.1f\n", 
            forest_opt$Target_Variable, 
            forest_opt$Direction,
            forest_opt$Optimal_Threshold, 
            forest_opt$Dynamic_Low_Budget, 
            forest_opt$RMSE))

# Run the optimization for Herbaceous (looking for values GREATER THAN the threshold)
herb_opt <- determine_optimal_threshold(
  req_data, 
  stats_data, 
  target_var = "Herbaceous", 
  direction = ">="
)

cat(sprintf("\nOptimal %s Threshold: %s %.3f%%\n-> Dynamic Low Budget: %.1f samples\n-> RMSE: %.1f\n", 
            herb_opt$Target_Variable, 
            herb_opt$Direction,
            herb_opt$Optimal_Threshold, 
            herb_opt$Dynamic_Low_Budget, 
            herb_opt$RMSE))

