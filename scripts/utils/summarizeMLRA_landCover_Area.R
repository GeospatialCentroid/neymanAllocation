library(dplyr)

#' Calculate Area-Weighted Land Cover Percentages
#'
#' This function calculates the weighted average of land cover columns based on grid area.
#' It automatically groups by year (or other specified columns) to track changes over time.
#'
#' @param df The dataframe containing the dataset (e.g., Master Dataset).
#' @param area_col The name of the column containing grid area (default: "grid_area").
#' @param group_cols A vector of column names to group by (default: "year").
#' @param exclude_cols A vector of columns to exclude from the summary (e.g., IDs).
#' @return A summary dataframe with area-weighted averages for each land cover class.
summarize_weighted_land_cover <- function(df, 
                                          area_col = "grid_area", 
                                          group_cols = "year",
                                          exclude_cols = c("id", "MLRA_ID")) {
  
  # Identify columns to summarize (numeric columns excluding grouping, area, and IDs)
  target_cols <- df %>%
    select(where(is.numeric)) %>%
    select(-any_of(c(area_col, group_cols, exclude_cols))) %>%
    names()
  
  message(paste("Summarizing weighted averages for:", paste(target_cols, collapse = ", ")))
  
  # Perform weighted calculation
  summary_df <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(across(all_of(target_cols), 
                     ~ weighted.mean(., w = .data[[area_col]], na.rm = TRUE)),
              .groups = "drop")
  
  return(summary_df)
}

# --- Usage Example ---
allFiles <- list.files("data/derived/dynamic_attributes", full.names = TRUE)
master_dfs <- allFiles[grepl(pattern = "master_dataset.csv", x = allFiles)]
for(i in master_dfs){
  print(i)
  master_df <- read_csv(i)
  summary_table <- summarize_weighted_land_cover(master_df)
  print(summary_table)
}

