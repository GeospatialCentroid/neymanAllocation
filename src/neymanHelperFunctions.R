# --- 2. HELPER FUNCTIONS ------------------------------------------------------

#' Merge Small Strata
#'
#' Iteratively merges the smallest stratum into its most similar neighbor
#' until all strata meet the minimum size threshold.
#'
#' @param df Dataframe containing 'strata' and the target variable column.
#' @param target_col Name of the column used for calculating means (e.g., 'Forest').
#' @param min_size Minimum number of grid cells required per stratum.
#' @return Dataframe with updated 'strata' column.
merge_small_strata <- function(df, target_col, min_size = 50) {
  # Loop until condition is met
  while (TRUE) {
    # 1. Calculate Stratum Sizes
    counts <- df %>%
      dplyr::count(strata) %>%
      dplyr::arrange(n)

    # 2. Check Exit Conditions
    if (nrow(counts) == 0) {
      return(df)
    } # Safety for empty df
    if (min(counts$n) >= min_size) {
      return(df)
    } # Smallest is big enough
    if (nrow(counts) <= 1) {
      return(df)
    } # Can't merge if only 1 left

    # 3. Identify the "Problem" Stratum (The smallest one)
    small_strat_id <- counts$strata[1]

    # 4. Find Best Merge Candidate
    # We want to merge with the neighbor (ID +/- 1) that has the closest Mean value.
    strat_means <- df %>%
      dplyr::group_by(strata) %>%
      dplyr::summarise(
        avg_val = mean(!!rlang::sym(target_col), na.rm = TRUE),
        .groups = "drop"
      )

    current_mean <- strat_means$avg_val[strat_means$strata == small_strat_id]

    # Identify potential neighbors (Previous and Next indices present in data)
    sorted_ids <- sort(strat_means$strata)
    idx <- which(sorted_ids == small_strat_id)

    neighbor_ids <- c()
    if (idx > 1) {
      neighbor_ids <- c(neighbor_ids, sorted_ids[idx - 1])
    }
    if (idx < length(sorted_ids)) {
      neighbor_ids <- c(neighbor_ids, sorted_ids[idx + 1])
    }

    neighbors <- strat_means %>% dplyr::filter(strata %in% neighbor_ids)

    if (nrow(neighbors) == 0) {
      return(df)
    }

    # Pick the neighbor with the absolute closest mean value
    target_strat_id <- neighbors$strata[which.min(abs(
      neighbors$avg_val - current_mean
    ))]

    # 5. Execute Merge
    df$strata[df$strata == small_strat_id] <- target_strat_id

    # 6. Re-Index Strata (1, 2, 3...)
    df$strata <- as.numeric(as.factor(df$strata))
  }
}

#' Calculate Neyman Statistics for a Stratification
calculate_neyman_stats <- function(df) {
  stats <- df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      Nh = n(),
      Sh = sd(TOF, na.rm = TRUE),
      Mean = mean(TOF, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Sh = ifelse(is.na(Sh), 0, Sh),
      Nh_Sh = Nh * Sh
    )

  optimization_metric <- sum(stats$Nh_Sh)

  return(data.frame(
    metric_score = optimization_metric,
    total_grids = sum(stats$Nh),
    n_strata = nrow(stats)
  ))
}

#' Generate Plot of Top 10 Performers
plot_top_performers <- function(results_df, m_id, out_dir) {
  # Filter top 10 per year
  top_10 <- results_df %>%
    dplyr::filter(method != "none") %>%
    dplyr::group_by(year) %>%
    dplyr::slice_max(order_by = efficiency_gain_pct, n = 10) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label = paste0(variable, " (", method, "-", k, ")"),
      # Reorder factor for plotting so bars are sorted
      label = factor(label, levels = unique(label[order(efficiency_gain_pct)]))
    )

  p <- ggplot2::ggplot(
    top_10,
    aes(x = efficiency_gain_pct, y = label, fill = as.factor(year))
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~year, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = paste("Top 10 Stratification Methods: MLRA", m_id),
      subtitle = "Ranked by Efficiency Gain (%) over Random Sampling",
      x = "Efficiency Gain (%)",
      y = "Stratification Strategy"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "lightgrey", color = NA)
    )

  ggsave(
    filename = file.path(out_dir, paste0("MLRA_", m_id, "_top_performers.png")),
    plot = p,
    width = 8,
    height = 10,
    bg = "white"
  )
}

#' Apply Stratification Logic
apply_stratification <- function(data, variable, method, k) {
  vals <- data[[variable]]

  # Safety: If variable is essentially constant/empty
  if (sum(vals, na.rm = TRUE) == 0 || var(vals, na.rm = TRUE) == 0) {
    return(NULL)
  }

  out_df <- data

  # --- METHOD 1: QUANTILE ---
  if (method == "quantile") {
    if (dplyr::n_distinct(vals) < k) {
      return(NULL)
    }
    tryCatch(
      {
        out_df$strata <- as.numeric(cut(
          vals,
          breaks = quantile(
            vals,
            probs = seq(0, 1, length.out = k + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE,
          labels = FALSE
        ))
      },
      error = function(e) return(NULL)
    )

    # --- METHOD 2: K-MEANS ---
  } else if (method == "kmeans") {
    tryCatch(
      {
        km <- kmeans(vals, centers = k, nstart = 10)
        rank_map <- rank(km$centers)
        out_df$strata <- rank_map[km$cluster]
      },
      error = function(e) return(NULL)
    )

    # --- METHOD 3: ZERO-INFLATED QUANTILE ---
  } else if (method == "zero_quantile") {
    if (k < 2) {
      return(NULL)
    }
    if (!any(vals == 0)) {
      return(apply_stratification(data, variable, "quantile", k))
    }

    tryCatch(
      {
        out_df$strata <- NA
        out_df$strata[vals == 0] <- 1

        non_zeros <- vals[vals > 0]
        if (length(non_zeros) > 0) {
          bins <- k - 1
          if (dplyr::n_distinct(non_zeros) >= bins) {
            nz_breaks <- quantile(
              non_zeros,
              probs = seq(0, 1, length.out = bins + 1),
              na.rm = TRUE
            )
            nz_strata <- as.numeric(cut(
              non_zeros,
              breaks = nz_breaks,
              include.lowest = TRUE,
              labels = FALSE
            ))
            out_df$strata[vals > 0] <- nz_strata + 1
          } else {
            return(NULL)
          }
        }
      },
      error = function(e) return(NULL)
    )

    # --- METHOD 4: ZERO-INFLATED K-MEANS ---
  } else if (method == "zero_kmeans") {
    if (k < 2) {
      return(NULL)
    }
    if (!any(vals == 0)) {
      return(apply_stratification(data, variable, "kmeans", k))
    }

    tryCatch(
      {
        out_df$strata <- NA
        out_df$strata[vals == 0] <- 1

        non_zeros <- vals[vals > 0]
        if (length(non_zeros) > 0) {
          bins <- k - 1
          # Run K-Means only on non-zeros
          km <- kmeans(non_zeros, centers = bins, nstart = 10)
          rank_map <- rank(km$centers) # 1..bins

          # Shift ranks by +1 so they start at 2
          out_df$strata[vals > 0] <- rank_map[km$cluster] + 1
        }
      },
      error = function(e) return(NULL)
    )
  }

  if ("strata" %in% names(out_df) && !any(is.na(out_df$strata))) {
    return(out_df)
  } else {
    return(NULL)
  }
}


library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

#' Combine Top Performing Strategies Across MLRAs
#'
#' Scans a directory for method comparison CSVs, extracts the top N performers
#' from each file based on Efficiency Gain, and combines them.
#'
#' @param use_unet Boolean. If TRUE, looks in the UNET directory for UNET files.
#' @param top_n Number of top rows to keep per file (default: 10).
#' @param base_dir Base derived data directory (default: "data/derived").
#' @return A single dataframe containing the best strategies for all MLRAs.
combine_top_performers <- function(
  use_unet = FALSE,
  top_n = 10,
  base_dir = "data/derived"
) {
  # 1. Set Paths Based on UNET Toggle
  if (use_unet) {
    input_dir <- file.path(base_dir, "method_testing_unet")
    file_pattern <- "_all_methods_comparison_unet.csv$"
  } else {
    input_dir <- file.path(base_dir, "method_testing")
    file_pattern <- "_all_methods_comparison.csv$"
  }

  # 2. Find Files
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)

  if (length(files) == 0) {
    warning(paste(
      "No files found in",
      input_dir,
      "matching pattern",
      file_pattern
    ))
    return(NULL)
  }

  message(paste("Found", length(files), "comparison files. Processing..."))

  # 3. Process Files
  combined_list <- lapply(files, function(f) {
    tryCatch(
      {
        df <- readr::read_csv(f, show_col_types = FALSE)

        if (!"Efficiency_Gain_Pct" %in% names(df)) {
          warning(paste(
            "Skipping",
            basename(f),
            "- Missing Efficiency_Gain_Pct column."
          ))
          return(NULL)
        }

        df_top <- df %>%
          dplyr::arrange(desc(Efficiency_Gain_Pct)) %>%
          head(top_n)

        return(df_top)
      },
      error = function(e) {
        warning(paste("Error reading", basename(f), ":", e$message))
        return(NULL)
      }
    )
  })

  # 4. Combine Results
  final_df <- dplyr::bind_rows(combined_list)

  message(paste("Combined dataframe created with", nrow(final_df), "rows."))
  return(final_df)
}

#' Find Universal Best Stratification Strategy
#'
#' Identifies the single strategy (Variable + Method + K) that yields the
#' highest average efficiency gain across ALL processed MLRAs.
#'
#' @param use_unet Boolean. If TRUE, looks in the UNET directory for UNET files.
#' @param base_dir Base derived data directory (default: "data/derived").
#' @return A list containing the Summary Table and the Detailed Comparison.
find_universal_strategy <- function(
  use_unet = FALSE,
  base_dir = "data/derived"
) {
  # 1. Set Paths Based on UNET Toggle
  if (use_unet) {
    input_dir <- file.path(base_dir, "method_testing_unet")
    file_pattern <- "_all_methods_comparison_unet.csv$"
  } else {
    input_dir <- file.path(base_dir, "method_testing")
    file_pattern <- "_all_methods_comparison.csv$"
  }

  # 2. Load ALL Data (Not just top 10)
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop(paste("No comparison files found in", input_dir))
  }

  full_df <- lapply(files, read_csv, show_col_types = FALSE) %>%
    bind_rows() %>%
    mutate(Strategy_ID = paste(Variable, Method, K, sep = " | "))

  # 3. Identify Valid Strategies (Must exist in ALL MLRAs)
  mlra_count <- n_distinct(full_df$MLRA)

  valid_strategies <- full_df %>%
    count(Strategy_ID) %>%
    filter(n == mlra_count) %>%
    pull(Strategy_ID)

  if (length(valid_strategies) == 0) {
    stop("No single strategy successfully ran across ALL MLRAs.")
  }

  common_df <- full_df %>% filter(Strategy_ID %in% valid_strategies)

  # 4. Calculate "Local Best" for Context
  local_best_df <- full_df %>%
    group_by(MLRA) %>%
    summarise(
      Local_Best_Gain = max(Efficiency_Gain_Pct),
      Local_Best_Strategy = Strategy_ID[which.max(Efficiency_Gain_Pct)],
      .groups = "drop"
    )

  # 5. Find the "Universal Best"
  strategy_stats <- common_df %>%
    group_by(Strategy_ID, Variable, Method, K) %>%
    summarise(
      Mean_Gain = mean(Efficiency_Gain_Pct),
      Min_Gain = min(Efficiency_Gain_Pct),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Gain))

  winner <- strategy_stats[1, ]
  winner_id <- winner$Strategy_ID

  # 6. Create Comparison Table (Winner vs Local Best)
  comparison_table <- common_df %>%
    filter(Strategy_ID == winner_id) %>%
    select(MLRA, Universal_Gain = Efficiency_Gain_Pct) %>%
    left_join(local_best_df, by = "MLRA") %>%
    mutate(Performance_Drop = Local_Best_Gain - Universal_Gain)

  # 7. Print Summary
  message(paste0("\n=== Universal Best Strategy ==="))
  message(paste("Strategy:  ", winner_id))
  message(paste("Avg Gain:  ", round(winner$Mean_Gain, 2), "%"))
  message(paste("Worst Case:", round(winner$Min_Gain, 2), "% (in one MLRA)"))

  message("\n=== Performance Drop (Regret) ===")
  avg_drop <- mean(comparison_table$Performance_Drop)
  message(paste0(
    "On average, choosing this universal strategy costs ",
    round(avg_drop, 2),
    "% efficiency compared to customizing per MLRA."
  ))

  return(list(
    Best_Strategy = winner,
    Comparison = comparison_table,
    All_Rankings = strategy_stats
  ))
}

#' Plot Universal vs. Local Strategy Comparison
#'
#' Generates a dumbbell chart showing the performance gap (regret) between
#' the best universal strategy and the absolute best local strategy for each MLRA.
#'
#' @param comparison_df The 'Comparison' dataframe returned by 'find_universal_strategy'.
#' @param use_unet Boolean. Updates the plot title to reflect if UNET data was used.
#' @return A ggplot object.
plot_universal_vs_local <- function(comparison_df, use_unet = FALSE) {
  plot_data <- comparison_df %>%
    pivot_longer(
      cols = c("Universal_Gain", "Local_Best_Gain"),
      names_to = "Strategy_Type",
      values_to = "Efficiency_Gain"
    ) %>%
    mutate(
      Strategy_Type = factor(
        Strategy_Type,
        levels = c("Local_Best_Gain", "Universal_Gain"),
        labels = c("Local Best (Custom)", "Universal (Standard)")
      )
    ) %>%
    arrange(desc(Performance_Drop)) %>%
    mutate(MLRA = factor(MLRA, levels = unique(MLRA)))

  # Dynamic title based on toggle
  plot_title <- ifelse(
    use_unet,
    "Cost of Standardization (UNET): Universal vs. Local",
    "Cost of Standardization: Universal vs. Local Optimization"
  )

  p <- ggplot(plot_data) +
    geom_segment(
      data = comparison_df,
      aes(
        x = Universal_Gain,
        xend = Local_Best_Gain,
        y = as.factor(MLRA),
        yend = as.factor(MLRA)
      ),
      color = "gray50",
      linewidth = 1
    ) +
    geom_point(
      aes(x = Efficiency_Gain, y = as.factor(MLRA), color = Strategy_Type),
      size = 4
    ) +
    scale_color_manual(
      values = c(
        "Local Best (Custom)" = "#E41A1C",
        "Universal (Standard)" = "#377EB8"
      )
    ) +
    labs(
      title = plot_title,
      subtitle = "Distance between dots represents the efficiency lost by using one standard strategy.",
      x = "Efficiency Gain (%)",
      y = "MLRA Region",
      color = "Strategy Type"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold")
    )

  return(p)
}


#' Plot Top Stratification Strategies
#'
#' Generates a faceted lollipop chart comparing the top performing
#' strategies for each MLRA. Sorts independently within each facet.
#'
#' @param df The dataframe created by 'combine_top_performers'.
#' @param use_unet Boolean. Updates the plot title to reflect if UNET data was used.
#' @return A ggplot object.
plot_strategy_comparison <- function(df, use_unet = FALSE) {
  # 1. Prepare Data for Plotting
  plot_data <- df %>%
    dplyr::mutate(
      Strategy_Label = paste0(Variable, " (", Method, ", k=", K, ")"),
      MLRA = as.factor(MLRA)
    ) %>%
    # Sort ascending so the highest value gets the highest factor level
    # (ggplot plots the highest factor level at the top of the y-axis)
    dplyr::arrange(MLRA, Efficiency_Gain_Pct) %>%
    dplyr::mutate(
      # Create a truly unique label per facet to prevent cross-facet interference
      Unique_Label = paste(Strategy_Label, MLRA, sep = "___"),
      # Lock in the factor levels in this perfectly sorted order
      Unique_Label = factor(Unique_Label, levels = Unique_Label)
    )

  # Dynamic title based on toggle
  plot_title <- ifelse(
    use_unet,
    "Top Stratification Strategies by MLRA (UNET)",
    "Top Stratification Strategies by MLRA"
  )

  # 2. Generate Plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Efficiency_Gain_Pct, y = Unique_Label)
  ) +

    # Draw the "Stick" of the lollipop
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = Efficiency_Gain_Pct,
        y = Unique_Label,
        yend = Unique_Label
      ),
      color = "gray60"
    ) +

    # Draw the "Pop" (Point)
    ggplot2::geom_point(ggplot2::aes(color = Variable), size = 4) +

    # Create separate panels for each MLRA
    ggplot2::facet_wrap(~MLRA, scales = "free_y") +

    # Strip the hidden '___MLRA' ID out of the labels so they look clean on the plot
    ggplot2::scale_y_discrete(labels = function(x) gsub("___.*$", "", x)) +

    # Styling and Labels
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = plot_title,
      subtitle = "Efficiency Gain (%) compared to Simple Random Sampling",
      x = "Efficiency Gain (%)",
      y = NULL,
      color = "Dominant Variable"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "gray95"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  return(p)
}


# --- 2. HELPER FUNCTION -------------------------------------------------------

#' Calculate Stratified Area-Weighted Mean and CI (Combined Ratio Estimator)
analyze_stratified_weighted_sample <- function(
  sample_df,
  pop_strata_sizes,
  N_total,
  area_col = "grid_area"
) {
  strata_stats <- sample_df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      n_h = dplyr::n(),
      sum_y_h = sum(TOF * !!rlang::sym(area_col)),
      sum_x_h = sum(!!rlang::sym(area_col)),
      s2_res_h = {
        # Prevent division by zero if sum_x_h is 0
        if (sum_x_h == 0) {
          0
        } else {
          r_h <- sum_y_h / sum_x_h
          y <- TOF * !!rlang::sym(area_col)
          x <- !!rlang::sym(area_col)
          res <- y - (r_h * x)
          ifelse(dplyr::n() > 1, sum(res^2) / (dplyr::n() - 1), 0)
        }
      },
      mean_x_h = mean(!!rlang::sym(area_col)),
      .groups = "drop"
    ) %>%
    dplyr::left_join(pop_strata_sizes, by = "strata")

  # Calculate point estimate
  numerator <- sum(
    (strata_stats$N_h / N_total) * (strata_stats$sum_y_h / strata_stats$n_h),
    na.rm = TRUE
  )
  X_bar_U <- sum(
    (strata_stats$N_h / N_total) * strata_stats$mean_x_h,
    na.rm = TRUE
  )

  if (X_bar_U == 0) {
    return(c(mean = NA, ci_l = NA, ci_u = NA))
  }

  R_hat_st <- numerator / X_bar_U

  # Calculate variance of the combined ratio estimator
  var_R_hat_st <- (1 / (X_bar_U^2)) *
    sum(
      ((strata_stats$N_h / N_total)^2) *
        (1 - (strata_stats$n_h / strata_stats$N_h)) *
        (strata_stats$s2_res_h / strata_stats$n_h),
      na.rm = TRUE
    )

  # Handle potential negative tiny variances due to floating point math
  var_R_hat_st <- max(0, var_R_hat_st)
  se <- sqrt(var_R_hat_st)

  # 95% Confidence Interval
  ci_lower <- R_hat_st - (1.96 * se)
  ci_upper <- R_hat_st + (1.96 * se)

  return(c(mean = R_hat_st, ci_l = ci_lower, ci_u = ci_upper))
}
