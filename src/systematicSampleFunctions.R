# --- 2. HELPER FUNCTIONS ------------------------------------------------------

# Systematic Sampling Implementation
draw_systematic_sample <- function(df, n_desired) {
  N_total <- nrow(df)

  # If we request more samples than population, return the whole population
  if (n_desired >= N_total) {
    return(df)
  }

  # Calculate the fractional step size
  k <- N_total / n_desired

  # Random initiation location between 1 and the step size
  start_idx <- runif(1, min = 1, max = k)

  # Generate the sequence of indices
  indices <- floor(start_idx + (0:(n_desired - 1)) * k)

  # Safety bounds (in case floating point math pushes the last index slightly over)
  indices[indices > N_total] <- N_total
  indices[indices < 1] <- 1

  return(df[indices, ])
}

analyze_weighted_sample <- function(
  sample_df,
  N_total,
  area_col = "grid_area"
) {
  n <- nrow(sample_df)
  x <- sample_df[[area_col]]
  y <- sample_df$TOF * x

  R_hat <- sum(y) / sum(x)
  fpc <- 1 - (n / N_total)
  x_bar <- mean(x)
  residuals <- y - (R_hat * x)
  s2_res <- sum(residuals^2) / (n - 1)
  se <- sqrt(fpc / n) * (1 / x_bar) * sqrt(s2_res)

  # 95% Confidence Interval
  ci_lower <- R_hat - (1.96 * se)
  ci_upper <- R_hat + (1.96 * se)

  return(c(mean = R_hat, ci_l = ci_lower, ci_u = ci_upper))
}


#' Decode Hex IDs and calculate a Z-Order (Morton) index for spatial sorting
add_spatial_sort_index <- function(df) {
  # Split the hex ID into X and Y character columns
  id_parts <- do.call(rbind, strsplit(df$id, "-"))

  # Convert Hex strings to numeric integers
  df$x_idx <- strtoi(id_parts[, 1], base = 16)
  df$y_idx <- strtoi(id_parts[, 2], base = 16)

  # Normalize to start at 0 (prevents integer overflow in R)
  x_norm <- df$x_idx - min(df$x_idx)
  y_norm <- df$y_idx - min(df$y_idx)

  # Interleave bits to create a Z-order curve index
  # (This translates the 2D grid into a spatially balanced 1D line)
  z_order <- rep(0, nrow(df))
  for (i in 0:14) {
    # 15 bits is safe for R's 32-bit signed integers
    z_order <- bitwOr(
      z_order,
      bitwOr(
        bitwShiftL(bitwAnd(x_norm, bitwShiftL(1, i)), i),
        bitwShiftL(bitwAnd(y_norm, bitwShiftL(1, i)), i + 1)
      )
    )
  }

  df$z_order <- z_order
  return(df)
}
