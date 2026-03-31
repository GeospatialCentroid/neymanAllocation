# ==============================================================================
# 09_k5_vs_k3_comparison.R
# Purpose: Compares the efficiency gain lost when dropping from k=5 to k=3
#          for the top 'zero_kmeans' predictor in each MLRA.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)
library(tidyr)

# --- 1. SETTINGS --------------------------------------------------------------
# We use the UNET directory as the basis for this check
INPUT_DIR <- file.path(DERIVED_DIR, "method_testing_unet")
OUTPUT_DIR <- file.path(DERIVED_DIR, "presentation_visuals")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

files <- list.files(
  INPUT_DIR,
  pattern = "_all_methods_comparison_unet\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No method comparison files found in ", INPUT_DIR)
}

results_list <- list()

# --- 2. EXTRACTION LOOP -------------------------------------------------------
for (f in files) {
  df <- readr::read_csv(f, show_col_types = FALSE) %>%
    dplyr::filter(Method == "zero_kmeans")

  if (nrow(df) == 0) {
    next
  }

  m_id <- unique(df$MLRA)[1]

  # Find the Absolute Best Variable for this MLRA (usually at k=5)
  best_row <- df %>%
    dplyr::arrange(desc(Efficiency_Gain_Pct)) %>%
    dplyr::slice(1)

  top_var <- best_row$Variable

  # Filter the dataset to only look at this top variable for k=3 and k=5
  comparison <- df %>%
    dplyr::filter(Variable == top_var, K %in% c(3, 5)) %>%
    dplyr::select(MLRA, Variable, K, Efficiency_Gain_Pct) %>%
    tidyr::pivot_wider(
      names_from = K,
      names_prefix = "K_",
      values_from = Efficiency_Gain_Pct
    )

  # Ensure both K_3 and K_5 exist in the data before calculating
  if ("K_5" %in% names(comparison) && "K_3" %in% names(comparison)) {
    comparison <- comparison %>%
      dplyr::mutate(
        Efficiency_Loss_Pct = K_5 - K_3,
        Relative_Loss_Pct = (Efficiency_Loss_Pct / K_5) * 100
      ) %>%
      dplyr::rename(
        Gain_at_K5 = K_5,
        Gain_at_K3 = K_3
      )

    results_list[[length(results_list) + 1]] <- comparison
  }
}

# --- 3. EXPORT RESULTS --------------------------------------------------------
final_comparison <- dplyr::bind_rows(results_list) %>%
  dplyr::arrange(desc(Efficiency_Loss_Pct))

out_file <- file.path(OUTPUT_DIR, "UNET_K5_vs_K3_Efficiency_Loss.csv")
readr::write_csv(final_comparison, out_file)

message("Comparison complete. Results saved to: ", out_file)
print(final_comparison)
