# ==============================================================================
# 00_config.R
# Purpose: Centralized library loading and file path definitions for workflow.
# Usage: Source this file at the start of every script to ensure consistency.
# ==============================================================================

# --- 1. LIBRARY MANAGEMENT ----------------------------------------------------
# Load necessary packages using pacman for efficient installation and loading
if (!require("pacman")) {
    install.packages("pacman")
}
pacman::p_load(
    "terra",
    "sf",
    "dplyr",
    "purrr",
    "tidyr",
    "readr",
    "stringr",
    "ggplot2",
    "tictoc",
    "exactextractr"
)

# --- 2. BASE DIRECTORIES & CRS ------------------------------------------------
# Define core folders
INPUT_DIR <- "data/raw"
DERIVED_DIR <- "data/derived"

# Coordinate Reference Systems (CRS)
ALBERS_CRS <- "EPSG:5070" # NAD83 / Conus Albers
WGS_CRS <- "EPSG:4326" # WGS 84

# --- 3. ROUTING & EXTENT LOGIC ------------------------------------------------

# Failsafes if run outside of main.R
if (!exists("PROCESSING_EXTENT")) {
    PROCESSING_EXTENT <- "NEBRASKA"
}
if (!exists("TARGET_LRR")) {
    TARGET_LRR <- NULL
}

if (PROCESSING_EXTENT == "NEBRASKA") {
    ACTIVE_MLRA_PATH <- file.path(DERIVED_DIR, "mlra/Nebraska_MLRA.gpkg")
    ACTIVE_GRID_PATH <- file.path(DERIVED_DIR, "grids/Nebraska_1km_mlra.gpkg")

    # Hardcoded Nebraska MLRAs
    ALL_MLRA_IDS <- c(
        "72",
        "77",
        "78",
        "79",
        "80",
        "81",
        "86",
        "87",
        "88",
        "89",
        "90",
        "142",
        "144A"
    )

    message("--> Config Loaded: Operating in NEBRASKA mode")
} else if (PROCESSING_EXTENT == "GREAT_PLAINS") {
    ACTIVE_MLRA_PATH <- file.path(DERIVED_DIR, "mlra/lower48MLRA.gpkg")
    ACTIVE_GRID_PATH <- file.path(
        DERIVED_DIR,
        "grids/GreatPlains_1km_mlra.gpkg"
    )

    # Dynamically fetch MLRA IDs based on the Target LRRs
    mlra_attrs <- sf::st_drop_geometry(sf::st_read(
        ACTIVE_MLRA_PATH,
        quiet = TRUE
    ))

    # Auto-detect the LRR column name (usually LRRSYM or LRR)
    lrr_col <- intersect(
        c("LRRSYM", "LRR_SYM", "LRR", "lrr", "lrr_sym"),
        names(mlra_attrs)
    )[1]

    if (!is.na(lrr_col) && !is.null(TARGET_LRR)) {
        ALL_MLRA_IDS <- unique(mlra_attrs$MLRA_ID[
            mlra_attrs[[lrr_col]] %in% TARGET_LRR
        ])
        message(sprintf(
            "--> Config Loaded: GREAT PLAINS mode. Filtered %d MLRAs inside LRRs: %s",
            length(ALL_MLRA_IDS),
            paste(TARGET_LRR, collapse = ", ")
        ))
    } else {
        ALL_MLRA_IDS <- unique(mlra_attrs$MLRA_ID)
        message("--> Config Loaded: GREAT PLAINS mode (No LRR filter applied).")
    }
} else {
    stop("PROCESSING_EXTENT must be either 'NEBRASKA' or 'GREAT_PLAINS'")
}

# If you ever want to test a single MLRA, you can override this here
TARGET_MLRA_IDS <- ALL_MLRA_IDS

# --- 4. STATIC INPUTS DIRECTORY -----------------------------------------------
STATIC_INPUTS <- list(
    # Grids
    grid_12mile = file.path(INPUT_DIR, "grid12M/twelve_mi_grid_uid.gpkg"),
    grid_1km = ACTIVE_GRID_PATH,

    # MLRA Boundaries
    mlra = ACTIVE_MLRA_PATH,

    # Riparian (Ensure this is not NULL)
    riparian = file.path(INPUT_DIR, "riparianArea/riparianArea10.tif"),

    # NLCD List
    nlcd = list(
        y2010 = file.path(
            INPUT_DIR,
            "nlcd/Annual_NLCD_LndCov_2010_CU_C1V1.tif"
        ),
        y2016 = file.path(
            INPUT_DIR,
            "nlcd/Annual_NLCD_LndCov_2016_CU_C1V1.tif"
        ),
        y2020 = file.path(INPUT_DIR, "nlcd/Annual_NLCD_LndCov_2020_CU_C1V1.tif")
    )
)

# --- 5. VALIDATION CHECK ---------------------------------------------------------
# Ensure inputs are not NULL
if (is.null(STATIC_INPUTS$riparian)) {
    stop("Error: STATIC_INPUTS$riparian is NULL. Check config.")
}
if (is.null(STATIC_INPUTS$nlcd$y2010)) {
    stop("Error: STATIC_INPUTS$nlcd$y2010 is NULL. Check config.")
}
