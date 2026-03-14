# Neyman Allocation TOF Sampling Pipeline

An R-based geospatial pipeline that determines optimal sample budgets for monitoring **Tree of Fire (TOF)** coverage across Multiple Land Resource Areas (MLRAs) using **Neyman optimal allocation** for stratified sampling.

---

## Overview

Simple Random Sampling (SRS) applied to heterogeneous landscapes requires large sample sizes to achieve acceptable accuracy. This pipeline replaces SRS with **stratified sampling** using Neyman allocation, which concentrates sampling effort in strata with high variance — reducing the required sample budget by up to 50–70% while maintaining the same accuracy targets.

The pipeline:

1. **Generates** a 1 km grid overlay for each MLRA.
2. **Extracts** static land cover attributes (NLCD classes, riparian percentage) per grid cell.
3. **Links** each grid cell to a TOF classification raster (from either a Random Forest or UNET model).
4. **Tests** hundreds of stratification strategies (variable × method × k combinations) and scores each by Neyman efficiency gain over SRS.
5. **Selects** a universal stratification strategy that works well across all MLRAs.
6. **Generates** variance profiles that describe stratum-level standard deviations for each Neyman group.
7. **Predicts** sample allocations — the number of samples per stratum needed to hit an accuracy target — for any MLRA.
8. **Validates** the predicted allocations through Monte Carlo simulation.

---

## Directory Structure

```
neymanAllocation/
├── main.r                          # Master execution script; toggle steps on/off here
├── neymanSampling.Rproj            # RStudio project file
├── requirements.txt                # (placeholder)
├── data/
│   ├── raw/                        # Input data (not tracked in git)
│   │   ├── grid12M/                # 12-mile tile grid (twelve_mi_grid_uid.gpkg)
│   │   ├── mlra/                   # MLRA boundary shapefiles
│   │   ├── nlcd/                   # Annual NLCD land cover rasters (2010, 2016, 2020)
│   │   ├── riparianArea/           # Riparian area raster (riparianArea10.tif)
│   │   ├── cot_10meter/            # Random Forest TOF classification tiles
│   │   └── cot_unet/               # UNET TOF classification tiles
│   ├── derived/                    # Intermediate outputs written by each pipeline step
│   │   ├── grids/                  # 1 km MLRA grids (gpkg)
│   │   ├── mlra/                   # Processed MLRA boundary files
│   │   ├── static_attributes/      # Per-grid NLCD + riparian percentages (CSV)
│   │   ├── dynamic_attributes/     # Per-grid TOF values joined to static data (CSV)
│   │   ├── dynamic_attributes_unet/
│   │   ├── method_testing/         # Stratification benchmark scores per MLRA (CSV)
│   │   ├── method_testing_unet/
│   │   ├── simulation_results/     # SRS baseline results from 100-iteration simulation
│   │   ├── simulation_results_unet/
│   │   ├── projected_sampleEstimates/  # Neyman sample budget estimates (CSV)
│   │   ├── projected_sampleEstimates_unet/
│   │   ├── variance_profiling/     # Universal variance profiles + MLRA group assignments
│   │   ├── variance_profiling_unet/
│   │   ├── final_allocations/      # Per-stratum sampling guides (CSV)
│   │   ├── final_allocations_unet/
│   │   ├── validation_results/     # Simulation validation reports (CSV)
│   │   └── novel_validation/       # Results from out-of-sample MLRA validation
│   └── temp/                       # Scratch space
├── scripts/
│   ├── 00_config.R                 # Library loading, path definitions, CRS constants
│   ├── 00a_generate_1kmGrids.R     # One-time grid generation from 100 km parent grids
│   ├── 01_static_processing.R      # Extract NLCD class percentages + riparian per cell
│   ├── 02_dynamic_processing.R     # Extract TOF (Random Forest model) per grid cell
│   ├── 02_dynamic_processing_unet.R# Extract TOF (UNET model) per grid cell
│   ├── 03_stratification_testing.R # Score 300+ stratification strategies per MLRA
│   ├── 04_random_sampling_test.R   # SRS baseline: 100-iteration simulation
│   ├── 05_gatherResults.R          # Aggregate SRS + Neyman results into sample budgets
│   ├── 06_mlra_grouping.R          # Cluster MLRAs by dominant Neyman stratification class
│   ├── 07_variance_profiling.R     # Build universal stratum-level variance profiles
│   ├── 08_predict_allocations.R    # Apply Neyman allocation formula to produce sampling guides
│   ├── 09_test_predicted_allocations.R  # 100-iteration validation of predicted budgets
│   ├── 10_novel_mlra_validation.R  # End-to-end validation on a never-seen MLRA
│   └── utils/
│       ├── generate_summary_stats.R
│       ├── comparingUNETandRFresults.R
│       ├── summarizeMLRA_landCover_Area.R
│       └── tempMapMLRA.R
├── src/
│   ├── neymanHelperFunctions.R     # Core reusable functions (stratification, stats, plots)
│   └── sampleGridsFunctions.R      # Grid generation utilities
└── outputs/
    ├── figures/                    # Saved ggplot visualizations
    ├── models/                     # Model output files
    └── reports/                    # Summary reports
```

---

## Pipeline Steps

The pipeline is controlled from `main.r` by toggling boolean flags. Steps must be run in order because each step consumes outputs from prior steps.

```
Step 00a  →  Step 01  →  Step 02  →  Step 03  →  Step 04
(grids)    (static)   (dynamic)  (strat test) (SRS baseline)
                                                     ↓
                                               Step 05 (gather results)
                                                     ↓
                                      Step 06  →  Step 07
                                   (MLRA groups) (variance profiles)
                                                     ↓
                                      Step 08  →  Step 09  →  Step 10
                                   (allocations) (validation) (novel MLRA)
```

| Step | Script | Inputs | Outputs |
|------|--------|--------|---------|
| 00a | `00a_generate_1kmGrids.R` | `grid100km_aea.gpkg`, MLRA boundaries | `data/derived/grids/*_1km_mlra.gpkg` |
| 01 | `01_static_processing.R` | 1 km grids, NLCD rasters, riparian raster | `static_attributes/MLRA_*_static_attributes.csv` |
| 02 | `02_dynamic_processing.R` | Static attributes, TOF tiles | `dynamic_attributes/MLRA_*_master_dataset.csv` |
| 03 | `03_stratification_testing.R` | Master dataset CSVs | `method_testing/MLRA_*_all_methods_comparison.csv` |
| 04 | `04_random_sampling_test.R` | Master dataset CSVs | `simulation_results/MLRA_*_baseline_weighted_results.csv` |
| 05 | `05_gatherResults.R` | Steps 03 + 04 outputs | `projected_sampleEstimates/MLRA_sample_estimates.csv` |
| 06 | `06_mlra_grouping.R` | Step 05 outputs | `variance_profiling/MLRA_Neyman_Groups.csv` |
| 07 | `07_variance_profiling.R` | Steps 02 + 06 outputs | `variance_profiling/Universal_Variance_Profiles.csv` |
| 08 | `08_predict_allocations.R` | Steps 06 + 07 outputs | `final_allocations/MLRA_*_predicted_allocations.csv` |
| 09 | `09_test_predicted_allocations.R` | Steps 02 + 08 outputs | `validation_results/Final_Validation_Report.csv` |
| 10 | `10_novel_mlra_validation.R` | Steps 06 + 07 outputs + raw tiles | `novel_validation/MLRA_*_Novel_Validation_Report.csv` |

---

## Setup and Configuration

### Prerequisites

- **R** (≥ 4.0)
- The [`pacman`](https://cran.r-project.org/package=pacman) package (installed automatically by `00_config.R` if missing)

Required packages are loaded via `pacman::p_load()` in `00_config.R`:

```r
terra, sf, dplyr, purrr, tidyr, readr, stringr, ggplot2, tictoc, exactextractr
```

### Configuration (`scripts/00_config.R`)

`00_config.R` is sourced automatically at the top of every script. It defines:

| Variable | Description |
|----------|-------------|
| `INPUT_DIR` | Path to raw input data (`data/raw`) |
| `DERIVED_DIR` | Path to intermediate outputs (`data/derived`) |
| `ALBERS_CRS` | Albers Equal Area projection (EPSG:5070) used for all spatial operations |
| `WGS_CRS` | WGS 84 (EPSG:4326) for geographic coordinates |
| `ALL_MLRA_IDS` | List of MLRAs to process, derived from `PROCESSING_EXTENT` |
| `STATIC_INPUTS` | Named list of paths to all input files |

### Running the Pipeline

1. Open `main.r` in RStudio and set the working directory:

   ```r
   setwd("path/to/neymanAllocation")
   ```

2. Choose the processing extent and target LRRs:

   ```r
   PROCESSING_EXTENT <- "GREAT_PLAINS"   # or "NEBRASKA"
   TARGET_LRR        <- c("F", "G", "H") # Land Resource Regions to include
   ```

3. Enable the steps you want to run by setting flags to `TRUE`:

   ```r
   RUN_00a_GENERATE_GRIDS      <- FALSE  # Run once only
   RUN_01_STATIC_PROCESSING    <- TRUE
   RUN_02_DYNAMIC_PROCESSING   <- FALSE
   RUN_03_STRATIFICATION_TESTING <- FALSE
   RUN_04_BASELINE_SAMPLING    <- FALSE
   RUN_05_GATHER_RESULTS       <- FALSE
   RUN_06_MLRA_GROUPING        <- FALSE
   RUN_07_VARIANCE_PROFILING   <- FALSE
   RUN_08_PREDICT_ALLOCATIONS  <- FALSE
   RUN_09_TEST_PREDICTIONS     <- FALSE
   RUN_10_NOVEL_VALIDATION     <- FALSE
   ```

4. Source `main.r`:

   ```r
   source("main.r")
   ```

Each step prints timing and progress messages. Intermediate CSV files are written after each step so the pipeline can be paused and resumed.

> **UNET vs. Random Forest:** Scripts 02–09 include a `USE_UNET` toggle at the top of each file. When `TRUE`, the script reads from `*_unet` input directories and writes to `*_unet` output directories, keeping the two model pipelines fully separated.

---

## Key Algorithms

### Neyman Optimal Allocation

Given a total sample budget *n* and *L* strata, the Neyman formula allocates samples to stratum *h* as:

```
n_h = n × (N_h × S_h) / Σ(N_i × S_i)
```

where *N_h* is the population size and *S_h* is the standard deviation of TOF within stratum *h*. Strata with greater variability receive proportionally more samples.

**Efficiency Gain** compared to SRS is calculated as:

```
Efficiency Gain (%) = (SRS_score - Neyman_score) / SRS_score × 100
```

where `score = Σ(N_h × S_h)` from `calculate_neyman_stats()`.

### Stratification Methods

Four stratification methods are tested in `03_stratification_testing.R`:

| Method | Description |
|--------|-------------|
| `quantile` | Divides the variable into *k* equal-frequency bins |
| `kmeans` | Clusters grid cells into *k* groups minimizing within-cluster variance |
| `zero_quantile` | Separates zero values into their own stratum; applies quantile bins to non-zeros |
| `zero_kmeans` | Separates zero values into their own stratum; applies k-means to non-zeros |

Each method is tested with *k* ∈ {3, 4, 5} strata across 9 stratification variables (`riparian_pct`, `Forest`, `Cultivated`, `Wetlands`, `Water`, `Developed`, `Barren`, `Shrubland`, `Herbaceous`) and 3 NLCD years (2010, 2016, 2020), yielding **324+ combinations per MLRA**.

### MLRA Grouping

`06_mlra_grouping.R` assigns each MLRA to a **Neyman Group** (Forest, Wetlands, Cultivated, or Mixed) based on which land cover variable drives the best stratification. This grouping allows variance profiles learned from well-sampled MLRAs to be transferred to data-sparse MLRAs.

### Universal Variance Profiles

`07_variance_profiling.R` computes stratum boundary thresholds and per-stratum standard deviations (*S_h*) for each Neyman Group. These profiles are stored in `Universal_Variance_Profiles.csv` and used in step 08 to allocate samples to any new MLRA without needing to re-run the full stratification search.

---

## Core Source Functions (`src/neymanHelperFunctions.R`)

| Function | Purpose |
|----------|---------|
| `apply_stratification(data, variable, method, k)` | Routes to the appropriate stratification algorithm and returns a dataframe with a `strata` column |
| `merge_small_strata(df, target_col, min_size)` | Iteratively merges the smallest stratum with its most similar neighbor until all strata meet a minimum cell count |
| `calculate_neyman_stats(df)` | Computes *N_h*, *S_h*, and *N_h × S_h* per stratum; returns the optimization metric Σ(*N_h* × *S_h*) |
| `combine_top_performers(use_unet, top_n, base_dir)` | Scans method-testing CSVs and returns the top *n* strategies per MLRA in a combined dataframe |
| `find_universal_strategy(use_unet, base_dir)` | Identifies the single stratification strategy with the highest average efficiency gain across all MLRAs |
| `analyze_stratified_weighted_sample(sample_df, pop_strata_sizes, N_total, area_col)` | Implements a combined ratio estimator with 95% CI (finite population correction) for area-weighted TOF mean estimation |
| `plot_top_performers(results_df, m_id, out_dir)` | Saves a bar chart of the top 10 strategies for one MLRA |
| `plot_strategy_comparison(df, use_unet)` | Generates a faceted lollipop chart comparing top strategies across all MLRAs |
| `plot_universal_vs_local(comparison_df, use_unet)` | Generates a dumbbell chart showing the efficiency cost of standardization vs. local optimization |

### Grid Functions (`src/sampleGridsFunctions.R`)

| Function | Purpose |
|----------|---------|
| `buildGrids(extent_object, cell_size)` | Creates a regular grid (Albers EPSG:5070) over the given extent |
| `buildSubGrids(grids, cell_size, aoi)` | Subdivides parent grid cells and clips to an area of interest |
| `cropGridsMLRA(mlra_Group, mlra, grids_84)` | Clips grid cells to MLRA boundaries |
| `cropGridsLRR(lrr_Group, lrr, grids_84)` | Clips grid cells to Land Resource Region boundaries |

---

## Data Flow Summary

```
Raw Inputs (data/raw/)
  NLCD rasters + Riparian raster
  MLRA boundaries + 12-mile tile grid
  TOF classification tiles (RF or UNET)
          │
          ▼
  Step 01: static_attributes (NLCD %, riparian % per 1 km cell)
          │
          ▼
  Step 02: dynamic_attributes (TOF % per 1 km cell)
          │
          ├──────────────────────────────────┐
          ▼                                  ▼
  Step 03: method_testing             Step 04: simulation_results
  (300+ strat. strategies scored)     (SRS baseline, 100 iterations)
          │                                  │
          └──────────────┬───────────────────┘
                         ▼
               Step 05: projected_sampleEstimates
               (Neyman vs. SRS budget comparison)
                         │
                         ▼
               Step 06: MLRA_Neyman_Groups
               (MLRA → Forest/Wetlands/Cultivated/Mixed)
                         │
                         ▼
               Step 07: Universal_Variance_Profiles
               (Sh per stratum per Neyman Group)
                         │
                         ▼
               Step 08: final_allocations
               (n_h per stratum per MLRA)
                         │
               ┌─────────┴─────────┐
               ▼                   ▼
       Step 09: validation    Step 10: novel MLRA
       (100× simulation)      (out-of-sample test)
```
