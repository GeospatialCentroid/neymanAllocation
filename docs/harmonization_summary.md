# Summary: Addressing Temporal Variability in NAIP for Stable TOF Monitoring

## 1. The Problem: Atmospheric & Sensor Noise in NAIP
Analysis of multi-temporal NAIP imagery (2011–2021) reveals significant radiometric instability. Because NAIP is often collected across different days, sensors, and atmospheric conditions, the raw digital numbers (DN) for the same location can vary drastically year-over-year.

### Evidence:
- **Spectral Drift:** Kolmogorov-Smirnov (KS) tests on NDVI and Blue band distributions show statistically significant shifts between years for identical ground locations (see `scripts/33_evaluatingHistogramNormalization.R`).
- **Phantom Change:** When a static machine learning model is applied to these uncorrected images, the resulting Tree Outside Forest (TOF) area estimates show "phantom change." In some cases, we observed a **>20% difference in area** purely due to image quality, not actual land cover change.

## 2. The Solution: Radiometric Harmonization (ARD)
To address this, we implemented a **Histogram Normalization** workflow, often referred to as producing "Analysis Ready Data" (ARD). This process aligns the spectral distribution of "target" years to a "baseline" year (typically the year the model was trained on).

### Key Components:
- **Reference Year Selection:** Establishing a radiometrically "clean" baseline.
- **Band-wise Alignment:** Adjusting Red, Green, Blue, and NIR distributions to match the baseline.
- **Improved Model Transferability:** Ensuring the classifier "sees" consistent input data across time.

## 3. Results & Validation
Comparing classifications from raw vs. harmonized imagery demonstrates a marked improvement in stability:
- **Visual Consistency:** Harmonized RGB imagery matches the baseline's color profile, reducing haze and sensor-related artifacts (Visuals from `scripts/35_evaluatingHarmonizationOnModelOutputs.R`).
- **Classification Accuracy:** Binary canopy masks generated from harmonized data show much higher spatial overlap with baseline classifications than those from unnormalized data.
- **Noise Reduction:** Our upcoming aggregate analysis (Visual 6) aims to show a significant reduction in the standard deviation of area estimates across years for stable forest plots.

## 4. Collaborative Presentation Strategy
We propose presenting this work to collaborators in three phases:
1.  **Identify the Issue:** Show the "Phantom Change" histogram and spectral shifts to establish the need for correction.
2.  **Demonstrate the Method:** Use the side-by-side RGB/Classification maps to show the visual and structural "fix."
3.  **Quantify the Benefit:** Present the reduction in relative error across all monitored grids to prove global effectiveness.

---
*Prepared by Gemini CLI for Neyman Sampling Project*
