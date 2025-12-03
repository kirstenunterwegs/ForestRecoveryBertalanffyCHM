# ForestRecoveryBertalanffyCHM
Analysis of recovery of forest height post disturbance in Bavaria


# Understanding the managament impact on forest height recovery dynamics in Central Europe 
This repository holds code for running the analysis that integrates a biologically grounded theoretical model of forest growth with remotely sensed canopy height data and a forest disturbance map to better understand post-disturbance forest canopy recovery. For a detailed description of the purpose, methodology, and results, see the following paper:

Krüger, K., Fischer, F.J., Stritih, A., Klemmt3. H.J., Seidl, R. (2026). Managing the future: Post-disturbance forest recovery across management types in Central Europe. Submitted to Forest Ecology and Management.

The data required to reproduce the analysis can be obtained via the following link: placeholder for the link

Most data needed to reproduce the results can be obtained from public and open access sources, which are indicated in the readme per data folder in the data repository. All other layers can be generated with the code provided in this repository. The forest ownership map used to classify forest management types cannot be made available due to data usage agreements with the respective regional authority. For reproducibility, please request the respective layers from the regional authority, while we provide links to those that are openly accessible.

## platforms

software used for development and implementation : R version 4.3.3

## Scripts:
Scripts are named in order (1_, etc.). 

- `01_synthetic_fit_bertalanffy`: Synthetic growth-model test; Checks whether the Bertalanffy model can reliably separate growth rate (g) and disturbance legacy (hOffset) using synthetic height data with noise.
- `02_compare_nDOM_CHM`: CHM vs. nDOM height comparison; Compares LiDAR-CHM and photogrammetric nDOM height estimates across forest types, elevation bins, and gap sizes using 2,400 paired samples.
- `03_filter_disturbance_patches`: Disturbance-patch preparation; Merges, filters, and masks disturbance patches across Bavaria, removes harvest disturbances and unsuitable areas, and generates patch polygons and buffer layers.
- `04_extract.chm_patch`: CHM extraction per disturbance patch;  Identifies and loads CHM tiles, masks them to disturbance patches and buffers, and extracts annual canopy height layers for all patches in 2017–2023.
- `05_prepare_auxiliary_data`: Prepares ownership, management, topography, climate, and soil layers needed for later CHM extraction and modelling.
- `06_extract_recovery_spatial_subsetting`: Extracts CHM height trajectories, links them with environmental/management data, aggregates to Landsat scale, cleans data, filters severity, and creates modelling subsets. 
- `07.1_model_fit_asym_laplace_full`: Fits the main hierarchical Bayesian canopy-height recovery models (brms).
- `07.2_model_fit_student_full`: Fits the main hierarchical Bayesian canopy-height recovery models (brms) with a student error distribution.
- `07.07.3_model_fit_asym_laplace_sub500`: Fits the main hierarchical Bayesian canopy-height recovery models (brms) with a different spatial subset on a 500 m grid.
- `07.07.4_model_fit_asym_laplace_sub1000`: Fits the main hierarchical Bayesian canopy-height recovery models (brms) with a different spatial subset on a 1000 m grid.
- `08_validation`: Validates model predictions by comparing early- vs. late-year CHM heights for the same pixels; computes R², RMSE, bias, and produces diagnostic plots.
- `09_model_analysis`: Summarises fixed/random effects, visualises posterior distributions, simulates growth curves, compares management types, and exports tables/plots.
- `10_map_samples_CHMS`: Builds study-region maps, sample distributions, and example CHM panels to create overview figure.
