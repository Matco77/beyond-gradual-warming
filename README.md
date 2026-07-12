# Beyond Gradual Warming — Thesis code & data pipeline

Estimating weather→crop-yield and weather→household-energy **response functions** from historical
European data, to be replayed on ISIMIP and AMOC-stressed climate scenarios. MSc thesis, Bocconi.

## Repository layout
- `Amoc/Code/` — the R pipeline. Start with `PIPELINE.md` (methodology) and `REVIEW_GUIDE.md`.
  - `Dataset_Preparation.R`, `agriculturegeografical.R`, `cropweather_eobs_nuts3.R`,
    `crop_panel_nuts3.R`, `energyweather_eobs_country.R`, `energy_panel.R` — data prep.
  - `regression.R` — response-function estimation (crop at NUTS3, energy at country), plus
    robustness (wild-cluster bootstrap, two-way clustering, snapped-region drops).
  - `crop_long_difference.R` — Burke & Emerick long-run bound (crop).
- `Amoc/datasets/` — **prepared** datasets only (numbered pipeline outputs `0`–`16` + source CSVs).
- `Amoc/results/` — coefficient tables.
- `Amoc/Methodology Explanation/`, `Amoc/Presentations and Proposal/`, `Amoc/figures/`,
  `Amoc/Plots_Anomaly/`, `CropCalendar/` — write-up material.

## Data NOT in this repo
The raw climate grids (~58 GB: E-OBS, CMIP6 piControl, NAHosMIP, ISIMIP, GISCO population grid,
shapefiles) are excluded via `.gitignore` — too large for GitHub. They must be downloaded from
their original sources; the prep scripts document where. The prepared CSVs in `Amoc/datasets/`
are sufficient to re-run `regression.R` and `crop_long_difference.R`.

## Requirements
R with `fixest` and `data.table`. Run each script with `Rscript <file>`.
