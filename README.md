# SJU-Hybrid · Summer Extreme (Linear / Statistical–Dynamical)

Statistical–dynamical (multiple linear regression, **MLR**) seasonal prediction of **Korean summer compound hot–dry extremes (CDHD)** over the July–August (JA) season.

This repository contains the research implementation of the **linear statistical–dynamical model** developed in the Korea Meteorological Administration (KMA) R&D project *“AI-based research on climate-extreme variability and prediction”* (Sejong University component). The model predicts a **compound dry–hot index (CDHD)** over the Korean Peninsula from physically-motivated, skill-screened predictors.

## Overview

Pure dynamical models have limited skill for regional seasonal extremes, while purely data-driven models struggle with the short observational record. This model takes a **statistical–dynamical (hybrid)** approach: a lean, interpretable linear regression built on physical predictors that are screened to those a dynamical multi-model ensemble forecasts skillfully.

- **Predictand:** JA compound hot–dry index **CDHD** over Korea, combining heatwave intensity and drought:
  `CDHD = ( log1p(HWD) / σ ) × ( − SPEI3 )`,
  where `HWD` = Korean heatwave-day count and `SPEI3` = standardized precipitation–evapotranspiration index.
- **Predictors:** teleconnection indices at seasonal lead; the final CDHD model is **lead-based** (observed early-summer precursors; no concurrent dynamical term was selected).
- **Estimator:** multiple linear regression (this repo). The non-linear (CatBoost) counterpart is in a sibling repository.

## Method (pipeline)

| step | script | description |
|---|---|---|
| 1 | `src/00_corr_leads.py` | Detrended correlation of the predictand against global climate fields at leads 0–3 months (teleconnection diagnosis). |
| 2 | `src/01_make_predictor_index.py` | Predictor-pool construction: (observation significance) ∩ (dynamical-ensemble forecast skill) mask → **8-connected clustering (≥100 grid cells, periodic in longitude)** → cosine-latitude area-mean index → greedy collinearity screening (pairwise \|Spearman\| ≤ 0.6) + Southern-Hemisphere exclusion. |
| 3 | `src/02_select_mlr.py` | Final predictor selection: **best-subset (≤3 predictors)** maximizing **leave-one-out cross-validated skill (ACC)**, with VIF ≤ 5 multicollinearity control. |
| 4 | `src/03_forecast_hybrid.py` | Forecast / cross-validated skill. (L0 terms, where selected, are supplied by a distribution-calibrated C3S-MME forecast; CDHD uses observed lead predictors.) |

## Selected predictors (final CDHD model)

- June precipitation over the Korean Peninsula (L1) · June 2-m temperature over the North Atlantic (L1) · April precipitation over the Arctic (L3)

The dominant predictor is **early-summer (June) Korean precipitation** (negative relationship: drier June → higher compound hot–dry risk).

## Repository structure

```
.
├── src/
│   ├── 00_corr_leads.py           # teleconnection correlation maps (predictand vs fields, leads 0–3)
│   ├── 01_make_predictor_index.py # predictor-pool construction (mask → cluster → index → screen)
│   ├── 02_select_mlr.py           # best-subset MLR selection (LOOCV + VIF)
│   └── 03_forecast_hybrid.py      # forecast / hybrid skill
├── data/                          # input data (not included — see data/README.md)
├── requirements.txt
└── README.md
```

## Data

Raw datasets are **not included** (size / licensing). Place inputs under `./data/` (or set the path at the top of each script). Required inputs and sources:

- **ERA5** monthly atmospheric fields — Copernicus Climate Data Store (CDS).
- **Heatwave-day index** `ERA5_JA_HW_index_1979-2024.nc` (ERA5-based), and **SPEI** (`SPEIbase` v2.11) — for the CDHD predictand.
- **C3S seasonal multi-model** hindcast/forecast (dynamical predictors) — CDS.

See `data/README.md` for the expected file list and format.

## Usage

```bash
pip install -r requirements.txt
# place inputs under ./data  (see data/README.md)
python src/00_corr_leads.py           # teleconnection correlations
python src/01_make_predictor_index.py # build predictor pool
python src/02_select_mlr.py CDHD      # final MLR selection
python src/03_forecast_hybrid.py      # forecast / skill
```

## Authors

Jee-Hoon Jeong (정지훈)¹, Gyu-Ri Lee (이규리)¹
¹ Sejong University, Republic of Korea

## Funding

This research was supported by the **Korea Meteorological Administration Research and Development Program** under Grant **RS-2025-02307979**.

## Citation

If you use this code, please cite the associated project report / manuscript (Jeong, Lee et al., 2026) and acknowledge the funding grant above.

## Legal Notice

This repository contains original research code developed for academic and scientific purposes.

Any use, reproduction, modification, distribution, or derivative work based on this code is strictly prohibited without prior written permission from the authors.

If this code is used in any form of academic, research, or commercial work, proper citation and explicit acknowledgment of all original authors and affiliated institutions is mandatory.

Unauthorized use or redistribution without permission may result in legal consequences.
