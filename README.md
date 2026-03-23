# InSAR Local-to-Global VLM Conversion with GNSS Integration

> Transform local InSAR line-of-sight velocities into a globally referenced vertical land motion (VLM) frame using MIDAS IGS20 GNSS data — with full uncertainty propagation.

---
Author: Mahmoud Reshadati (mahmoudreshadati@vt.edu)

## Overview

This repository provides a MATLAB workflow for model-based conversion of **local InSAR-derived vertical land motion (VLM)** into the **IGS20 global reference frame** by integrating sparse GNSS measurements from the [MIDAS dataset](https://geodesy.unr.edu). The pipeline handles everything from downloading global GNSS data to spatially interpolating station velocities over dense InSAR pixel grids and applying a rigorous polynomial correction.

The method is designed for large-scale InSAR datasets (millions of pixels) and includes memory-managed processing, outlier rejection, and covariance-based uncertainty propagation.

---

## Methodology

```
MIDAS IGS20 GNSS Velocities
        │
        ▼
 Regional Station Filtering
 (temporal + spatial criteria)
        │
        ▼
 IDW Oversampling onto InSAR Grid
 (inverse-distance × uncertainty weighting)
        │
        ▼
 Polynomial WLS Correction (Δ-field)
 (local InSAR → global IGS20 frame)
        │
        ▼
 Global VLM Map + Propagated Uncertainties
```

### Step 1 — GNSS Data Download (`Code0_GNSS_MidasDownloader.m`)
Downloads the latest MIDAS IGS20 velocity file from the University of Nevada, Reno (UNR) and parses all stations into a structured MATLAB table. Velocities are stored in **m/yr** for East, North, and Up components along with their uncertainties.

### Step 2 — Regional Station Filtering (`Code1_FindingRegionalStations.m` + `fnc_GNSS_near_Region.m`)
Filters the global GNSS dataset to retain only stations that satisfy:
- A minimum and maximum observation **duration** (years)
- A minimum acceptable **start year**
- A maximum **radial distance** from the region of interest (haversine-based)

### Step 3 — IDW Oversampling (`InSAR_OverSampling_IDW_optimized.m`)
Spatially interpolates sparse GNSS vertical velocities onto the dense InSAR pixel grid using **Inverse Distance Weighting (IDW)** combined with GNSS uncertainty weighting:

$$w_{ij} = \frac{1}{d_{ij}^p \cdot \sigma_i^2}, \qquad \hat{V}_j = \frac{\sum w_{ij} V_i}{\sum w_{ij}}, \qquad \hat{\sigma}_j = \sqrt{\sum \left(\frac{w_{ij}}{\sum w_{ij}}\right)^2 \sigma_i^2}$$

Processing is chunked to stay within a user-defined RAM budget, making it suitable for datasets with **>10⁶ pixels**.

### Step 4 — Local-to-Global Transformation (`InSAR_local2global_WLS_optimized.m`)
Fits a **weighted polynomial correction surface** (degree 1–3) to the difference field Δ = GNSS − InSAR using Weighted Least Squares (WLS), and applies it to every pixel:

$$V_{\text{global},i} = V_{\text{local},i} + \mathbf{B}_i \hat{\boldsymbol{\beta}}$$

Uncertainty propagation is complete and accounts for the covariance between the InSAR input and the fitted correction:

$$\text{Var}(V_{\text{global},i}) = \sigma^2_{\text{InSAR},i} + \text{Var}(\text{correction}_i) + 2\,\text{Cov}(V_{\text{local},i}\, \text{correction}_i)$$

---

## Repository Structure

```
├── Code0_GNSS_MidasDownloader.m          # Download and parse MIDAS IGS20 data
├── Code1_FindingRegionalStations.m       # Extract GNSS stations per region
├── Code2_ConvertLocal2Global.m           # Main conversion pipeline
├── fnc_GNSS_near_Region.m                # Function: spatial/temporal GNSS filtering
├── InSAR_OverSampling_IDW_optimized.m    # Function: IDW interpolation (chunked)
├── InSAR_local2global_WLS_optimized.m    # Function: WLS local-to-global correction
└── InSAR_haversine_distance_multi.m      # Utility: vectorized haversine distance
```

---

## Getting Started

### Prerequisites

- MATLAB R2020a or later
- Mapping Toolbox (for `latlon2local`)
- Active internet connection (for initial MIDAS download)

### Usage

**1. Download the global GNSS dataset**

```matlab
run('Code0_GNSS_MidasDownloader.m')
```

This creates `GlobalMaps/MIDAS_IGS20_Velocities.mat`.

**2. Extract regional GNSS stations**

Edit the region parameters in `Code1_FindingRegionalStations.m`:

```matlab
CityName             = {'NewYork', 'LosAngeles'};
Region_Reflat        = [40.7,  34.05];
Region_Reflon        = [-74.0, -118.24];
MaxDistanceRef_meter = [0.5e6, 0.5e6];   % 500 km radius
minYears_period      = [10, 10];
```

Then run:
```matlab
run('Code1_FindingRegionalStations.m')
```

**3. Run the full local-to-global conversion**

Edit `Code2_ConvertLocal2Global.m` to set your InSAR data paths and parameters:

```matlab
IDW_power        = 1;    % IDW distance power
PN_Degree        = 1;    % Polynomial correction degree (1, 2, or 3)
available_memory = 8;    % RAM budget in GB for chunked processing
```

Then run:
```matlab
run('Code2_ConvertLocal2Global.m')
```

---

## Input Data Format

| Variable | Source | Units | Description |
|---|---|---|---|
| `lat_GNSS`, `lon_GNSS` | MIDAS IGS20 | degrees | GNSS station coordinates |
| `Vu_GNSS`, `Su_GNSS` | MIDAS IGS20 | m/yr → mm/yr | Vertical velocity and uncertainty |
| `elpx_ll` | InSAR product | degrees / degrees | Pixel coordinates and incidence angle |
| `MVel`, `MVelq` | InSAR product | cm/yr | LOS velocity and uncertainty |

> InSAR LOS velocities are converted to vertical using the incidence angle: `Vvert = VLOS / cos(θ)`

---

## Outputs

All outputs are saved in **mm/yr** per city under the configured output directory:

| File | Contents |
|---|---|
| `{City}_GNSS_Oversampled_P{p}.mat` | IDW-interpolated GNSS surface at InSAR pixel locations |
| `{City}_GlobalInSAR_L{l}_P{p}.mat` | Globally corrected InSAR VLM + propagated uncertainties |

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `IDW_power` | `1` | Power for inverse-distance weighting (higher = more local) |
| `PN_Degree` | `1` | Polynomial surface degree for WLS correction |
| `available_memory` | `8` | RAM limit (GB) for chunked IDW interpolation |
| `minYears_period` | `10` | Minimum GNSS observation duration (years) |
| `MaxDistanceRef_meter` | `500,000` | Search radius around region center (m) |

---

## Reference

GNSS data source:
> Blewitt, G., et al. (2016). GNSS Velocity Field for the United States and Globally: The MIDAS Solution. *Journal of Geophysical Research: Solid Earth*.

---

## Author

**Mahmoud Reshadati**
*2025*
*mahmoudreshadati@vt.edu*

---
## How to Cite
If you use this software, please cite:

Reshadati, M. (2026). *InSAR-Local2Global* (Version 1.0.1) [MATLAB code]. Zenodo. https://doi.org/10.5281/zenodo.19187890

---
## License
This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.