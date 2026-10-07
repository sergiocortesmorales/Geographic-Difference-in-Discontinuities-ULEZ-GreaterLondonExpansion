# ULEZ 2023 Expansion and House Prices

R code to estimate the effect of the **London-wide ULEZ expansion** (29 August 2023) on house prices at the Greater London boundary, using a **geographic difference-in-discontinuities** design on **repeat sales**.

---

## Abstract

This project estimates how the 2023 expansion of London's Ultra Low Emission Zone capitalised into house prices near the Greater London boundary. The boundary separates a dense urban area from the Green Belt, so a cross-sectional RD or a difference-in-differences would be biased. Instead, a geographic difference-in-discontinuities is estimated on properties sold both before and after the expansion, which removes pre-existing jumps at the boundary and time-invariant property characteristics. The boundary effect is close to zero overall, but becomes positive and significant once properties within 100m of the boundary are excluded, and is larger for homes near rail and metro stations. There is no clear link between the price effect and measured changes in NO2 or PM2.5.

---

## Approach

- **Outcome**: change in log price for properties sold both before and after the expansion.
- **Running variable**: signed distance to the Greater London boundary in metres (positive = inside).
- **Estimation**: local polynomial RD with `rdrobust`, triangular kernel, forced bandwidths, clustered by LSOA.

---

## Scripts

| Script                              | Description                                                                                         |
| ----------------------------------- | --------------------------------------------------------------------------------------------------- |
| `scripts/1-clean.R`                 | Cleans Land Registry sales (2013-2025), geocodes them with ONSPD postcode centroids, keeps southern England regions and joins DfT journey times. |
| `scripts/2-spatialtransformation.R` | Computes distance to the boundary, treatment and period variables, and attaches Defra PCM pollution and distance to the nearest station. |
| `scripts/3-preanalysis.R`           | Descriptive plots and maps of prices, volumes and pollution near the boundary.                      |
| `scripts/4-estimation.R`            | Main estimation on repeat sales, with robustness checks. Exports tables to Word.                    |
| `scripts/5-estimation2.R`           | Same estimation on all sales (repeated cross-section).                                              |

Scripts 1 and 2 save intermediate `.rds` files to `data/`, so scripts 3 to 5 can be rerun without rebuilding the data.

---

## Setup

Packages are installed and loaded with `pacman` at the top of each script.

Raw data is not included. Download it from the sources below and place it as follows:

```
rawdata/
  pp-complete.csv
  ONSPD_NOV_2025_UK.csv
  jts0501.ods
  mapno22022.csv
  mapno22024.csv
  mappm252022g.csv
  mappm252024g.csv
  Stops.csv
data/
  Greater_London_Authority_(GLA).shp   (with its .shx, .dbf and .prj)
```

Two files need editing by hand first:

- `jts0501.ods`: variable names in the first row, 2017 sheet only.
- PCM csvs: metadata rows removed, keeping the 4 columns `gridcode`, `x`, `y` and the concentration.

---

## Usage

1. Clone this repository.
2. Open `R_ULEZ.Rproj` in RStudio.
3. Run the scripts in `scripts/` in order, from `1-clean.R` to `5-estimation2.R`.

Figures and tables are written to `output/`. The complete Price Paid file is several GB, so `1-clean.R` needs plenty of RAM.

---

## Data sources

- **House prices**: [HM Land Registry Price Paid Data](https://www.gov.uk/government/statistical-data-sets/price-paid-data-downloads)
- **Postcode coordinates**: [ONS Postcode Directory](https://geoportal.statistics.gov.uk/), November 2025
- **Boundary**: [Greater London boundary, London Datastore](https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london)
- **Air quality**: [Defra Pollution Climate Mapping](https://uk-air.defra.gov.uk/data/pcm-data)
- **Stations**: [DfT NaPTAN](https://beta-naptan.dft.gov.uk/download)
- **Journey times**: [DfT Journey Time Statistics, JTS0501](https://www.gov.uk/government/statistical-data-sets/journey-time-statistics-data-tables-jts)

Contains HM Land Registry data © Crown copyright and database right 2026. This data is licensed under the Open Government Licence v3.0.