# ThermoSSM

**ThermoSSM** is an R package for state-space analysis of temperature time series.
It implements linear Gaussian state-space models and performs inference using
Kalman filtering and smoothing.

The core implementation is adapted from the supplementary code provided in
Baba et al. (2024), which is publicly available on GitHub:
https://github.com/logics-of-blue/sea-temperature-trend-jogashima

### Key features

- Designed for **monthly temperature time series data**
- Applies **linear Gaussian state-space models** to estimate latent states
  using Kalman filtering and smoothing
- Represents temperature dynamics as a sum of latent state components:
  long-term trend, seasonal cycle, autoregressive structure,
  and optional exogenous effects
- Implement of time series cross-validation (under development)

# How to Install
```r
if(!require("devtools"))
	install.packages("devtools")
devtools::install_github("akihirao/ThermoSSM")

# install with vignettes
devtools::install_github("akihirao/ThermoSSM", build_vignettes=TRUE)
```

# Quick Tutorial
```r
browseVignettes("ThermoSSM")
```
https://github.com/akihirao/ThermoSSM/blob/main/vignettes/quick_tutorial.md


# Example Data

The package includes three example datasets for demonstration and testing of state-space temperature analyses. Two datasets are provided in `.rda` format
and consist of publicly available temperature data released by the Japan Meteorological Agency (JMA) (https://www.jma.go.jp/jma/indexe.html).
```r
data(ibaraki_sst)  # zoo object of daily sea surface temperature (SST) off the southern coast of Ibaraki Prefecture, Japan
data(fuji_temp)    # ts object of monthly air temperature at the summit of Mt. Fuji, Japan
```
The third dataset is provided in .csv format and included in inst/extdata. It contains a monthly temperature time series measured at Mount Akadake, Hokkaido, Japan (elevation 1,840 m; 43.6766°N, 142.9423°E). The data originate from the Monitoring Sites 1000 Project conducted by the Ministry of the Environment of Japan (KOZ01.zip, downladed from https://www.biodic.go.jp/moni1000/findings/data/index.html).
```r
path <- system.file("extdata", "example_monthly_temp.csv", package = "ThermoSSM")
readr::read_csv(path)
```

# References
Baba, S., Ishii, H., and Yoshiyama, T. (2024).  
Estimating sea temperature trends using a linear Gaussian state-space model in Jogashima, Kanagawa, Japan.
*Bulletin of the Japanese Society of Fisheries Oceanography*, 88(3), 190–199 (in Japanese with English abstruct)
https://doi.org/10.34423/jsfo.88.3_190

Baba, S. (2024).  
Supplementary code and test data for estimating sea temperature trends using a linear Gaussian state-space model.
GitHub repository:  
https://github.com/logics-of-blue/sea-temperature-trend-jogashima

