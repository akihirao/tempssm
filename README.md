# ThermoSSM

- R package for state-space analysis of temperature time series.  
The implementation is based on linear Gaussian state-space models and makes use of Kalman filtering and smoothing.  
Main parts of the implementation were adapted from the supplementary code provided in Baba et al. (2024), available on GitHub: https://github.com/logics-of-blue/sea-temperature-trend-jogashima


# How to Install
```
if(!require("devtools"))
	install.packages("devtools")
devtools::install_github("akihirao/ThermoSSM")
```

# Quick Tutorial
https://github.com/akihirao/ThermoSSM/blob/main/quick_tutorial.md


# Future Development Considerations
* Time Series Cross-Validation


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
Supplementary code for estimating sea temperature trends using a linear Gaussian state-space model.
GitHub repository:  
https://github.com/logics-of-blue/sea-temperature-trend-jogashima
