# ThermoSSM

**ThermoSSM** is an R package for state-space analysis of temperature time series.
It implements linear Gaussian state-space models and performs inference using
Kalman filtering and smoothing.

The core implementation is adapted from the supplementary code provided in
Baba et al. (2024), which is publicly available on GitHub:
https://github.com/logics-of-blue/sea-temperature-trend-jogashima

### Key features

- Designed for temperature time series data with **arbitrary seasonal frequency**;
  currently **validated mainly on monthly data**
- Applies **linear Gaussian state-space models** to estimate latent states
  using Kalman filtering and smoothing
- Represents temperature dynamics as a sum of interpretable latent components:
  long-term trend, seasonal cycle, autoregressive structure,
  and optional exogenous effects
- Implements **time-series cross-validation** for model evaluation


# How to Install
```r
if(!require("devtools"))
	install.packages("devtools")
devtools::install_github("akihirao/ThermoSSM")

# install with vignettes
devtools::install_github("akihirao/ThermoSSM", build_vignettes=TRUE)
```

# Simple usage

### 1. Prepare a CSV file  

For monthly temperature data, prepare a CSV file with a header row.
By default, the column names should be `Year`, `Month`, and `Temp`.

```text
Year,Month,Temp
2010,8,13.6
2010,9,6.8
2010,10,NA
2010,11,-1.4
...
```
* Use NA for missing temperature values, and always keep the corresponding Year and Month entries.
* The CSV file must be comma-separated and UTF-8 encoded.

### 2. Load the CSV file into R and covert it to a time series oject     

If the CSV file is named `hogehoge.csv`, load the data and convert it to a time series object using the following commands:

```r
library(ThermoSSM)

hogehoge_ts <- ThermoSSM::monthly_temp_csv2ts("hogehoge.csv")
```

### 3. Execute the state-space modelling    
The function `lgsmm()` fits a state-space model to the monthly temperature
time series using a Gaussian structural model.

```r
res <- lgsmm(hogehoge_ts)

summary(res)  # summarise results
plot(res) # plot results
```
The returned object `res` is an S3 object of class `'ThermoSSM`, which contains the fitted state-space model results. The object `res` can be analyzed using methods such as `summary()` and `plot()`.
``

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

