# tempssm

R package **tempssm** provides tools for fitting and analyzing 
linear Gaussian state-space models to temperature time series, 
with a focus on climate and environmental applications.
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


# To Install
```r
if(!require("devtools"))
	install.packages("devtools")
devtools::install_github("akihirao/tempssm")

# install with vignettes
devtools::install_github("akihirao/tempssm", build_vignettes=TRUE)
```

# Quick Tutorial
```r
browseVignettes("tempssm")
```
https://github.com/akihirao/tempssm/blob/main/vignettes/quick_tutorial.md


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
library(tempssm)

hogehoge_ts <- monthly_temp_csv2ts("hogehoge.csv")
```
The function `monthly_temp_csv2ts()` internally uses the base R function [`ts()`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ts) to create a time series object. Users may alternatively convert their data manually using `ts()` if finer control over the time series structure is required.

### 3. Execute the state-space modelling    
The function `lgssm()` fits a state-space model to the monthly temperature
time series using a Gaussian structural model.

```r
res <- lgssm(hogehoge_ts)

summary(res)  # summarise results
plot(res) # plot results
```
The returned object `res` is an S3 object of class `'tempssm`, which contains the fitted state-space model results. The object `res` can be analyzed using methods such as `summary()` and `plot()`.


# Example Data

The package includes four example datasets for demonstration and testing of state-space temperature analyses. Three out of the four datasets are provided in `.rda` format.
```r
data(niigata_sst)  # ts object of monthly sea surface temperature off Niigata, Japan
data(fuji_temp)    # ts object of monthly air temperature at the summit of Mt. Fuji, Japan
data(hmo_temp)     # ts object of monthly air temperature at the Hohenpeissenberg Meteorological Observatory (HMO), Germany
```
The fourth dataset is provided in .csv format and included in inst/extdata. It contains a monthly temperature time series measured at Mount Akadake, Hokkaido, Japan (elevation 1,840 m; 43.6766°N, 142.9423°E). The data originate from the Monitoring Sites 1000 Project conducted by the Ministry of the Environment of Japan (KOZ01.zip, downladed from https://www.biodic.go.jp/moni1000/findings/data/index.html).
```r
path <- system.file("extdata", "example_monthly_temp.csv", package = "tempssm")
akadake_temp_info <- readr::read_csv(path)
head(akadake_temp_info)

# conver from data frame to monthly ts object
akadake_temp <- monthly_temp_csv2ts(akadake_temp_info)
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

