# tempssm

The `tempssm` R package provides tools for state-space analysis of temperature 
time series, focusing on linear Gaussian state-space models estimated via Kalman 
filtering and smoothing as implmented in the `KFAS` package (Helske, 2017).

The package was previously named `ThermoSSM` and has been renamed to `tempssm`.

The core implementation is adapted from the supplementary code provided in
Baba et al. (2024), which is publicly available on GitHub:
https://github.com/logics-of-blue/sea-temperature-trend-jogashima

### Key features

- Designed for temperature time series with **arbitrary seasonal frequency**;
  currently **validated primarily on monthly data**
- Estimates latent states using **linear Gaussian state-space models** combined 
  with Kalman filtering and smoothing
- Models temperature dynamics as a sum of interpretable latent components, 
  including a long-term trend (level), seasonal cucle, autoregressive structure,
  and optional exogenous effects
- Allows users to specify an arbitrary order of the autoregressive component
  (default: AR(1))
- Implements **time-series cross-validation** for model evaluation


# To Install
```r
if(!require("devtools"))
	install.packages("devtools")

devtools::install_github("akihirao/tempssm")
```

# Tutorial
https://github.com/akihirao/tempssm/blob/main/vignettes/tutorial.pdf


# Simple usage

## 1. Load the package and example data

This example uses monthly sea surface temperature (SST) data off Niigata, Japan (February 2002–December 2023).

```r
library(tempssm)
data(niigata_sst)
```

## 2. Fit a state-space model    
The function `tempssm()` fits a state-space model to the monthly temperature
time series using a Gaussian structural model. 

```r
res <- tempssm(hogehoge_ts) # AR(1) is used by default
```
The returned object `res` is an S3 object of class `'tempssm`. You can inspect the results with standard methods.

```r
summary(res)  # summarise results
autoplot(res) # plot results
```


## 3. Using your own data
### 3.1 Quick example
If you already have a ts object, you can pass it directly:
```r
res <- tempssm(my_ts)
```


### 3.2 Preparing your own dataset

#### CSV format 
Prepare a CSV file of monthly temperature time series with the following columns:
- Year
- Month
- Temp

Example:
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

#### Convert CSV to a time series
```{r, eval = FALSE}
my_ts <- tempssm::read_monthly_temp_ts("hogehoge.csv")
```
The function `read_monthly_temp_ts()` internally uses the base R function [`ts()`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/ts) to create a time series object. You may alternatively construct a ts object manually if you need finer control. 

# References
Baba, S., Ishii, H., and Yoshiyama, T. (2024).  
Estimating sea temperature trends using a linear Gaussian state-space model in Jogashima, Kanagawa, Japan.
*Bulletin of the Japanese Society of Fisheries Oceanography*, 88(3), 190–199 (in Japanese with English abstruct)
https://doi.org/10.34423/jsfo.88.3_190

Baba, S. (2024).  
Supplementary code and test data for estimating sea temperature trends using a linear Gaussian state-space model.
GitHub repository:  
https://github.com/logics-of-blue/sea-temperature-trend-jogashima

