# Set Ennvironment

``` r
## Set libraries
library(ThermoSSM)
library(forecast)

library(purrr)   # list_rbind
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)

library(future)
library(parallelly)

has_future <- requireNamespace("future", quietly = TRUE)

if (has_future) {
  if (future::nbrOfWorkers() == 1) {

    # minimqnized setting for vignette
    workers <- min(2, future::availableCores())

    if (future::supportsMulticore()) {
      future::plan(future::multicore, workers = workers)
    } else {
      future::plan(future::multisession, workers = workers)
    }

  }
} else {
  message("future not available; running sequentially")
}
```

# Practice I: Applying a linear Gaussian state-space model to univariate temperature time series

## Objective:

Executing this package with a simple model.

## Loading dataset: Sea surface temperature (SST)

- Data: Monthly sea surface temperature (SST) off Niigata, Japan  
- Unit: Celsius  
- Duration: Feb 2002 to December 2023

This dataset provided by Japan Oceanographic Data Center (JODC),
Hydrographic and Oceanographic Department, Japan Coast Guard. Original
daily SST data were obtained from
<https://www.jodc.go.jp/jodcweb/JDOSS/index.html> and aggregated into
monthly means.

``` r
data(niigata_sst) # load a ts object of SST off Niigata

head(niigata_sst)
```

    ##            Jan       Feb       Mar       Apr       May       Jun
    ## 2002  9.951613  8.332143  9.348387 11.713333 14.529032 18.906667

## Plotting Monthly Time Series of Air Temperature

``` r
plt_niigata_sst <- forecast::autoplot(niigata_sst) +
  labs(y = expression(Temperature~(degree*C)), 
       x = "Time") +
  ggtitle("Monthly SST off Niigata, Japan") +
  theme_classic()

plot(plt_niigata_sst)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# The example data includes missing observation values
```

## Plotting Time Series of Temperature Anomaly

``` r
# Generate temperature anomalies
niigata_sst_anomaly <- ThermoSSM::monthly_anomaly(niigata_sst)

niigata_sst_anomaly_tidy <- tibble(time=time(niigata_sst_anomaly),
                           anomaly = as.numeric(niigata_sst_anomaly)) %>%
  mutate(sign = ifelse(anomaly >= 0, "Positive","Negative"))


plt_niigata_sst_anomaly <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linetype="dotted") +
  geom_line(data = niigata_sst_anomaly_tidy,
            aes(x = time, y = anomaly, fill=""), 
            color="black",alpha = 0.5,
            linewidth = 0.5) +
  guides(fill = "none") +
  labs(
    y = expression(Anomaly~(degree*C)), 
    x = "Year"
  )  +
  theme_classic()


plot(plt_niigata_sst_anomaly)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Plotting Monthly Seasonal Cycle

``` r
monthly_seasonal_cycle_niigata_sst <- ThermoSSM::mean_seasonal_cycle(niigata_sst) 
summary(monthly_seasonal_cycle_niigata_sst)
```

    ##      Month        Temperature    
    ##  Min.   : 1.00   Min.   : 9.365  
    ##  1st Qu.: 3.75   1st Qu.:11.169  
    ##  Median : 6.50   Median :16.318  
    ##  Mean   : 6.50   Mean   :17.036  
    ##  3rd Qu.: 9.25   3rd Qu.:22.294  
    ##  Max.   :12.00   Max.   :26.873

``` r
plt_monthly_seasonal_cycle_niigata_sst <- ggplot(
  data = monthly_seasonal_cycle_niigata_sst,
  aes(x = Month, y = Temperature)
) +
  geom_point(size = 2) +
  geom_line(linetype = "dashed") +
  labs(
    title = "Monthly seasonal cycle of temperature",
    x = "Month",
    y = expression(Temperature~(degree*C))
  ) +
  scale_x_continuous(
    breaks = 1:12,
    labels = sprintf("%02d", 1:12)
  ) +
  theme_classic()


plot(plt_monthly_seasonal_cycle_niigata_sst)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Applying a Linear Gaussian State-Space Model

``` r
# a model with a first-order autoregressive component
res <- lgssm(niigata_sst) # first order of auto-regressive model (ar_order=1: default)
summary(res)
```

    ## ThermoSSM summary
    ## -----------------
    ## Call:
    ## lgssm_seasonal(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
    ##     inits = inits, maxit = maxit, reltol = reltol)
    ## 
    ## Model fit:
    ##   Log-likelihood : -277.95 
    ##   k              : 5 
    ##   AIC            : 565.9 
    ##   Converged      : TRUE 
    ## 
    ## Variance parameters:
    ##   Observation (H): 0.00598564 
    ##   State (Q trend): 1.268118e-07 
    ##   State (Q season): 0.001346139 
    ##   State (Q ar): 0.4097882 
    ## 
    ## Components of auto-regression:
    ##   Order of AR: 1 
    ##   Coefficient of AR1: 0.7442999

``` r
res_ar2 <- lgssm(niigata_sst,ar_order=2) 
summary(res_ar2)
```

    ## ThermoSSM summary
    ## -----------------
    ## Call:
    ## lgssm_seasonal(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
    ##     inits = inits, maxit = maxit, reltol = reltol)
    ## 
    ## Model fit:
    ##   Log-likelihood : -277.95 
    ##   k              : 6 
    ##   AIC            : 567.89 
    ##   Converged      : TRUE 
    ## 
    ## Variance parameters:
    ##   Observation (H): 0.04528328 
    ##   State (Q trend): 1.411893e-07 
    ##   State (Q season): 0.001338972 
    ##   State (Q ar): 0.3463124 
    ## 
    ## Components of auto-regression:
    ##   Order of AR: 2 
    ##   Coefficient of AR1: 0.8278959 
    ##   Coefficient of AR2: -0.0651634

``` r
res_ar3 <- lgssm(niigata_sst,ar_order=3) 
summary(res_ar3)
```

    ## ThermoSSM summary
    ## -----------------
    ## Call:
    ## lgssm_seasonal(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
    ##     inits = inits, maxit = maxit, reltol = reltol)
    ## 
    ## Model fit:
    ##   Log-likelihood : -277.94 
    ##   k              : 7 
    ##   AIC            : 569.89 
    ##   Converged      : TRUE 
    ## 
    ## Variance parameters:
    ##   Observation (H): 0.0002537319 
    ##   State (Q trend): 1.418092e-07 
    ##   State (Q season): 0.001336738 
    ##   State (Q ar): 0.418686 
    ## 
    ## Components of auto-regression:
    ##   Order of AR: 3 
    ##   Coefficient of AR1: 0.7335228 
    ##   Coefficient of AR2: 0.0118387 
    ##   Coefficient of AR3: -0.005368551

## Model selection based on AIC

``` r
# Extract AIC
AIC_ar1 <- extract_AIC(res)
AIC_ar2 <- extract_AIC(res_ar2)
AIC_ar3 <- extract_AIC(res_ar3)

AIC_table_res <- tibble(model=c("AR1","AR2","AR3"),
                        AIC = c(AIC_ar1,AIC_ar2,AIC_ar3)) %>% 
  mutate(delta_AIC = AIC - min(AIC))

AIC_table_res %>% knitr::kable()
```

| model |      AIC | delta_AIC |
|:------|---------:|----------:|
| AR1   | 565.8955 |  0.000000 |
| AR2   | 567.8903 |  1.994836 |
| AR3   | 569.8889 |  3.993417 |

``` r
## AR1 model is better than the other models.
```

## Plotting Level, Drift, Seasonal, and Auto-Regressive Components

``` r
# plot each of components at once
plot(res)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Simple Model Diagnostics

``` r
resid_test_output <- ThermoSSM::wrapper_checkresiduals(res)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
print(resid_test_output)
```

    ## $Ljung_Box
    ## 
    ##  Box-Ljung test
    ## 
    ## data:  std_obs_resid
    ## X-squared = 24.225, df = 24, p-value = 0.4488
    ## 
    ## 
    ## $kurtosis
    ## [1] 3.057601

## Estimated Parameters and Components

``` r
# Smoothing estimates
alpha_hat <- res$kfs$alphahat
head(alpha_hat)
```

    ##            level       slope sea_dummy1 sea_dummy2 sea_dummy3 sea_dummy4
    ## Jan 2002 16.3944 0.005600280  -6.624678  -3.337435  0.6067452  4.7593086
    ## Feb 2002 16.4000 0.005600421  -7.676822  -6.624678 -3.3374348  0.6067452
    ## Mar 2002 16.4056 0.005600415  -7.346316  -7.676822 -6.6246785 -3.3374348
    ## Apr 2002 16.4112 0.005600311  -5.476554  -7.346316 -7.6768220 -6.6246785
    ## May 2002 16.4168 0.005600335  -2.217000  -5.476554 -7.3463157 -7.6768220
    ## Jun 2002 16.4224 0.005600467   2.468021  -2.217000 -5.4765544 -7.3463157
    ##          sea_dummy5 sea_dummy6 sea_dummy7 sea_dummy8 sea_dummy9 sea_dummy10
    ## Jan 2002  8.5280152  9.9007961  6.4159191  2.4680213  -2.217000   -5.476554
    ## Feb 2002  4.7593086  8.5280152  9.9007961  6.4159191   2.468021   -2.217000
    ## Mar 2002  0.6067452  4.7593086  8.5280152  9.9007961   6.415919    2.468021
    ## Apr 2002 -3.3374348  0.6067452  4.7593086  8.5280152   9.900796    6.415919
    ## May 2002 -6.6246785 -3.3374348  0.6067452  4.7593086   8.528015    9.900796
    ## Jun 2002 -7.6768220 -6.6246785 -3.3374348  0.6067452   4.759309    8.528015
    ##          sea_dummy11       arima1
    ## Jan 2002   -7.346316  0.175231707
    ## Feb 2002   -5.476554 -0.377441196
    ## Mar 2002   -2.217000  0.286840259
    ## Apr 2002    2.468021  0.767966399
    ## May 2002    6.415919  0.330186853
    ## Jun 2002    9.900796  0.009042494

``` r
#　Smoothing estimate of level component
level_ts <- ThermoSSM::extract_level_ts(res)

#　Smoothing estimate of drift component
drift_ts <- ThermoSSM::extract_drift_ts(res)

# Average drift rate per year across the full period
mean_drift_year <- mean(drift_ts) 
print(mean_drift_year)
```

    ## [1] 0.05259212

``` r
# Average drift rate per year from Jan 2006 to Dec 2010
ave_drift_2006_2010 <- window(drift_ts,
                              start=c(2006,1),
                               end=c(2010,12)
                               ) %>%  mean()
print(ave_drift_2006_2010)
```

    ## [1] 0.05878758

``` r
# Average drift rate per year from Jan 2011 to Dec 2015
ave_drift_2011_2015 <- window(drift_ts,
                              start=c(2011,1),
                               end=c(2015,12)
                               ) %>%  mean()
print(ave_drift_2011_2015)
```

    ## [1] 0.04891312

``` r
# Average drift rate per year from Jan 2016 to Dec 2020
ave_drift_2016_2020 <- window(drift_ts,
                              start=c(2016,1),
                               end=c(2020,12)
                               ) %>%  mean()
print(ave_drift_2016_2020)
```

    ## [1] 0.04482296

## Plotting Level Component with 95% Confidence Interval

``` r
plt_level_ci <- plot(res,
                     components = "level",
                     ci = TRUE,
                     ci_level = 0.95
                     ) +
  theme_classic()

plot(plt_level_ci)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Plotting Drift Component with 95% Confidence Interval

``` r
plt_drift_ci <- plot(res,
                     components = "drift",
                     ci = TRUE,
                     ci_level = 0.95
                     ) +
  geom_hline(yintercept = 0, lty="dotted") +
  theme_classic()

plot(plt_drift_ci)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# Practice II: Applying a linear Gaussian state-space model to temperature time series with an exogenous variable.

## Objective:

Using the long record of temperature time series data to isolate and
verify the effect of of a large scale climate mode on long-term
temperature variations.

## Loading dataset 1: the Long Record of Air Temperature

- Data: Monthly air temperature at Hohenpeissenberg Meteorological
  Observatory, German  
- Observed Location: 47.801 °N; 11.011 °E  
- Unit: Celsius  
- Duration: Jan 1781 to December 2025

This dataset is the longest record of air temperature observed in
Mountain region across the world. Original dataset provided by Global
Historical Climatology Network monthly (GHCNm) ver.4 was obtained from
<https://www.ncei.noaa.gov/data/global-historical-climatology-network-monthly/v4>

``` r
data(hmo_temp) # loading a ts object of air temperature at the HMO station

head(hmo_temp)
```

    ##        Jan   Feb   Mar   Apr   May   Jun
    ## 1781 -1.76 -1.03  2.37  8.71 12.25 14.54

## Loading dataset 2: NAO index as an exogenous valiable

- Data: Monthly North Atlantic Oscillation (NAO) index (Hurrell)  
- Duration: Jan 1865 to Jun 2023

The winter (December thru March) station-based index of the NAO is based
on the difference of normalized sea level pressure (SLP) between Lisbon,
Portugal and Stykkisholmur/Reykjavik, Iceland since 1864. Positive
values of the NAO index are typically associated with
stronger-than-average westerlies over the middle latitudes, more intense
weather systems over the North Atlantic and wetter/milder weather over
western Europe. Monthly, seasonal and annual indices using slightly
different data sources for the southern station are also available.

Source dataset provided by the Climate Analysis Section, NCAR, Boulder,
USA, Hurrell (2003) was obtained from
<https://www.jodc.go.jp/jodcweb/JDOSS/index.html>.

``` r
data(nao) # load a ts object of NAO index

head(nao)
```

    ##       Jan  Feb  Mar  Apr  May  Jun
    ## 1865 -0.6 -1.2  0.2 -0.2 -0.4  0.0

## Intersecting the two ts objects of temperature and NAO

``` r
# Generate an object on a shared timeline
hmo_nao_ts <- ts.intersect(hmo_temp, nao)
colnames(hmo_nao_ts) <- c("Temp", "NAO")

start(hmo_nao_ts)
```

    ## [1] 1865    1

``` r
end(hmo_nao_ts)
```

    ## [1] 2023    6

``` r
# Generate an ts object of HMO air-temperature with the shared timeline
hmo_temp_common <- hmo_nao_ts[,"Temp"]
start(hmo_temp_common)
```

    ## [1] 1865    1

``` r
end(hmo_temp_common)
```

    ## [1] 2023    6

``` r
# Generate an ts object of NAO with the shared timeline
nao_common <- hmo_nao_ts[,"NAO"]

# Execute label_ts_mono() function to label exogenous variable(s)!!
nao_common <- label_ts_mono(nao_common, label="NAO")
start(nao_common)
```

    ## [1] 1865    1

``` r
end(nao_common)
```

    ## [1] 2023    6

## Plotting Time Series of air temperature and NAO index

``` r
plt_homo_temp <- forecast::autoplot(hmo_temp_common) +
  labs(x = "Time", y = expression(Temperature~(degree*C))) +
  ggtitle("Air temperature at Hohenpeissenberg Meteorological Observatory") +
  theme_classic()

plt_nao <- forecast::autoplot(nao_common) +
  labs(x = "Time", y = "NAO index") +
  ggtitle("NAO index") +
  theme_classic()

cowplot::plot_grid(plt_homo_temp,
                   plt_nao,
                   ncol=1)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Applying a Linear Gaussian State-Space Model without exogenous variables

``` r
res_without <- lgssm(hmo_temp_common) 
summary(res_without)
```

    ## ThermoSSM summary
    ## -----------------
    ## Call:
    ## lgssm_seasonal(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
    ##     inits = inits, maxit = maxit, reltol = reltol)
    ## 
    ## Model fit:
    ##   Log-likelihood : -4116.01 
    ##   k              : 5 
    ##   AIC            : 8242.02 
    ##   Converged      : TRUE 
    ## 
    ## Variance parameters:
    ##   Observation (H): 2.198694 
    ##   State (Q trend): 1.115939e-08 
    ##   State (Q season): 0.0001822552 
    ##   State (Q ar): 2.085541 
    ## 
    ## Components of auto-regression:
    ##   Order of AR: 1 
    ##   Coefficient of AR1: 0.2284355

``` r
AIC_without <- extract_AIC(res_without)
plot(res_without)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Applying a Linear Gaussian State-Space Model with exogenous variables

``` r
res_with <- lgssm(temp_data = hmo_temp_common,exo_data = nao_common) 
summary(res_with)
```

    ## ThermoSSM summary
    ## -----------------
    ## Call:
    ## lgssm_seasonal(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
    ##     inits = inits, maxit = maxit, reltol = reltol)
    ## 
    ## Model fit:
    ##   Log-likelihood : -4063.67 
    ##   k              : 6 
    ##   AIC            : 8139.34 
    ##   Converged      : TRUE 
    ## 
    ## Variance parameters:
    ##   Observation (H): 2.623818 
    ##   State (Q trend): 9.64448e-09 
    ##   State (Q season): 0.0002605165 
    ##   State (Q ar): 1.402478 
    ## 
    ## Components of auto-regression:
    ##   Order of AR: 1 
    ##   Coefficient of AR1: 0.2676844 
    ## Exogenous variable    NAO 
    ## Estimated coefficient     0.2911478 
    ## Lower CI  0.2376886 
    ## Upper CI  0.3446071

``` r
AIC_with <- extract_AIC(res_with)
plot(res_with)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Model comparison based on AICs

``` r
models_AICs <- tibble(
  model = c("Without","With"),
  AIC = c(AIC_without,AIC_with),
  delta_AIC = min(AIC)-AIC
  )

models_AICs %>% knitr::kable()
```

| model   |      AIC | delta_AIC |
|:--------|---------:|----------:|
| Without | 8242.017 | -102.6742 |
| With    | 8139.343 |    0.0000 |

``` r
## Model with the exogenous variable of NAO is remarkably better than model without one. 
```

## time series Cross-Validation

### In this tutorial, the number of tsCV iterations has been set to eigth to reduce execution time. Please ensure you perform a sufficient number of iterations during actual validation.

``` r
# ts cross-validation of the model without exogenous variables

## Generate a list of training and test datasets with their indices
# Procedure for constructing year-based time-series cross-validation folds:
# 
# First training data: January 1865–December 1950;
# First test data: one year starting from January 1951
# 
# Second training data: January 1875–December 1960;
# Second test data: one year starting from January 1961
# ...
# 
# Eighth training data: January 1875–December 2020;
# 15th test data: one year starting from January 2021
#
# These folds are automatically generated by the ts_cv_folds() function.
# Generate training and test dataset
folds_without <- ts_cv_folds(
  temp_data = hmo_temp_common,
  exo_data = NULL,
  initial = 1032, # 228 monthly observations from Jan 1865 to Dec 1950
  horizon = 12, # forecast 12 monthly time-series
  step = 120, # training data is moved by 120 months (10 years) steps
  fixed_window = FALSE,
  allow_partial = FALSE
  )

folds_with <- ts_cv_folds(
  temp_data = hmo_temp_common,
  exo_data = nao_common,
  initial = 1032, # 228 monthly observations from Jan 1865 to Dec 1950
  horizon = 12, # forecast 12 monthly time-series
  step = 120, # training data is moved by 12 months (10 years) steps
  fixed_window = FALSE,
  allow_partial = FALSE
)

## Performe tsCV for the model without NAO index
cv_meta_without <- rolling_origin_tsCV(folds_without,
                                       fold_ids=seq(1:8), 
                                       ar_order=1) 
# Hold on a few minutes 

cv_without <- map(cv_meta_without, pluck, "CV", .default = NULL)

cv_without_tidy <- cv_without %>%
  map(~ as.data.frame(.x) %>%tibble::rownames_to_column(var = "set")) %>%
  purrr::list_rbind(names_to = "cv_id") %>%
  mutate(Model="Without") %>%
  dplyr::relocate(Model,cv_id, set)  # 識別列を先頭に

print(cv_without_tidy)
```

    ##     Model cv_id      set         ME     RMSE      MAE       MPE      MAPE
    ## 1 Without fold1 Test set  0.4889957 1.474649 1.148685 58.787802 103.61067
    ## 2 Without fold2 Test set  1.1904456 2.606346 2.009215 24.273406  32.54334
    ## 3 Without fold3 Test set  0.1135489 2.037836 1.857256 96.961995 117.51916
    ## 4 Without fold4 Test set -0.1159486 1.575582 1.099652 23.311897  25.49696
    ## 5 Without fold5 Test set -0.3779929 1.945921 1.665165  7.593377  33.60758
    ## 6 Without fold6 Test set  0.1382091 2.651390 2.270733       Inf       Inf
    ## 7 Without fold7 Test set  1.1811413 2.019612 1.629615 24.218516  34.10742
    ## 8 Without fold8 Test set -0.8494083 1.993760 1.683971  2.033586  34.11603
    ##           ACF1 Theil's U MASE_naive MASE_snaive
    ## 1  0.202919625 0.5186926  0.3277553   0.4858329
    ## 2  0.052893544 0.8037302  0.5754156   0.8526872
    ## 3 -0.315093475 0.5675867  0.5308151   0.7917307
    ## 4 -0.013508238 0.4340856  0.3169881   0.4734971
    ## 5 -0.007453192 0.5675805  0.4791300   0.7219113
    ## 6 -0.435854697       NaN  0.6541575   0.9848380
    ## 7  0.099740227 0.5470913  0.4688726   0.7049821
    ## 8 -0.277090168 0.5301606  0.4855664   0.7242332

``` r
## Performe tsCV for the model with NAO index
cv_meta_with <- rolling_origin_tsCV(folds_with,
                                    fold_ids=seq(1:8),
                                    ar_order=1)
# Hold on a few minutes 

cv_with <- map(cv_meta_with, pluck, "CV", .default = NULL)


cv_with_tidy <- cv_with %>%
  map(~ as.data.frame(.x) %>%tibble::rownames_to_column(var = "set")) %>%
  purrr::list_rbind(names_to = "cv_id") %>%
  mutate(Model="With") %>%
  dplyr::relocate(Model,cv_id, set)  # Reorder columns to place identifier variables first

print(cv_with_tidy)
```

    ##   Model cv_id      set          ME     RMSE      MAE        MPE      MAPE
    ## 1  With fold1 Test set  0.48543195 1.295995 1.026690  52.583598  83.53820
    ## 2  With fold2 Test set  1.11609061 2.499648 1.992019   9.128501  37.62516
    ## 3  With fold3 Test set  0.16374746 2.079393 1.845031 107.098558 129.75021
    ## 4  With fold4 Test set -0.04621108 1.738008 1.270618  24.235775  27.14959
    ## 5  With fold5 Test set -0.43349108 2.040727 1.740272   7.616265  35.40459
    ## 6  With fold6 Test set  0.32903989 2.527949 2.168827        Inf       Inf
    ## 7  With fold7 Test set  0.81409922 1.752448 1.368590   9.449692  23.09277
    ## 8  With fold8 Test set -0.61119759 1.738435 1.471531   3.981070  31.65661
    ##          ACF1 Theil's U MASE_naive MASE_snaive
    ## 1 -0.01964670 0.3750799  0.2929461   0.4342352
    ## 2  0.10503306 0.7799049  0.5704907   0.8453892
    ## 3 -0.32985210 0.6535652  0.5273213   0.7865196
    ## 4 -0.01644048 0.4551757  0.3662711   0.5471130
    ## 5 -0.01897327 0.6042778  0.5007411   0.7544730
    ## 6 -0.48586337       NaN  0.6248004   0.9406408
    ## 7 -0.11047192 0.4093789  0.3937705   0.5920609
    ## 8 -0.33494646 0.5134535  0.4243103   0.6328684

``` r
cv_tidy <- bind_rows(cv_without_tidy,cv_with_tidy)
cv_tidy$Model <- factor(cv_tidy$Model,
                        levels=c("Without","With"))

plot_MAE <- ggplot(data=cv_tidy,
                   aes(x=Model,y=MAE)) +
  geom_boxplot()

plot_MAPE <- ggplot(data=cv_tidy,
                   aes(x=Model,y=MAPE)) +
  geom_boxplot()

plot_MASE_naive <- ggplot(data=cv_tidy,
                   aes(x=Model,y=MASE_naive)) +
  geom_boxplot()

plot_MASE_snaive <- ggplot(data=cv_tidy,
                   aes(x=Model,y=MASE_snaive)) +
  geom_boxplot()

cowplot::plot_grid(plot_MAE,plot_MAPE,plot_MASE_naive,plot_MASE_snaive,nrow=1)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

# Appendix: Utility Function

## Utility function of sst_jma2ts()

### This function downloads publicly available daily mean sea-surface temperature (SST) data

### for Japanese coastal waters provided by the Japan Meteorological Agency (JMA),

### and returns the corresponding monthly mean time series as an object of class ts.

``` r
#' @param sea_area_id
#' Numeric sea area ID. The default is NULL.
#' For example, 138 corresponds to the coastal waters off southern Ibaraki, Japan.
#' A list of sea area IDs and their corresponding regions is available at:
#' \url{https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html}
```

``` r
# Download SST data for the southern coastal waters of Ibaraki Prefecture, Japan
sst_138_ts <- ThermoSSM::sst_jma2ts(sea_area_id = 138) 
head(sst_138_ts)
```

    ##           Jan      Feb      Mar      Apr      May      Jun
    ## 1982 15.04419 14.22500 13.63903 15.31933 17.52258 19.52300
