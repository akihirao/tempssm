R package **tempssm** provides tools for state-space analysis of
temperature time series. It provides tools for fitting linear Gaussian
state-space models and conducting inference via Kalman filtering and
smoothing, implemented using the KFAS R package (Helske, 2017).

## Key features

- Designed for temperature time series with arbitrary seasonal
  frequencies; currently validated primarily on monthly data
- Estimates latent states using linear Gaussian state-space models
  combined with Kalman filtering and smoothing
- Models temperature dynamics as a sum of interpretable latent
  components, including a long-term trend, seasonal cycle,
  autoregressive structure, and optional exogenous effects
- Allows users to specify an arbitrary order of the autoregressive
  component (default: AR(1))
- Implements time-series cross-validation for model evaluation

## Input Data Format

Input data for **tempssm** must be supplied as an R `ts` object, which
represents a regularly spaced time series (see `?stats::ts` or
<https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ts.html>).

To support data preparation, the package includes utility functions that
convert external observational data into `ts` objects suitable for model
fitting (see Appendix).

## To Install

``` r
if(!require("devtools"))
    install.packages("devtools")
devtools::install_github("akihirao/tempssm")
```

## Set Ennvironment for Practices

``` r
## Set libraries
library(tempssm)
library(forecast)

library(purrr) 
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

## Practice I: Applying a Linear Gaussian State-Space Model

## to a Univariate Temperature Time Series

### Objective

The objective of this practice is to demonstrate the basic application
of a linear Gaussian state-space model to a univariate temperature time
series. This example serves as an introduction to the modeling framework
and highlights the role of autoregressive dynamics without the inclusion
of exogenous variables.

### Loading the Dataset

A sample sea surface temperature (SST) dataset is included in the
package.

- **Dataset**: Monthly sea surface temperature (SST) off Niigata,
  Japan  
- **Unit**: Degrees Celsius  
- **Period**: February 2002 to December 2023  
  \`\`

This dataset is derived from observations archived at the Japan
Oceanographic Data Center (JODC), Hydrographic and Oceanographic
Department, Japan Coast Guard. The original daily SST data were obtained
from  
<https://www.jodc.go.jp/jodcweb/JDOSS/index.html>  
and subsequently aggregated into monthly means.

``` r
data(niigata_sst) # load a ts object of SST off Niigata
head(niigata_sst)
```

    ##            Jan       Feb       Mar       Apr       May       Jun
    ## 2002  9.951613  8.332143  9.348387 11.713333 14.529032 18.906667

``` r
summary(niigata_sst)
```

    ##       Temp       
    ##  Min.   : 7.707  
    ##  1st Qu.:11.217  
    ##  Median :16.345  
    ##  Mean   :17.033  
    ##  3rd Qu.:22.787  
    ##  Max.   :28.897  
    ##  NA's   :2

The dataset includes two missing observations, which are retained and
handled explicitly within the state-space modeling framework.

### Plotting the Monthly SST Time Series

We begin by visualizing the monthly SST time series to examine its
overall structure, including apparent trends, seasonal variability, and
missing observations.

``` r
plt_niigata_sst <- forecast::autoplot(niigata_sst) +
  labs(y = expression(Temperature~(degree*C)), 
       x = "Time") +
  ggtitle("Monthly SST off Niigata, Japan") +
  theme_classic()

plot(plt_niigata_sst)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Applying a Linear Gaussian State-Space Model

Here, we apply linear Gaussian state-space models to the univariate SST
time series and examine the effect of the autoregressive (AR) order on
model behavior.  
Specifically, three models are fitted with the order of the
autoregressive component varying from 1 to 3, while all other model
components, including the explicit seasonal cycle, are kept the same.

By comparing models with different AR orders, we assess how short- and
longer-term temporal dependencies are represented within the state-space
framework.  
The model summaries provide parameter estimates and diagnostics that can
be used to evaluate the adequacy of different autoregressive orders.
Model comparison and selection are discussed in the following section.

``` r
# model with first-order autoregressive component
res <- tempssm(niigata_sst) # first order of auto-regressive model (ar_order=1: default)
summary(res)
```

    ## tempssm summary
    ## -----------------
    ## Call:
    ## tempssm_season(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
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
# model with second-order autoregressive component
res_ar2 <- tempssm(niigata_sst,ar_order=2) 
summary(res_ar2)
```

    ## tempssm summary
    ## -----------------
    ## Call:
    ## tempssm_season(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
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
# model with third-order autoregressive component
res_ar3 <- tempssm(niigata_sst,ar_order=3) 
summary(res_ar3)
```

    ## tempssm summary
    ## -----------------
    ## Call:
    ## tempssm_season(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
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

### Model selection based on AIC

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

AR1 model is better than the other models.

### Plotting Level, Drift, Seasonal, and Auto-Regressive Components

``` r
# plot each of components at once
plot(res)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Simple Model Diagnostics

``` r
resid_test_output <- tempssm::wrapper_checkresiduals(res)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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

Autocorrelation of residuals was not significant by Ljung-Box test (P \>
0.05).

### Estimated Parameters and Components

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
level_ts <- tempssm::extract_level_ts(res)

#　Smoothing estimate of drift component
drift_ts <- tempssm::extract_drift_ts(res)

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

### Plotting Level and Drift Components with 95% Confidence Interval

We visualize the estimated long-term evolution of temperature levels and
their rates of change (drift) by extracting the corresponding latent
components from the state-space model. Seasonal variability and
autoregressive dependence are separated out, allowing the underlying
trend behavior to be examined more clearly.

``` r
plt_level_ci <- plot(res,
                     components = "level",
                     ci = TRUE,
                     ci_level = 0.95
                     ) +
  theme_classic()

plt_drift_ci <- plot(res,
                     components = "drift",
                     ci = TRUE,
                     ci_level = 0.95
                     ) +
  geom_hline(yintercept = 0, lty="dotted") +
  theme_classic()


plt_level_drift_ci <- plt_level_ci + plt_drift_ci + patchwork::plot_layout(ncol=1)

plot(plt_level_drift_ci)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

The level component shows a persistent upward trend in sea surface
temperature over the study period. The drift component indicates a
relatively stable positive rate of change, with an estimated average
annual increase of approximately 0.0525 °C. The shaded gray areas
represent 95% confidence intervals for the estimated latent states,
illustrating the uncertainty associated with the inferred long-term
trend and its rate of change.

## Practice II: Applying a Linear Gaussian State-Space Model

## to a Temperature Time Series with an Exogenous Variable

### Objective

The objective of this practice is to use a long record of temperature
time series data to investigate and quantify the influence of a
large-scale climate mode on long-term temperature variations.
Specifically, we examine the effect of the North Atlantic Oscillation
(NAO) as an exogenous variable within a state-space modeling framework.

### Loading Dataset 1: Long-Term Air Temperature Record

- **Data**: Monthly air temperature at the Hohenpeissenberg
  Meteorological Observatory, Germany  
- **Location**: 47.801°N, 11.011°E  
- **Unit**: Degrees Celsius  
- **Period**: January 1781 to December 2025

This dataset represents one of the longest continuous instrumental
records of air temperature available worldwide, observed at a mountain
meteorological station. The original data were obtained from the Global
Historical Climatology Network Monthly (GHCNm) version 4, distributed by
the National Centers for Environmental Information (NCEI), NOAA:
<https://www.ncei.noaa.gov/data/global-historical-climatology-network-monthly/v4>.

``` r
data(hmo_temp) # loading a ts object of air temperature at the HMO station
head(hmo_temp)
```

    ##        Jan   Feb   Mar   Apr   May   Jun
    ## 1781 -1.76 -1.03  2.37  8.71 12.25 14.54

### Loading Dataset 2: NAO Index as an Exogenous Variable

- **Data**: Monthly North Atlantic Oscillation (NAO) index (Hurrell)  
- **Period**: January 1865 to June 202

Monthly NAO indices derived using related data　sources are available
from the Climate Analysis Section, National Center for Atmospheric
Research (NCAR), Boulder, USA (Hurrell, 2003). The dataset used in this
package was obtained from:
<https://www.ncei.noaa.gov/access/monitoring/nao/>.

``` r
data(nao) # load a ts object of NAO index
head(nao)
```

    ##       Jan  Feb  Mar  Apr  May  Jun
    ## 1865 -0.6 -1.2  0.2 -0.2 -0.4  0.0

### Intersecting the Temperature and NAO Time Series

For state-space modeling with exogenous variables, all input time series
must share a common and aligned time index. In this step, the
temperature and NAO time series are restricted to their overlapping
period by trimming the leading and trailing portions, ensuring that both
datasets cover an identical time span.

The function `ts.intersect()` is used to align the two `ts` objects on a
shared timeline, returning a multivariate time series containing only
the common period.

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

### Plotting Time Series of Air Temperature and NAO Index

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

![](quick_tutorial_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### Applying a Model Without an Exogenous Variable

We first fit a baseline state-space model that does not include any
exogenous variables. This model serves as a reference case in which
temperature variability is explained solely by the latent trend,
seasonal cycle, and autoregressive dependence.

``` r
res_without <- tempssm(hmo_temp_common) 
summary(res_without)
```

    ## tempssm summary
    ## -----------------
    ## Call:
    ## tempssm_season(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
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

The fitted model converges successfully and includes a first-order
autoregressive \[AR(1)\] component. Prior testing of autoregressive
orders from AR(1) to AR(3) indicated that AR(1) provided the best model
fit for both the baseline model and the model including the exogenous
PDO index. For brevity, the detailed results of this comparison are
omitted here; interested users are encouraged to explore alternative AR
orders within their own analyses.

This baseline model provides a useful benchmark for evaluating the
additional explanatory power of the PDO index introduced in the
following section.

### Applying Model With an Exogenous Variable

``` r
res_with <- tempssm(temp_data = hmo_temp_common,
                    exo_data = nao_common) 
summary(res_with)
```

    ## tempssm summary
    ## -----------------
    ## Call:
    ## tempssm_season(temp_data = temp_data, exo_data = exo_data, ar_order = ar_order, 
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

The estimated coefficient for the exogenous NAO index is positive
(0.29), and its 95% confidence interval does not include zero
``` math
0.24,
0.34
```
, indicating a statistically significant relationship between NAO
variability and local air temperature at the Hohenpeissenberg station.

Specifically, the model suggests that a one-unit increase in the NAO
index is associated with an average increase of approximately 0.29 °C in
monthly air temperature, after accounting for the underlying trend,
seasonal cycle, and autoregressive dependence. This result implies that
positive phases of the NAO contribute systematically to warmer
temperature conditions at the study site.

Importantly, this effect is identified in addition to the internal
dynamics of the temperature time series, rather than being an artifact
of an improved representation of the trend or temporal dependence.
Compared with the baseline model without exogenous variables, the
statistical significance of the NAO coefficient indicates that the NAO
acts as an independent large-scale climate driver influencing local
temperature variability.

We next compare the overall goodness of fit between the two models using
the Akaike Information Criterion (AIC).

### Model comparison based on AICs

Model selection criteria such as the AIC provide a quantitative measure
of model adequacy that balances goodness of fit against model
complexity. Lower AIC values indicate a more parsimonious model with
better support from the data.

``` r
AIC_without <- extract_AIC(res_without)
AIC_with <- extract_AIC(res_with)

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

The model including the NAO index as an exogenous variable exhibits a
substantially lower AIC than the baseline model without exogenous
variables, indicating a markedly better overall fit. This result
supports the conclusion that explicitly accounting for NAO variability
improves the statistical description of the temperature time series
beyond what is achieved by internal dynamics alone.

It should be noted that AIC-based comparisons are meaningful only among
models fitted to the same time series over an identical period.
Moreover, in state-space models, AIC differences reflect both regression
effects and changes in the stochastic structure of latent components.
Therefore, AIC should be used as a complementary diagnostic alongside
parameter estimates and their uncertainty, rather than as a sole basis
for inference.

To further assess the robustness and predictive performance of the
models, we employ time-series cross-validation as described bellow.

### Plotting Level and Drift Components with 95% CI

``` r
plt_level_without_ts <- plot(res_without, 
                             components=c("level"),
                             ci=TRUE) +
  labs(title="Model without the NAO index") +
  theme_classic()

plt_drift_without_ts <- plot(res_without, 
                             components=c("drift"),
                             ci=TRUE) + 
  labs(title="Model without the NAO index") +
  theme_classic()

plt_level_drift_without_ts <- plt_level_without_ts + 
  plt_drift_without_ts + 
  patchwork::plot_layout(ncol=1)

plot(plt_level_drift_without_ts)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plt_level_with_ts <- plot(res_with,
                          components=c("level"),
                          ci=TRUE)+ 
  labs(title="Model with the NAO index") +
  theme_classic()

plt_drift_with_ts <- plot(res_with,
                          components=c("drift"),
                          ci=TRUE) + 
  labs(title="Model with the NAO index") + 
  theme_classic()

plt_level_drift_with_ts <- plt_level_with_ts + 
  plt_drift_with_ts + 
  patchwork::plot_layout(ncol=1)

plot(plt_level_drift_with_ts)
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

Gray areas in the above graph show 95% confidence interval.

### Time Series Cross-Validation (tsCV)

Time series cross-validation (tsCV) approach evaluates model performance
based on out-of-sample prediction errors while respecting the temporal
ordering of the data, and thus provides an additional, independent
perspective on model adequacy.

Time-series cross-validation is conducted by repeatedly fitting the
model to an expanding training window and evaluating its predictive
performance on subsequent observations. Unlike random cross-validation,
this procedure avoids information leakage from the future to the past
and is therefore well suited for temporal data.

By comparing cross-validation metrics for models with and without the
exogenous PDO variable, we assess whether the improvement suggested by
AIC is also reflected in out-of-sample predictive skill. \`\`

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
# Eighth test data: one year starting from January 2021
#
# These folds are automatically generated by the ts_cv_folds() function.
# Generate training and test dataset

folds_without <- ts_cv_folds(
  temp_data = hmo_temp_common,
  exo_data = NULL,
  initial = 1032, # 1032 monthly observations from Jan 1865 to Dec 1950
  horizon = 12, # forecast 12 monthly time-series
  step = 120, # training data is moved by 120 months (10 years) steps
  fixed_window = FALSE,
  allow_partial = FALSE
  )

folds_with <- ts_cv_folds(
  temp_data = hmo_temp_common,
  exo_data = nao_common,
  initial = 1032, # 1032 monthly observations from Jan 1865 to Dec 1950
  horizon = 12, # forecast 12 monthly time-series
  step = 120, # training data is moved by 12 months (10 years) steps
  fixed_window = FALSE,
  allow_partial = FALSE
)

## Performe tsCV for the model without the exogenous variable of NAO
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
## Performe tsCV for the model with the exogenous variable of NAO
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

![](quick_tutorial_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Analysis of tsCV shows that the model with the exogenous variable is
better than the model without one. In this tutorial, the number of tsCV
iterations has been set to eigth to reduce execution time. Please ensure
you perform a sufficient number of iterations during actual validation.

Taken together, the results consistently support the inclusion of the
NAO index as an exogenous variable in the state-space model. The NAO
coefficient is statistically significant, indicating a clear local
effect on temperature variability, while the improvement in AIC
demonstrates enhanced overall model fit. Furthermore, time-series
cross-validation confirms that this improvement translates into superior
out-of-sample predictive performance.

These complementary lines of evidence provide a robust and multifaceted
justification for adopting the exogenous-variable model.

## Appendix: Utility Functions

### Utility function: `monthly_temp_csv2ts()`

The function `monthly_temp_csv2ts()` converts monthly temperature data
stored in a CSV file into an R `ts` object. It is intended to facilitate
the ingestion of externally prepared time-series data into tempssm by
enforcing a simple and consistent data format.

#### Example

For monthly temperature data, prepare a CSV file with a header row. By
default, the column names must be `Year`, `Month`, and `Temp`, where
`Temp` represents the observed temperature value.

``` text
Year,Month,Temp
2010,8,13.6
2010,9,6.8
2010,10,NA
2010,11,-1.4
...
```

- Use NA for missing temperature values, and always keep the
  corresponding Year and Month entries.
- The CSV file must be comma-separated and UTF-8 encoded.

The following example uses a sample CSV file included in the package.

``` r
path <- system.file("extdata", "example_monthly_temp.csv", package = "tempssm")
example_temp <- monthly_temp_csv2ts(path)
head(example_temp)
```

    ##        Jan Feb Mar Apr May Jun Jul   Aug   Sep   Oct   Nov   Dec
    ## 2010                                13.6   6.8   0.2  -6.8 -12.5
    ## 2011 -18.8

### Utility function: `monthly_temp_df2ts()`

The function `monthly_temp_df2ts()` converts a data frame containing
monthly temperature data into an R `ts` object. It is designed to
support workflows where temperature data are first imported or prepared
as a data frame before being used for time-series analysis.

The input data frame is expected to contain at least two columns: a date
column (`Date`) and a temperature column (`Temp`).

#### Example

``` r
# Load example CSV file included in the package
path <- system.file("extdata", "example_monthly_temp.csv", package = "tempssm")

# Create a data frame with Date and Temp columns
example_temp_df <- readr::read_csv(path) %>%
  dplyr::mutate(Date = as.Date(paste0(Year, "-", Month, "-01"))) %>%
  dplyr::select(Date, Temp)

# Convert the data frame to a ts object
example_temp_ts <- monthly_temp_df2ts(example_temp_df)

head(example_temp_ts)
```

    ##        Jan Feb Mar Apr May Jun Jul   Aug   Sep   Oct   Nov   Dec
    ## 2010                                13.6   6.8   0.2  -6.8 -12.5
    ## 2011 -18.8

### Utility function: `sst_jma2ts()`

The function `sst_jma2ts()` downloads publicly available daily mean
sea-surface temperature (SST) data for Japanese coastal waters provided
by the Japan Meteorological Agency (JMA). It aggregates the daily values
into monthly means and returns the resulting time series as an object of
class `ts`.

#### Example

The following example demonstrates how to download SST data for the
southern coastal waters of Ibaraki Prefecture, Japan.

The argument `sea_area_id` specifies the numeric identifier of a sea
area defined by JMA. For example, `138` corresponds to the coastal
waters off southern Ibaraki. A list of available sea area IDs and their
corresponding regions is provided by JMA at:

<https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html>

``` r
sst_138_ts <- tempssm::sst_jma2ts(sea_area_id = 138)
head(sst_138_ts)
```

    ##           Jan      Feb      Mar      Apr      May      Jun
    ## 1982 15.04419 14.22500 13.63903 15.31933 17.52258 19.52300

### Utility function: `monthly_anomaly()`

The function `monthly_anomaly()` computes monthly anomalies from a time
series provided as an R `ts` object. Monthly anomalies are calculated by
subtracting the corresponding long-term monthly mean from each
observation, thereby removing the seasonal cycle while preserving
interannual variability.

The reference period used to compute the monthly climatology can be
specified via a function argument; by default, the climatology is
computed over the entire available time series (see `?monthly_anomaly`).

This transformation is useful for exploratory analysis and modeling
applications that focus on departures from typical seasonal conditions.

#### Example

``` r
# Generate temperature anomalies
data(niigata_sst)
niigata_sst_anomaly <- monthly_anomaly(niigata_sst)
plot(niigata_sst_anomaly, ylab=expression(Anomaly~(degree*C))) 
```

![](quick_tutorial_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### Utility function: `mean_seasonal_cycle()`

The function `mean_seasonal_cycle()` computes the climatological mean
seasonal cycle from a `ts` object by averaging each seasonal value
across years. It is primarily intended for exploratory analysis and
visualization of the seasonal structure in temperature time series.

#### Example

``` r
data(niigata_sst)
monthly_seasonal_cycle_niigata_sst <- mean_seasonal_cycle(niigata_sst) 
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

![](quick_tutorial_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

These utility functions are provided to support data preparation
and　exploratory analysis and are not required for the core modeling
functionality of tempssm.

## References

The statistical modeling framework implemented in `tempssm` is based on
the methodology described in Baba et al. (2024), with implementation
details adapted from the accompanying supplementary materials and code
repository.

Baba, S., Ishii, H., and Yoshiyama, T. (2024). Estimating sea
temperature trends using a linear Gaussian state-space model in
Jogashima, Kanagawa, Japan. *Bulletin of the Japanese Society of
Fisheries Oceanography*, 88(3), 190–199. (In Japanese with an English
abstract.) <https://doi.org/10.34423/jsfo.88.3_190>

Baba, S. (2024). Supplementary code and test data for estimating sea
temperature trends using a linear Gaussian state-space model. GitHub
repository:
<https://github.com/logics-of-blue/sea-temperature-trend-jogashima>

Helske, J. (2017). KFAS: Exponential Family State Space Models in R.  
*Journal of Statistical Software*, 78(10), 1–39.
<https://doi.org/10.18637/jss.v078.i10>
