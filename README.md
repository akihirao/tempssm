# ThermoSSM

- R package for state-space analysis of temperature time series.  
The implementation is based on linear Gaussian state-space models and makes use of Kalman filtering and smoothing.  
Main parts of the implementation were adapted from the supplementary code provided in Baba et al. (2024), available on GitHub: https://github.com/logics-of-blue/sea-temperature-trend-jogashima



# Quick tutorial
https://github.com/akihirao/ThermoSSM/blob/main/quick_tutorial.md

# How to install
```
install.packages("devtools")
devtools::install_github("akihirao/ThermoSSM")
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
