#' tempssm: State space models for temperature time series
#'
#' \pkg{tempssm} provides tools for analyzing temperature time series
#' (e.g., air temperature and sea surface temperature) using
#' linear Gaussian state-space models.
#'
#' The package is designed for long-term environmental monitoring data
#' and focuses on trend extraction, seasonal components, and uncertainty
#' estimation under a unified state-space modeling framework.
#'
#' @details
#' Typical applications include:
#' \itemize{
#'   \item Analysis of long-term air temperature records
#'   \item Sea surface temperature trend estimation
#'   \item Decomposition into level, trend, and seasonal components
#' }
#'
#' The implementation is based on linear Gaussian state-space models
#' and makes use of Kalman filtering and smoothing.
#'
#' Main parts of the implementation were adapted from the supplementary code
#' provided in Baba et al. (2024), available on GitHub:
#' https://github.com/logics-of-blue/sea-temperature-trend-jogashima
#'
#'
#' @section Life Cycle and Development Status:
#' `tempssm` is under active development. The current design focuses on
#' temperature time series represented as base R `ts` objects with integer
#' seasonal frequencies greater than 1. The examples and validation primarily
#' use monthly data, but the core modelling path preserves the input
#' frequency. This scope is intentional: it keeps model fitting and
#' cross-validation feasible on typical personal computing environments used
#' by R users.
#'
#' The package aims to stabilize this workflow before expanding the supported
#' data structures. Future versions may consider direct modelling support for
#' daily or irregular time series, such as `zoo` objects, but such extensions
#' are not part of the current design goals and would require careful
#' consideration of computational costs, input validation, and documentation.
#'
#' @section Statistical Terminology:
#' In this package, a linear Gaussian state-space model refers to a model in
#' which observed temperature values are represented as noisy observations of
#' unobserved latent states, with Gaussian observation and state disturbances.
#' The latent states may include a level component, a seasonal component, an
#' autoregressive component, and optional exogenous effects.
#'
#' The level component represents the slowly varying baseline temperature. The
#' seasonal component represents recurring within-year variation. The
#' autoregressive component represents short-term serial dependence not captured
#' by the level or seasonal components. The AR order is the number of lagged
#' autoregressive terms included in this component.
#'
#' Exogenous variables are external time series supplied by the user and aligned
#' with the temperature series. They are included as regression effects in the
#' state-space model. Kalman filtering refers to sequential estimation of the
#' latent state using observations up to each time point, whereas Kalman
#' smoothing refers to estimation of latent states using the full observed time
#' series.
#'
#' Residual diagnostics in this package use standardized recursive residuals
#' extracted from the fitted state-space model. Time-series cross-validation
#' refers to rolling-origin evaluation, where models are fitted on earlier
#' observations and evaluated on later held-out observations. MAE denotes mean
#' absolute error, and MASE denotes mean absolute scaled error.
#'
#' @section Time Index and Calendar Conventions:
#' Core modelling functions use base R `ts` objects. The `frequency` attribute
#' defines the number of observations per seasonal cycle, and the core model
#' preserves the input frequency. For example, `frequency = 12` represents
#' monthly observations, `frequency = 24` represents twice-monthly
#' observations, `frequency = 36` represents three observations per month, and
#' `frequency = 4` represents four seasonal observations per year. The model
#' treats these observations as regularly spaced periods and does not convert
#' periods to a fixed number of days.
#'
#' Daily SST data are handled only in conversion utilities using `zoo` objects
#' indexed by `Date` or `POSIXt`. Daily observations are assigned to calendar
#' months using `zoo::as.yearmon()` before aggregation to monthly `ts` objects.
#' No assumption that a year has 365 or 365.2422 days is used in this
#' conversion path.
#'
#' Cross-validation arguments such as `initial`, `horizon`, and `step` are
#' counts of observations, not calendar durations.
#'
#' @references
#' Helske, J. (2017). KFAS: Exponential family state space models in R.
#' \emph{Journal of Statistical Software}, 78(10), 1--39.
#' \doi{10.18637/jss.v078.i10}
#'
#' Baba, S., Ishii, H., and Yoshiyama, T. (2024). Estimating sea temperature
#' trends using a linear Gaussian state-space model in Jogashima, Kanagawa,
#' Japan (in Japanese with English abstract, tables, and figures).
#' \emph{Bulletin of the Japanese Society of Fisheries Oceanography},
#' 88(3), 190--199. \doi{10.34423/jsfo.88.3_190}
#'
#' Supplementary code:
#' https://github.com/logics-of-blue/sea-temperature-trend-jogashima
#'
#' @seealso
#' \url{https://github.com/akihirao/tempssm}
#'
#' @name tempssm
#' @aliases tempssm-package
NULL
