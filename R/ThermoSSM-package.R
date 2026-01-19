#' ThermoSSM: State space analysis of temperature time series
#'
#' \pkg{ThermoSSM} provides tools for analyzing temperature time series
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
#' @author
#' Akira Hirao
#'
#' @references
#' Baba, S., Ishii, H., and Yoshiyama, T. (2024).
#' Estimating sea temperature trends using a linear Gaussian state-space model
#' in Jogashima, Kanagawa, Japan.
#' Bulletin of the Japanese Society of Fisheries Oceanography, 88(3), 190–199
#' (in Japanese with English abstract, tables, and figures).
#'
#' Supplementary code:
#' https://github.com/logics-of-blue/sea-temperature-trend-jogashima
#'
#' @seealso
#' \url{https://github.com/akihirao/ThermoSSM}
#'
#' @name ThermoSSM
#' @aliases ThermoSSM-package
NULL
