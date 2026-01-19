#' Daily mean temperature time series off southern Ibaraki Prefecture
#'
#' A daily mean temperature time series observed in the coastal waters
#' off southern Ibaraki Prefecture, Japan. The data are provided as a
#' \code{zoo} object and serve as an example dataset for state-space
#' analysis of daily temperature time series.
#'
#' @format
#' A \code{zoo} object with:
#' \describe{
#'   \item{index}{Date from 1982-01-01 to 2025-12-31}
#'   \item{Temp}{Daily mean temperature (°C)}
#' }
#'
#' @details
#' The variable name is unified as \code{Temp} for consistency with other
#' example datasets in the package. In this dataset, \code{Temp}
#' represents daily mean sea surface temperature. The observation period
#' spans from January 1, 1982 to December 31, 2025. Missing values, if any,
#' are encoded as \code{NA}.
#'
#' @source
#' Japan Meteorological Agency (JMA). Data obtained from the JMA website:
#' \url{https://www.jma.go.jp/jma/indexe.html}
#'
#' @examples
#' data(ibaraki_sst)
#' plot(ibaraki_sst,
#'      ylab = "Temperature (°C)",
#'      main = "Daily mean sea surface temperature off southern Ibaraki Prefecture")
#'
"ibaraki_sst"
