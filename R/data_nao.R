#' North Atlantic Oscillation：NAO (Hurrell)
#'
#' The winter (December thru March) station-based index of the NAO is based on
#' the difference of normalized sea level pressure (SLP) between Lisbon, Portugal
#' and Stykkisholmur/Reykjavik, Iceland since 1864. Positive values of the NAO
#' index are typically associated with stronger-than-average westerlies over
#' the middle latitudes, more intense weather systems over the North Atlantic
#' and wetter/milder weather over western Europe. Monthly, seasonal and annual
#' indices using slightly different data sources for the southern station are
#' also available.
#'
#' @format
#' A \code{ts} object with:
#' \describe{
#'   \item{frequency}{12 (monthly data)}
#'   \item{start}{January 1865}
#'   \item{end}{June 2023}
#' }
#'
#' @details
#' This time series represents the NAO index provided by the Climate
#' Analysis Section, NCAR, Boulder, USA, Hurrell (2003). The observation
#' period spans from January 1865 to June 2023.
#'
#' @source
#' NSF NCAR. Data obtained from the NCAR website:
#' \url{https://climatedataguide.ucar.edu/climate-data/}
#'
#' @examples
#' data(nao)
#' plot(nao)
#'
"nao"
