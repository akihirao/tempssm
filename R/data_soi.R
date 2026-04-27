#' Southern Oscillation Index (SOI)
#'
#' The Southern Oscillation Index (SOI) is a standardized index 
#' based on the observed sea level pressure (SLP) differences between 
#' Tahiti and Darwin, Australia. The SOI is one measure of the large-scale 
#' fluctuations in air pressure occurring between the western and 
#' eastern tropical Pacific (i.e., the state of the Southern Oscillation) 
#' during El Niño and La Niña episodes. 
#'
#' @format
#' A \code{ts} object with:
#' \describe{
#'   \item{frequency}{12 (monthly data)}
#'   \item{start}{January 1951}
#'   \item{end}{December 2025}
#' }
#'
#' @details
#' This time series represents the SOI index provided by National Centers 
#' for Environmental Information, National Oceanic and Atmospheric 
#' Administration (NOAA). The observation period spans from
#' January 1951 to December 2025.
#'
#' @source
#' National Centers for Environmental Information, National Oceanic and 
#' Atmospheric Administration (NAOO). Data obtained from the below website:
#' \url{https://www.ncei.noaa.gov/access/monitoring/enso/soi}
#'
#' @examples
#' data(soi)
#' plot(soi, ylab = "SOI index")
#'
"soi"
