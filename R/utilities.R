
#' Convert a monthly temperature CSV file to a monthly time series object
#'
#' @details
#' The input CSV file must have the following format:
#'
#' \preformatted{
#' Year,Month,Temp
#' 2001,1,10.4
#' 2001,2,8.2
#' 2001,3,NA
#' 2001,4,13.6
#' ...
#' }
#'
#' @param csv Path to a CSV file containing monthly temperature data.
#'
#' @importFrom readr read_csv
#'
#' @encoding UTF-8
#'
#' @export
monthly_csv2ts <- function(csv) {

  raw_data <- readr::read_csv(csv)
  message("Columns in a loading file: ", paste(names(raw_data), collapse=", "))

  Temp_ts <- ts(
    as.matrix(raw_data["Temp"]),
    start = c(raw_data$Year[1], raw_data$Month[1]),
    frequency = 12
  )

  colnames(Temp_ts) <- "Temp"
  Temp_ts
}



#' Convert dairy zoo object to monthly ts object
#'
#' @import zoo
#'
#' @param zoo zoo object
#'
#' @encoding UTF-8
#'
#' @export

zoo_daily2ts_monthly <- function(zoo, var = "Temp", na.rm = TRUE) {

  if (!inherits(zoo, "zoo")) {
    stop("Input must be a zoo object.")
  }

  if (!var %in% colnames(zoo)) {
    stop(paste0("Variable '", var, "' not found in zoo object."))
  }

  zoo_monthly <- aggregate(
    zoo[, var, drop = FALSE],
    as.yearmon,
    mean,
    na.rm = na.rm
  )
  colnames(zoo_monthly) <- var

  start_ym <- start(zoo_monthly)
  start_year  <- floor(start_ym)
  start_month <- round((start_ym %% 1) * 12 + 1)

  ts_monthly <- ts(
    coredata(zoo_monthly),
    start = c(start_year, start_month),
    frequency = 12
  )

  return(ts_monthly)

}



#' Generate a monthly temperature anomaly time series
#'
#' Temperature anomalies are calculated by subtracting the long-term
#' monthly mean (climatology) from each observation.
#'
#' @param temp_ts A monthly temperature time series object of class \code{ts}.
#'
#' @details
#' Monthly climatology is calculated as the mean temperature for each calendar
#' month. Missing values are ignored.
#'
#' @return A \code{ts} object of monthly temperature anomalies.
#'
#' @examples
#' temp_ts <- ts(
#'   c(5.2, 6.1, 9.3, 13.5, 17.8, 21.0,
#'     24.3, 25.1, 22.0, 17.1, 11.2, 7.0),
#'   start = c(2000, 1),
#'   frequency = 12
#' )
#'
#' temp_dev <- generate_temp_dev(temp_ts)
#' tapply(as.numeric(temp_dev), cycle(temp_dev), mean)
#'
#' @export

generate_temp_dev <- function(temp_ts) {

  stopifnot(inherits(temp_ts, "ts"))
  stopifnot(frequency(temp_ts) == 12)

  temp <- as.numeric(temp_ts)
  mnum <- cycle(temp_ts)
  mfac <- factor(mnum, levels = 1:12)

  monthly_ave_temp <- tapply(temp, mfac, mean, na.rm = TRUE)
  clim <- monthly_ave_temp[as.integer(mfac)]

  ts(
    temp - clim,
    start = start(temp_ts),
    frequency = frequency(temp_ts)
  )
}



#' Retrieve daily mean sea-surface temperature as a zoo object
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the data as a \code{zoo} object indexed by date.
#'
#' @importFrom readr read_csv
#' @import dplyr
#' @import zoo
#'
#' @param sea_area_id
#' Numeric sea area ID.
#' The default is 138, corresponding to the coastal sea off southern Ibaraki.
#' A list of sea area IDs and their corresponding regions is available at:
#' \url{https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html}
#'
#' @details
#' The function retrieves a text-format dataset from the JMA website,
#' parses daily observations, and constructs a \code{zoo} object
#' with calendar dates as the index.
#' Missing values in the original dataset are represented as \code{NA}.
#'
#' @return
#' A \code{zoo} object of daily mean sea-surface temperature
#' with a single column named \code{Temp}.
#'
#' @examples
#' \dontrun{
#' sst_138 <- sst_jma2zoo(sea_area_id = 138)
#' head(sst_138)
#' }
#'
#' @export

sst_jma2zoo <- function(sea_area_id = 138) {

  url_head <- "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area"
  url <- paste0(url_head, sea_area_id, ".txt")

  sst_tidy <- readr::read_csv(url, show_col_types = FALSE) |>
    dplyr::slice(-n()) |>
    dplyr::rename(Temp = "Temp.") |>
    dplyr::mutate(
      Temp = as.numeric(Temp),
      date = as.Date(paste0(yyyy, "-", mm, "-", dd))
    ) |>
    dplyr::select(date, Temp, flag)

  sst_zoo <- zoo::zoo(
    x = data.frame(Temp = sst_tidy$Temp),
    order.by = sst_tidy$date
  )

  return(sst_zoo)
}

