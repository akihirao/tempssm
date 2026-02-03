
#' Read and convert a monthly temperature CSV file to a \code{ts} object
#'
#' @details
#' The input CSV file must contain monthly temperature data with
#' three columns: \code{Year}, \code{Month}, and \code{Temp}.
#' The expected format is:
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
#' @param csv A character string specifying the path to a CSV file
#'   containing monthly temperature data.
#'
#' @importFrom readr read_csv
#' @importFrom stats ts
#'
#' @encoding UTF-8
#'
#' @export
monthly_temp_csv2ts <- function(csv) {

  raw_data <- readr::read_csv(csv)

  message(
    "Columns in the input file: ",
    paste(names(raw_data), collapse = ", ")
  )

  required_cols <- c("Year", "Month", "Temp")
  if (!all(required_cols %in% names(raw_data))) {
    stop(
      "The CSV file must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  ts(
    as.matrix(raw_data$Temp),
    start = c(raw_data$Year[1], raw_data$Month[1]),
    frequency = 12
  )

}




#' Convert a daily zoo object to a monthly \code{ts} object
#'
#' @param zoo_obj A \code{zoo} object with daily observations.
#'   The index must be of class \code{Date} or \code{POSIXt}.
#' @param var A character string specifying the name of the variable
#'   to be aggregated (default: \code{"Temp"}).
#' @param na.rm Logical; should missing values be removed before averaging?
#'
#' @return A monthly \code{ts} object with frequency = 12.
#'
#' @examples
#' data(ibaraki_sst)
#' ts_monthly <- zoo_daily2ts_monthly(ibaraki_sst)
#' head(ts_monthly)
#'
#' @importFrom zoo index coredata as.yearmon
#' @importFrom stats ts start aggregate
#'
#' @export
zoo_daily2ts_monthly <- function(zoo_obj, var = "Temp", na.rm = TRUE) {

  if (!inherits(zoo_obj, "zoo")) {
    stop("Input must be a zoo object.", call. = FALSE)
  }

  if (!var %in% colnames(zoo_obj)) {
    stop(
      paste0("Variable '", var, "' not found in the zoo object."),
      call. = FALSE
    )
  }

  idx <- zoo::index(zoo_obj)
  if (!inherits(idx, c("Date", "POSIXt"))) {
    stop("Index of the zoo object must be Date or POSIXt.", call. = FALSE)
  }

  zoo_monthly <- stats::aggregate(
    zoo_obj[, var, drop = FALSE],
    zoo::as.yearmon,
    mean,
    na.rm = na.rm
  )

  ym_start <- stats::start(zoo_monthly)

  ts_monthly <- ts(
    zoo::coredata(zoo_monthly),
    start = c(
      as.integer(format(ym_start, "%Y")),
      as.integer(format(ym_start, "%m"))
    ),
    frequency = 12
  )

  colnames(ts_monthly) <- var
  ts_monthly
}





#' Estimate monthly climatology (mean seasonal cycle)
#'
#' @param temp_ts Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12 (monthly data).
#'
#' @return A tibble with one row per month (January–December) containing
#'   the climatological mean temperature.
#'
#' @importFrom stats cycle
#' @importFrom tibble tibble
#'
#' @encoding UTF-8
#'
#' @export
mean_seasonal_cycle <- function(temp_ts){

  if (!inherits(temp_ts, "ts")) {
    stop("Input must be a 'ts' object.", call. = FALSE)
  }

  if (frequency(temp_ts) != 12) {
    stop("Time series must be monthly (frequency = 12).", call. = FALSE)
  }

  temp  <- as.numeric(temp_ts)
  month <- factor(cycle(temp_ts), levels = 1:12)

  monthly_mean <- tapply(temp, month, mean, na.rm = TRUE)

  tibble(
    Month = 1:12,
    Temperature = as.numeric(monthly_mean)
  )
}



#' Calculate monthly temperature anomalies
#'
#' Monthly temperature anomalies are calculated by subtracting a
#' monthly climatology from each observation. The climatology can be
#' computed either from the full time series or from a user-defined
#' baseline period.
#'
#' @importFrom stats cycle window frequency start
#' @importFrom ThermoSSM mean_seasonal_cycle
#'
#' @param temp_ts Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12 (monthly data).
#'
#' @param baseline Optional numeric vector of length 2 specifying
#'   the reference period for climatology in years
#'   (e.g., \code{c(1981, 2010)}). If \code{NULL}, the full period
#'   is used.
#'
#' @details
#' Monthly climatology is computed using
#' \code{mean_seasonal_cycle()}.
#' If \code{baseline} is provided, climatology is estimated using
#' only data within the specified reference period.
#' Missing values are ignored when calculating climatological means.
#'
#' @return A \code{ts} object of monthly temperature anomalies.
#'
#' @examples
#' temp_ts <- ts(
#'   rnorm(12 * 30, mean = 10),
#'   start = c(1981, 1),
#'   frequency = 12
#' )
#'
#' # Full-period climatology
#' anom_all <- monthly_anomaly(temp_ts)
#'
#' # Baseline climatology (1981–2010)
#' anom_base <- monthly_anomaly(temp_ts, baseline = c(1981, 2010))
#'
#' @export
monthly_anomaly <- function(temp_ts, baseline = NULL) {

  if (!inherits(temp_ts, "ts")) {
    stop("Input must be a 'ts' object.", call. = FALSE)
  }

  if (frequency(temp_ts) != 12) {
    stop("Time series must be monthly (frequency = 12).", call. = FALSE)
  }

  ## ---- Select baseline period --------------------------------------------
  if (is.null(baseline)) {

    ts_base <- temp_ts

  } else {

    if (!is.numeric(baseline) || length(baseline) != 2) {
      stop("baseline must be a numeric vector of length 2 (start_year, end_year).",
           call. = FALSE)
    }

    ts_base <- window(
      temp_ts,
      start = c(baseline[1], 1),
      end   = c(baseline[2], 12)
    )

    if (length(ts_base) == 0) {
      stop("No data available in the specified baseline period.",
           call. = FALSE)
    }
  }

  ## ---- Monthly climatology -----------------------------------------------
  clim_tbl <- ThermoSSM::mean_seasonal_cycle(ts_base)
  clim_vec <- clim_tbl$Temperature

  clim <- clim_vec[cycle(temp_ts)]

  ## ---- Anomaly -----------------------------------------------------------
  ts(
    as.numeric(temp_ts) - clim,
    start = start(temp_ts),
    frequency = frequency(temp_ts)
  )
}





#' Retrieve daily mean sea-surface temperature as a zoo object from JMA 
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the data as a \code{zoo} object indexed by date.
#'
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @import dplyr
#' @import zoo
#' @import lubridate
#'
#' @param sea_area_id
#' Numeric sea area ID. The default is NULL
#' For example, 138 corresponding to the coastal sea off southern Ibaraki.
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
#' sst_138_zoo <- sst_jma2zoo(sea_area_id = 138)
#' head(sst_138_zoo)
#' }
#'
#' @export

sst_jma2zoo <- function(sea_area_id = NULL) {
  url_head <- "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area"
  url <- paste0(url_head, sea_area_id, ".txt")
  
  sst_tidy <- readr::read_csv(url, show_col_types = FALSE) %>%
    dplyr::slice(-dplyr::n()) # filtering out the last line of fudder
  
  colnames(sst_tidy) <- c("Year", "Month", "Day", "areaNo", "flag", "Temp")
  
  # Use '.data' to explicitly reference data-frame columns and avoid R CMD check notes”
  sst_tidy <- sst_tidy %>%
    dplyr::mutate(
      Temp = as.numeric(.data$Temp),
      date = lubridate::make_date(
        year  = as.integer(.data$Year),
        month = as.integer(.data$Month),
        day   = as.integer(.data$Day)
      )
    ) %>%
    dplyr::select(.data$date, .data$Temp, .data$flag) %>%
    dplyr::mutate(
      Temp = dplyr::if_else(.data$Temp <= -999, NA_real_, .data$Temp)
    )
  
  # convert tidy object to zoo object
  sst_zoo <- zoo::zoo(
    x = data.frame(Temp = sst_tidy[["Temp"]]),
    order.by = sst_tidy[["date"]]
  )
  
  return(sst_zoo)
}


#' Retrieve daily mean sea-surface temperature as a monthly ts object from JMA 
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the monthly average data as a \code{ts} object.
#'
#' @importFrom ThermoSSM zoo_daily2ts_monthly sst_jma2zoo
#'
#' @param sea_area_id
#' Numeric sea area ID. The default is NULL
#' For example, 138 corresponding to the coastal sea off southern Ibaraki.
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
#' A \code{ts} object of monthly mean sea-surface temperature
#' with a single column named \code{Temp}.
#'
#' @examples
#' \dontrun{
#' sst_138_ts <- sst_jma2ts(sea_area_id = 138)
#' head(sst_138_ts)
#' }
#'
#' @export
sst_jma2ts <- function(sea_area_id = NULL) {
  
  sst_zoo <- ThermoSSM::sst_jma2zoo(sea_area_id)
  monthly_sst_ts <- ThermoSSM::zoo_daily2ts_monthly(sst_zoo)
  
  return(monthly_sst_ts)
}
