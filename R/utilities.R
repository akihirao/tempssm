#' Convert a data frame of monthly temperature time series to a \code{ts} object
#'
#' @details
#' This function converts a data frame containing a monthly temperature
#' time series into an R \code{ts} object with a frequency of 12.
#'
#' The input data frame must contain the following columns:
#' \describe{
#'   \item{\code{Date}}{A date column indicating the time index.
#'   If the exact day of the month is not uniquely defined, any arbitrary
#'   day (e.g., the first day of the month) may be used.}
#'   \item{\code{Temp}}{Monthly temperature values. Missing values should be
#'   represented as \code{NA}.}
#' }
#'
#' The function assumes that the input data represent a regularly spaced
#' monthly time series. If missing months are detected in the \code{Date}
#' column, a warning is issued, but the \code{ts} object is still created.
#'
#' The expected format is:
#'
#' \preformatted{
#' Date        Temp
#' 2001-01-01  10.4
#' 2001-02-01   8.2
#' 2001-03-01   NA
#' 2001-04-01  13.6
#' ...
#' }
#'
#' @param df A data frame containing monthly temperature data with columns
#'   \code{Date} and \code{Temp}.
#'
#' @return A univariate \code{ts} object representing the monthly temperature
#'   time series.
#'
#' @importFrom stats ts
#'
#' @encoding UTF-8
#'
#' @examples
#' ## Create a data frame of monthly temperature data
#' df <- data.frame(
#'   Date = as.Date(c(
#'     "2001-01-01",
#'     "2001-02-01",
#'     "2001-03-01",
#'     "2001-04-01",
#'     "2001-05-01"
#'   )),
#'   Temp = c(10.4, 8.2, NA, 13.6, 16.1)
#' )
#'
#' ## Convert to a monthly ts object
#' temp_ts <- convert_monthly_df_to_ts(df)
#'
#' ## Inspect the result
#' temp_ts
#'
#' @export
convert_monthly_df_to_ts <- function(df) {

  # Check required columns
  required_cols <- c("Date", "Temp")
  if (!all(required_cols %in% names(df))) {
    stop(
      "The data frame must contain the following columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Ensure Date is of class Date
  if (!inherits(df$Date, "Date")) {
    stop("The 'Date' column must be of class Date.", call. = FALSE)
  }

  # Sort by Date to ensure correct temporal order
  df <- df[order(df$Date), ]

  # Check regular monthly spacing
  ym_index <- as.integer(format(df$Date, "%Y")) * 12 +
              as.integer(format(df$Date, "%m"))

  if (any(diff(ym_index) != 1)) {
    warning(
      "The input data do not form a strictly regular monthly time series. ",
      "Some months may be missing.",
      call. = FALSE
    )
  }

  # Create ts object (monthly)
  ts(
    as.numeric(df$Temp),
    start = c(
      as.integer(format(df$Date[1], "%Y")),
      as.integer(format(df$Date[1], "%m"))
    ),
    frequency = 12
  )
}


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
#' @examples
#' ## Create a temporary CSV file with monthly temperature data
#' tmp_csv <- tempfile(fileext = ".csv")
#'
#' writeLines(
#'   c(
#'     "Year,Month,Temp",
#'     "2001,1,10.4",
#'     "2001,2,8.2",
#'     "2001,3,NA",
#'     "2001,4,13.6",
#'     "2001,5,16.1"
#'   ),
#'   tmp_csv
#' )
#'
#' ## Read the CSV file and convert to a monthly ts object
#' temp_ts <- read_monthly_temp_ts(tmp_csv)
#'
#' ## Inspect the result
#' temp_ts
#'
#' @export
read_monthly_temp_ts <- function(csv) {

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
#' @param na_prop_max Maximum allowed proportion of NA values within a month.
#'   If the proportion of missing values exceeds this threshold, the monthly
#'   mean is set to NA. Default is \code{1} (no additional filtering).
#'
#' @return A monthly \code{ts} object with frequency = 12.
#'
#' @examples
#' \dontrun{
#' sst_zoo <- get_jma_sst_zoo(sea_area_id = 138)
#' sst_ts_monthly <- aggregate_daily_zoo_to_monthly_ts(sst_138_zoo)
#' }
#'
#' @importFrom zoo index coredata as.yearmon
#' @importFrom stats ts start aggregate
#'
#' @export
aggregate_daily_zoo_to_monthly_ts <- function(zoo_obj,
                                 var = "Temp",
                                 na.rm = TRUE,
                                 na_prop_max = 1) {
  
  if (!inherits(zoo_obj, "zoo")) {
    stop("Input must be a zoo object.", call. = FALSE)
  }
  
  if (!var %in% colnames(zoo_obj)) {
    stop(
      paste0("Variable '", var, "' not found in the zoo object."),
      call. = FALSE
    )
  }
  
  if (!is.numeric(na_prop_max) || na_prop_max < 0 || na_prop_max > 1) {
    stop("na_prop_max must be between 0 and 1.", call. = FALSE)
  }
  
  idx <- zoo::index(zoo_obj)
  if (!inherits(idx, c("Date", "POSIXt"))) {
    stop("Index of the zoo object must be Date or POSIXt.", call. = FALSE)
  }
  
  # --- custom monthly summary ---
  monthly_fun <- function(x) {
    
    # Check for all of NA in a month
    if (all(is.na(x))) {
      return(NA_real_)
    }
    
    na_prop <- mean(is.na(x))
    
    if (na_prop > na_prop_max) {
      return(NA_real_)
    }
    
    mean(x, na.rm = na.rm)
  }
  
  zoo_monthly <- stats::aggregate(
    zoo_obj[, var, drop = FALSE],
    zoo::as.yearmon,
    monthly_fun
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



#' Compute monthly climatology (mean seasonal cycle)
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
#' @examples
#' temp_ts <- ts(
#'   rnorm(12 * 30, mean = 10),
#'   start = c(1981, 1),
#'   frequency = 12
#' )
#'
#' monthly_climatology <- compute_monthly_climatology(temp_ts)
#' 
#' @export
compute_monthly_climatology <- function(temp_ts){

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



#' Compute monthly temperature anomalies
#'
#' Monthly temperature anomalies are calculated by subtracting a
#' monthly climatology from each observation. The climatology can be
#' computed either from the full time series or from a user-defined
#' baseline period.
#'
#' @importFrom stats cycle window frequency start
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
#' \code{compute_temp_anomaly()}.
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
#' anom_all <- compute_temp_anomaly(temp_ts)
#'
#' # Baseline climatology (1981–2010)
#' anom_base <- compute_temp_anomaly(temp_ts, baseline = c(1981, 2010))
#' 
#' @export
compute_temp_anomaly <- function(temp_ts, baseline = NULL) {

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
  clim_tbl <- tempssm::compute_monthly_climatology(ts_base)
  clim_vec <- clim_tbl$Temperature

  clim <- clim_vec[cycle(temp_ts)]

  ## ---- Anomaly -----------------------------------------------------------
  ts(
    as.numeric(temp_ts) - clim,
    start = start(temp_ts),
    frequency = frequency(temp_ts)
  )
}



# internal utilities ------------------------------------
.user_agent <- function() {
  getOption(
    "tempssm.user_agent",
    paste0(
      "tempssm/",
      utils::packageVersion("tempssm"),
      " (https://github.com/yourname/tempssm)"
    )
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
#' @importFrom httr2
#'   request
#'   req_user_agent
#'   req_error
#'   req_perform
#'   resp_body_raw
#'   resp_status
#' @import dplyr
#' @import zoo
#' @import lubridate
#'
#' @param sea_area_id
#' Character string giving the JMA sea area ID
#' (numeric values are accepted and internally coerced to character).
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
#' To comply with good API usage practices, HTTP requests sent by this
#' function include a custom \emph{User-Agent} string identifying the
#' \pkg{tempssm} package. Users may override the default User-Agent by
#' setting the \code{"tempssm.user_agent"} option, for example:
#' \code{options(tempssm.user_agent = "my-analysis/1.0")}.
#'
#' @return
#' A \code{zoo} object of daily mean sea-surface temperature
#' with a single column named \code{Temp}.
#'
#' @examples
#' \dontrun{
#' sst_138_zoo <- get_jma_sst_zoo(sea_area_id = 138)
#' head(sst_138_zoo)
#' }
#'
#' @export
get_jma_sst_zoo <- function(sea_area_id) {

  if (!is.character(sea_area_id)) {
    sea_area_id <- as.character(sea_area_id)
  }

  url <- paste0(
    "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area",
    sea_area_id,
    ".txt"
  )

  # ---- HTTP request ----
  resp <- httr2::request(url) |>
    httr2::req_user_agent(.user_agent()) |>
    httr2::req_error(body = function(resp) {
      paste0(
        "Failed to retrieve SST data for sea_area_id = '",
        sea_area_id,
        "'.\nHTTP status: ",
        httr2::resp_status(resp)
      )
    }) |>
    httr2::req_perform()

  # ---- parse ----
  sst_raw <- readr::read_csv(
    httr2::resp_body_raw(resp),
    show_col_types = FALSE
  )

  if (nrow(sst_raw) < 2) {
    stop("Retrieved data has unexpected format.", call. = FALSE)
  }

  sst_tidy <- sst_raw |>
    dplyr::slice(-dplyr::n())

  colnames(sst_tidy) <- c("Year", "Month", "Day", "areaNo", "flag", "Temp")

  sst_tidy <- sst_tidy |>
    dplyr::mutate(
      Temp = as.numeric(.data$Temp),
      date = lubridate::make_date(
        year  = as.integer(.data$Year),
        month = as.integer(.data$Month),
        day   = as.integer(.data$Day)
      ),
      Temp = dplyr::if_else(.data$Temp <= -999, NA_real_, .data$Temp)
    ) |>
    dplyr::select(.data$date, .data$Temp, .data$flag)

  zoo::zoo(
    x = data.frame(Temp = sst_tidy$Temp),
    order.by = sst_tidy$date
  )
}



#' Retrieve daily mean sea-surface temperature as a monthly ts object from JMA
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the monthly average data as a \code{ts} object.
#'
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @importFrom httr2
#'   request
#'   req_user_agent
#'   req_error
#'   req_perform
#'   resp_body_raw
#'   resp_status
#' @import dplyr
#' @import zoo
#' @import lubridate
#'
#' @param sea_area_id
#' Character string giving the JMA sea area ID
#' (numeric values are accepted and internally coerced to character).
#' For example, "138" corresponding to the coastal sea off southern Ibaraki.
#' A list of sea area IDs and their corresponding regions is available at:
#' \url{https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html}

#' @param na_prop_max Maximum allowed proportion of NA values within a month.
#'   If the proportion of missing values exceeds this threshold, the monthly
#'   mean is set to NA. Default is \code{1} (no additional filtering).
#'
#' @details
#' The function retrieves a text-format dataset from the JMA website,
#' parses daily observations, and constructs a \code{zoo} object
#' with calendar dates as the index.
#' Missing values in the original dataset are represented as \code{NA}.
#'
#' To comply with good API usage practices, HTTP requests sent by this
#' function include a custom \emph{User-Agent} string identifying the
#' \pkg{tempssm} package. Users may override the default User-Agent by
#' setting the \code{"tempssm.user_agent"} option, for example:
#' \code{options(tempssm.user_agent = "my-analysis/1.0")}.
#'
#' @return
#' A \code{ts} object of monthly mean sea-surface temperature
#' with a single column named \code{Temp}.
#'
#' @examples
#' \dontrun{
#' sst_138_ts <- get_jma_sst_ts(sea_area_id = 138)
#' head(sst_138_ts)
#' }
#'
#' @export
get_jma_sst_ts <- function(sea_area_id,
                       na_prop_max = 1) {

  if (!is.character(sea_area_id)) {
    sea_area_id <- as.character(sea_area_id)
  }

  url <- paste0(
    "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area",
    sea_area_id,
    ".txt"
  )

  # ---- HTTP request ----
  resp <- httr2::request(url) |>
    httr2::req_user_agent(.user_agent()) |>
    httr2::req_error(body = function(resp) {
      paste0(
        "Failed to retrieve SST data for sea_area_id = '",
        sea_area_id,
        "'.\nHTTP status: ",
        httr2::resp_status(resp)
      )
    }) |>
    httr2::req_perform()

  # ---- parse ----
  sst_raw <- readr::read_csv(
    httr2::resp_body_raw(resp),
    show_col_types = FALSE
  )

  if (nrow(sst_raw) < 2) {
    stop("Retrieved data has unexpected format.", call. = FALSE)
  }

  sst_tidy <- sst_raw |>
    dplyr::slice(-dplyr::n())

  colnames(sst_tidy) <- c("Year", "Month", "Day", "areaNo", "flag", "Temp")

  sst_tidy <- sst_tidy |>
    dplyr::mutate(
      Temp = as.numeric(.data$Temp),
      date = lubridate::make_date(
        year  = as.integer(.data$Year),
        month = as.integer(.data$Month),
        day   = as.integer(.data$Day)
      ),
      Temp = dplyr::if_else(.data$Temp <= -999, NA_real_, .data$Temp)
    ) |>
    dplyr::select(.data$date, .data$Temp, .data$flag)

  sst_zoo <- zoo::zoo(
    x = data.frame(Temp = sst_tidy$Temp),
    order.by = sst_tidy$date
  )

  monthly_sst_ts <- tempssm::aggregate_daily_zoo_to_monthly_ts(sst_zoo,
    na_prop_max = na_prop_max)
  return(monthly_sst_ts)
}

