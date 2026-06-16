#' Assign variable names to a ts object
#'
#' @description
#' `set_ts_name()` assigns variable (column) names to a time series object
#' of class \code{ts}. The function supports both univariate and multivariate
#' time series and ensures that variable names are handled consistently within
#' the \pkg{tempssm} framework.
#'
#' For univariate series, the input is converted to a single-column \code{ts}
#' object with the specified name. For multivariate series, column names are
#' assigned directly.
#'
#' @param ts_in
#' A time series object of class \code{ts}. May be univariate or multivariate.
#'
#' @param label
#' A character string or character vector specifying variable names.
#' For a univariate series, must be a length-one character string.
#' For a multivariate series, must be either length one (recycled) or the same
#' length as the number of columns in \code{ts_in}.
#'
#' @details
#' This function does not modify the time attributes of the input series
#' (e.g., start time, frequency). It only assigns variable names while preserving
#' the internal structure required by downstream functions such as
#' \code{tempssm()}.
#'
#' @return
#' A \code{ts} object with assigned variable names.
#'
#' @seealso
#' \code{\link{trim_ts_overlap}},
#' \code{\link{split_multi_ts}},
#' \code{\link{tempssm}}
#'
#' @examples
#' ## Univariate example
#' ts_uni <- ts(
#'   rnorm(12 * 10),
#'   start = c(2000, 1),
#'   frequency = 12
#' )
#'
#' ts_uni_named <- set_ts_name(ts_uni, label = "temperature")
#'
#' ## Multivariate example
#' ts_multi <- ts(
#'   matrix(rnorm(240), ncol = 3),
#'   start = c(2000, 1),
#'   frequency = 12
#' )
#'
#' ts_multi_named <- set_ts_name(
#'   ts_multi,
#'   label = c("precip", "solar", "wind")
#' )
#'
#' @export
set_ts_name <- function(ts_in, label) {
  ## ---- input check ----------------------------------------------------
  if (!inherits(ts_in, "ts")) {
    cli::cli_abort(
      "`ts_in` must be an object of class {.cls ts}."
    )
  }

  n_col <- NCOL(ts_in)

  if (!is.character(label)) {
    cli::cli_abort("`label` must be a character vector.")
  }

  if (!(length(label) == 1L || length(label) == n_col)) {
    cli::cli_abort(
      "Length of {.arg label} must be 1 or equal to the number of series in {.arg ts_in}."
    )
  }

  .tempssm_cli_debug("Assigning variable names to {.cls ts} object")

  ## ---- recycle label --------------------------------------------------
  if (length(label) == 1L) {
    .tempssm_cli_debug("Recycling label to {n_col} column{?s}")
    label <- rep(label, n_col)
  }

  ## ---- ensure matrix form --------------------------------------------
  x <- if (is.null(dim(ts_in))) {
    .tempssm_cli_debug("Converting univariate series to matrix form")
    matrix(ts_in, ncol = 1)
  } else {
    ts_in
  }

  ## ---- assign names ---------------------------------------------------
  colnames(x) <- label

  ## ---- reconstruct ts -------------------------------------------------
  ts_out <- stats::ts(
    x,
    start = stats::start(ts_in),
    frequency = stats::frequency(ts_in)
  )

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Assigned {length(label)} variable name{?s} to time series"
  )

  return(ts_out)
}


#' Trim and align temperature and exogenous time series over their shared period
#'
#' @description
#' `trim_ts_overlap()` aligns a temperature time series and one or more exogenous
#' time series by trimming them to their shared (overlapping) time period.
#' The function returns the trimmed series as a named list of `ts` objects,
#' with consistent variable labels applied for downstream modeling in
#' \code{tempssm()}.
#'
#' Univariate `ts` objects are labeled using \code{set_ts_name()} to ensure a
#' consistent handling of variable names within the \pkg{tempssm} framework.
#'
#' @param temp_ts
#' A univariate \code{ts} object representing the temperature time series.
#'
#' @param exo_ts
#' A multivariate or univariate \code{ts} object of exogenous variables.
#' Each column represents a distinct exogenous covariate.
#'
#' @param temp_name
#' Character string giving the variable name for the temperature series.
#' This name is applied using \code{set_ts_name()}.
#'
#' @param exo_name
#' Optional character vector giving variable names for the exogenous variables.
#' If \code{NULL}, default names of the form \code{var1}, \code{var2}, \dots
#' are assigned with a warning.
#'
#' @details
#' The shared time window is determined using \code{ts.intersect()}, and both
#' temperature and exogenous series are trimmed accordingly.
#'
#' For univariate series, variable names are assigned via
#' \code{set_ts_name()} rather than \code{colnames()} in order to maintain
#' internal consistency required by \code{tempssm()}.
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'   \item{temperature}{A univariate \code{ts} object of the trimmed temperature series.}
#'   \item{exogenous}{A \code{ts} object of the trimmed exogenous variables.}
#' }
#'
#' @seealso
#' [ts.intersect()], [set_ts_name()], \code{\link{tempssm}}
#'
#' @examples
#' temp_ts <- ts(rnorm(100), start = c(2000, 1), frequency = 12)
#' exo_ts <- ts(matrix(rnorm(200), ncol = 2),
#'   start = c(2001, 1), frequency = 12
#' )
#'
#' trim_ts_overlap(
#'   temp_ts,
#'   exo_ts,
#'   temp_name = "mean_temp",
#'   exo_name  = c("precip", "solar")
#' )
#'
#' @importFrom stats ts.intersect
#' @export
trim_ts_overlap <- function(
  temp_ts,
  exo_ts,
  temp_name = "temp",
  exo_name = NULL
) {
  ## ---- basic checks ---------------------------------------------------
  if (!inherits(temp_ts, "ts")) {
    cli::cli_abort("`temp_ts` must be an object of class {.cls ts}.")
  }

  if (!inherits(exo_ts, "ts")) {
    cli::cli_abort("`exo_ts` must be an object of class {.cls ts}.")
  }

  num_exo_variable <- NCOL(exo_ts)

  ## ---- exo name handling ----------------------------------------------
  if (is.null(exo_name)) {
    cli::cli_warn(
      "`exo_name` not supplied; assigning default names: var1, var2, ..."
    )
    exo_name <- paste0("var", seq_len(num_exo_variable))
  } else {
    if (length(exo_name) != num_exo_variable) {
      cli::cli_abort(
        "Length of {.arg exo_name} must equal the number of exogenous variables."
      )
    }
  }

  .tempssm_cli_debug(
    "Trimming time series overlap (temp length={length(temp_ts)}, exo vars={num_exo_variable})"
  )

  ## ---- overlap --------------------------------------------------------
  temp_exo_ts_overlap <- stats::ts.intersect(temp_ts, exo_ts)

  if (is.null(temp_exo_ts_overlap) || NROW(temp_exo_ts_overlap) == 0) {
    cli::cli_abort(
      "No overlapping time period found between {.arg temp_ts} and {.arg exo_ts}."
    )
  }

  temp_ts_overlap <- temp_exo_ts_overlap[, 1, drop = FALSE]
  exo_ts_overlap <- temp_exo_ts_overlap[, -1, drop = FALSE]

  .tempssm_cli_debug("Overlap length: {NROW(temp_ts_overlap)}")

  ## ---- labels ---------------------------------------------------------
  temp_ts_overlap <- set_ts_name(temp_ts_overlap, label = temp_name)

  if (num_exo_variable == 1) {
    exo_ts_overlap <- set_ts_name(exo_ts_overlap, label = exo_name)
  } else {
    colnames(exo_ts_overlap) <- exo_name
  }

  ## ---- output ---------------------------------------------------------
  out <- list(
    temperature = temp_ts_overlap,
    exogenous   = exo_ts_overlap
  )

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Trimmed series to {NROW(temp_ts_overlap)} overlapping observation{?s}"
  )

  return(out)
}


#' Split a multivariate ts object into univariate ts objects
#'
#' @description
#' `split_multi_ts()` splits a multivariate \code{ts} object into a list of
#' univariate \code{ts} objects, one for each variable (column).
#'
#' Each resulting univariate series preserves the original time attributes and
#' is labeled using \code{set_ts_name()} with the corresponding column name.
#' This ensures consistent handling of variable names within the
#' \pkg{tempssm} framework.
#'
#' @param multi_ts
#' A multivariate \code{ts} object.
#' Each column represents a distinct variable.
#'
#' @details
#' The function requires \code{multi_ts} to have valid column names.
#' If a univariate \code{ts} object is supplied, the function stops with an error.
#'
#' @return
#' A named list of univariate \code{ts} objects.
#' Each list element corresponds to one column of \code{multi_ts}, and the list
#' names are taken from \code{colnames(multi_ts)}.
#'
#' @seealso
#' \code{\link{trim_ts_overlap}},
#' \code{\link{set_ts_name}},
#' \code{\link{tempssm}}
#'
#' @examples
#' multi_ts <- ts(
#'   matrix(rnorm(200), ncol = 2),
#'   start = c(2001, 1),
#'   frequency = 12
#' )
#' colnames(multi_ts) <- c("precip", "solar")
#'
#' split_multi_ts(multi_ts)
#'
#' @export
split_multi_ts <- function(multi_ts) {
  ## ---- input checks ---------------------------------------------------
  if (!inherits(multi_ts, "ts")) {
    cli::cli_abort(
      "`multi_ts` must be an object of class {.cls ts}."
    )
  }

  num_variable <- NCOL(multi_ts)

  if (num_variable == 1) {
    cli::cli_abort(
      "`multi_ts` must be multivariate; a univariate series cannot be split."
    )
  }

  if (is.null(colnames(multi_ts))) {
    cli::cli_abort(
      "`multi_ts` must have column names."
    )
  }

  names_var <- colnames(multi_ts)

  .tempssm_cli_debug(
    "Splitting multivariate ts with {num_variable} variable{?s}"
  )

  ## ---- split ----------------------------------------------------------
  out <- vector("list", num_variable)

  for (i in seq_len(num_variable)) {
    target_name <- names_var[i]

    .tempssm_cli_debug("Processing variable: {target_name}")

    out[[i]] <- set_ts_name(
      multi_ts[, i, drop = FALSE],
      label = target_name
    )
  }

  names(out) <- names_var

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Split time series into {num_variable} univariate series"
  )

  return(out)
}


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
  ## ---- input check ----------------------------------------------------
  if (!is.data.frame(df)) {
    cli::cli_abort(
      "`df` must be a data frame."
    )
  }

  required_cols <- c("Date", "Temp")

  if (!all(required_cols %in% names(df))) {
    cli::cli_abort(
      "The data frame must contain columns: {paste(required_cols, collapse = ', ')}."
    )
  }

  if (!inherits(df$Date, "Date")) {
    cli::cli_abort(
      "The {.col Date} column must be of class {.cls Date}."
    )
  }

  .tempssm_cli_debug("Validating input data frame with {nrow(df)} rows")

  ## ---- sort -----------------------------------------------------------
  df <- df[order(df$Date), ]

  ## ---- check regularity ----------------------------------------------
  ym_index <- as.integer(format(df$Date, "%Y")) * 12 +
    as.integer(format(df$Date, "%m"))

  if (any(diff(ym_index) != 1)) {
    cli::cli_warn(
      "Input data are not strictly monthly; some months may be missing."
    )
  }

  ## ---- construct ts ---------------------------------------------------
  start_year <- as.integer(format(df$Date[1], "%Y"))
  start_month <- as.integer(format(df$Date[1], "%m"))

  .tempssm_cli_debug(
    "Creating ts object starting at {start_year}-{start_month}"
  )

  ts_out <- stats::ts(
    as.numeric(df$Temp),
    start = c(start_year, start_month),
    frequency = 12
  )

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Converted data frame to monthly time series with {length(ts_out)} observation{?s}"
  )

  return(ts_out)
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
#' @return
#' \code{ts} object representing the monthly time series
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
  ## ---- input check ----------------------------------------------------
  if (!is.character(csv) || length(csv) != 1) {
    cli::cli_abort("`csv` must be a file path (character scalar).")
  }

  if (!file.exists(csv)) {
    cli::cli_abort("File does not exist: {csv}")
  }

  .tempssm_cli_debug("Reading CSV file: {csv}")

  ## ---- read -----------------------------------------------------------
  raw_data <- tryCatch(
    readr::read_csv(csv, show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort("Failed to read CSV file: {conditionMessage(e)}")
    }
  )

  .tempssm_cli_debug(
    "Columns in input: {paste(names(raw_data), collapse = ', ')}"
  )

  ## ---- column check ---------------------------------------------------
  required_cols <- c("Year", "Month", "Temp")

  if (!all(required_cols %in% names(raw_data))) {
    cli::cli_abort(
      "The CSV file must contain columns: {paste(required_cols, collapse = ', ')}."
    )
  }

  ## ---- basic validation -----------------------------------------------
  if (nrow(raw_data) == 0) {
    cli::cli_abort("The input CSV file contains no rows.")
  }

  ## ---- optional warning -----------------------------------------------
  if (any(raw_data$Month < 1 | raw_data$Month > 12, na.rm = TRUE)) {
    cli::cli_warn("Detected invalid month values outside 1-12.")
  }

  ## ---- construct ts ---------------------------------------------------
  start_year <- raw_data$Year[1]
  start_month <- raw_data$Month[1]

  .tempssm_cli_debug(
    "Creating ts object starting at {start_year}-{start_month}"
  )

  ts_out <- stats::ts(
    as.numeric(raw_data$Temp),
    start = c(start_year, start_month),
    frequency = 12
  )

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Loaded monthly temperature series with {length(ts_out)} observation{?s}"
  )

  return(ts_out)
}


#' Convert a daily zoo object to a monthly \code{ts} object
#'
#' @param zoo_obj A \code{zoo} object with daily observations.
#'   The index must be of class \code{Date} or \code{POSIXt}.
#'
#' @param var A character string specifying the name of the variable
#'   to be aggregated (default: \code{"Temp"}).
#'
#' @param na.rm Logical; should missing values be removed before averaging?
#'
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
  ## ---- input check ----------------------------------------------------
  if (!inherits(zoo_obj, "zoo")) {
    cli::cli_abort("Input must be a {.cls zoo} object.")
  }

  if (!var %in% colnames(zoo_obj)) {
    cli::cli_abort(
      "Variable {.val {var}} not found in the zoo object."
    )
  }

  if (!is.numeric(na_prop_max) || na_prop_max < 0 || na_prop_max > 1) {
    cli::cli_abort("`na_prop_max` must be between 0 and 1.")
  }

  idx <- zoo::index(zoo_obj)
  if (!inherits(idx, c("Date", "POSIXt"))) {
    cli::cli_abort(
      "Index of the zoo object must be {.cls Date} or {.cls POSIXt}."
    )
  }

  .tempssm_cli_debug(
    "Aggregating daily zoo data to monthly ts (var={var})"
  )

  ## ---- custom monthly summary ----------------------------------------
  monthly_fun <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }

    na_prop <- mean(is.na(x))

    if (na_prop > na_prop_max) {
      return(NA_real_)
    }

    mean(x, na.rm = na.rm)
  }

  ## ---- aggregation ----------------------------------------------------
  zoo_monthly <- stats::aggregate(
    zoo_obj[, var, drop = FALSE],
    zoo::as.yearmon,
    monthly_fun
  )

  if (NROW(zoo_monthly) == 0) {
    cli::cli_abort(
      "No valid observations available after aggregation."
    )
  }

  ## ---- construct ts ---------------------------------------------------
  ym_start <- stats::start(zoo_monthly)

  start_year <- as.integer(format(ym_start, "%Y"))
  start_month <- as.integer(format(ym_start, "%m"))

  .tempssm_cli_debug(
    "Creating monthly ts starting at {start_year}-{start_month}"
  )

  ts_monthly <- stats::ts(
    zoo::coredata(zoo_monthly),
    start = c(start_year, start_month),
    frequency = 12
  )

  colnames(ts_monthly) <- var

  ## ---- optional warning (data quality) -------------------------------
  na_ratio <- mean(is.na(ts_monthly))
  if (na_ratio > 0.3) {
    cli::cli_warn(
      "More than 30% of aggregated monthly values are missing."
    )
  }

  ## ---- inform ---------------------------------------------------------
  .tempssm_cli_inform(
    "Aggregated daily data to monthly series with {NROW(ts_monthly)} observation{?s}"
  )

  return(ts_monthly)
}


#' Compute monthly climatology (mean seasonal cycle)
#'
#' @param temp_ts Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12 (monthly data).
#'
#' @return A tibble with one row per month (January-December) containing
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
compute_monthly_climatology <- function(temp_ts) {
  ## ---- input check ----------------------------------------------------
  if (!inherits(temp_ts, "ts")) {
    cli::cli_abort("Input must be a {.cls ts} object.")
  }

  if (stats::frequency(temp_ts) != 12) {
    cli::cli_abort("Time series must be monthly (frequency = 12).")
  }

  .tempssm_cli_debug(
    "Computing monthly climatology for {length(temp_ts)} observations"
  )

  ## ---- compute --------------------------------------------------------
  temp <- as.numeric(temp_ts)
  month <- factor(stats::cycle(temp_ts), levels = 1:12)

  monthly_mean <- tapply(temp, month, mean, na.rm = TRUE)

  ## ---- output ---------------------------------------------------------
  out <- tibble::tibble(
    Month = 1:12,
    Temperature = as.numeric(monthly_mean)
  )

  .tempssm_cli_debug("Computed climatology for 12 months")

  return(out)
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
#' # Baseline climatology (1981-2010)
#' anom_base <- compute_temp_anomaly(temp_ts, baseline = c(1981, 2010))
#'
#' @export
compute_temp_anomaly <- function(temp_ts, baseline = NULL) {
  ## ---- input check ----------------------------------------------------
  if (!inherits(temp_ts, "ts")) {
    cli::cli_abort("Input must be a {.cls ts} object.")
  }

  if (stats::frequency(temp_ts) != 12) {
    cli::cli_abort("Time series must be monthly (frequency = 12).")
  }

  ## ---- baseline selection ---------------------------------------------
  if (is.null(baseline)) {
    .tempssm_cli_debug("Using full period for climatology")
    ts_base <- temp_ts
  } else {
    if (!is.numeric(baseline) || length(baseline) != 2) {
      cli::cli_abort(
        "`baseline` must be a numeric vector of length 2 (start_year, end_year)."
      )
    }

    .tempssm_cli_debug(
      "Using baseline period: {baseline[1]}-{baseline[2]}"
    )

    ## ---- baseline validation --------------------------------------------
    ts_start <- stats::start(temp_ts)[1]
    ts_end <- stats::end(temp_ts)[1]

    if (baseline[1] > baseline[2]) {
      cli::cli_abort("`baseline` start year must be <= end year.")
    }

    if (baseline[2] < ts_start || baseline[1] > ts_end) {
      cli::cli_abort(
        "No data available in the specified baseline period."
      )
    }

    # Check for overlap
    if (max(baseline[1], ts_start) > min(baseline[2], ts_end)) {
      cli::cli_abort("No data available in the specified baseline period.")
    }

    ts_base <- stats::window(
      temp_ts,
      start = c(baseline[1], 1),
      end   = c(baseline[2], 12)
    )

    if (length(ts_base) == 0) {
      cli::cli_abort(
        "No data available in the specified baseline period."
      )
    }
  }

  ## ---- climatology ----------------------------------------------------
  .tempssm_cli_debug("Computing monthly climatology")

  clim_tbl <- tempssm::compute_monthly_climatology(ts_base)
  clim_vec <- clim_tbl$Temperature

  clim <- clim_vec[stats::cycle(temp_ts)]

  ## ---- anomaly --------------------------------------------------------
  .tempssm_cli_debug("Computing anomalies")

  ts_out <- stats::ts(
    as.numeric(temp_ts) - clim,
    start = stats::start(temp_ts),
    frequency = stats::frequency(temp_ts)
  )

  return(ts_out)
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


#' Parse JMA SST CSV raw data into tidy format
#'
#' Internal helper to parse raw CSV data retrieved from the JMA SST
#' endpoint. The function converts the raw content into a tidy data
#' frame with proper date handling and missing value treatment.
#'
#' @param raw A raw vector containing CSV data.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{date}{Date object}
#'   \item{Temp}{Numeric temperature (NA for missing values)}
#'   \item{flag}{Original flag column}
#' }
#'
#' @keywords internal
#' @noRd
.parse_jma_csv <- function(raw) {
  sst_raw <- readr::read_csv(
    raw,
    show_col_types = FALSE
  )

  if (nrow(sst_raw) < 2) {
    cli::cli_abort("Retrieved data has unexpected or empty format.")
  }

  ## remove last row (END or footer)
  sst_tidy <- sst_raw %>%
    dplyr::slice(-dplyr::n())

  colnames(sst_tidy) <- c(
    "Year", "Month", "Day", "areaNo", "flag", "Temp"
  )

  sst_tidy <- sst_tidy %>%
    dplyr::mutate(
      Temp = as.numeric(.data$Temp),
      date = lubridate::make_date(
        year  = as.integer(.data$Year),
        month = as.integer(.data$Month),
        day   = as.integer(.data$Day)
      ),
      Temp = dplyr::if_else(.data$Temp <= -999, NA_real_, .data$Temp)
    ) %>%
    dplyr::select(date, Temp, flag)

  sst_tidy
}


#' Construct JMA SST endpoint URL
#'
#' Internal helper to generate the HTTP endpoint URL for retrieving
#' JMA sea-surface temperature data based on a sea area ID.
#'
#' @param sea_area_id Character string or numeric JMA sea area ID.
#'
#' @return A character string representing the full URL.
#'
#' @details
#' This function performs only string construction and does not validate
#' the existence or accessibility of the endpoint.
#'
#' @keywords internal
#' @noRd
.build_jma_url <- function(sea_area_id) {
  paste0(
    "https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/txt/area",
    sea_area_id,
    ".txt"
  )
}


#' Fetch raw SST data from JMA endpoint
#'
#' Internal helper to perform an HTTP request to the JMA SST endpoint
#' and return the raw response body. This function encapsulates
#' request construction, User-Agent handling, and error management.
#'
#' The returned value is a raw vector intended to be parsed by
#' \code{.parse_jma_csv()}.
#'
#' @param sea_area_id Character string or numeric JMA sea area ID.
#'
#' @return A raw vector containing CSV data retrieved from the endpoint.
#'
#' @details
#' This function is responsible only for data retrieval and does not
#' perform any parsing or validation of the content.
#'
#' @keywords internal
#' @noRd
.fetch_jma_raw <- function(sea_area_id) {
  if (!is.character(sea_area_id)) {
    sea_area_id <- as.character(sea_area_id)
  }

  url <- .build_jma_url(sea_area_id)

  .tempssm_cli_debug("Fetching JMA SST data (area_id={sea_area_id})")

  resp <- tryCatch(
    httr2::request(url) %>%
      httr2::req_user_agent(.user_agent()) %>%
      httr2::req_error(body = function(resp) {
        paste0(
          "Failed to retrieve SST data for sea_area_id = '",
          sea_area_id,
          "'. HTTP status: ",
          httr2::resp_status(resp)
        )
      }) %>%
      httr2::req_perform(),
    error = function(e) {
      cli::cli_abort(
        "Failed to access JMA SST endpoint: {conditionMessage(e)}"
      )
    }
  )

  .tempssm_cli_debug("HTTP request successful")

  httr2::resp_body_raw(resp)
}


#' Check proportion of missing values and issue warning
#'
#' Internal helper to compute the proportion of missing values in a vector
#' and issue a warning if it exceeds a specified threshold.
#'
#' @param x A numeric vector.
#' @param threshold Numeric threshold for the proportion of missing values.
#'   A warning is triggered if \code{mean(is.na(x))} exceeds this value.
#' @param msg Character string; warning message to display.
#'
#' @return Invisibly returns \code{NULL}. This function is called for its
#'   side effect (issuing a warning).
#'
#' @keywords internal
#' @noRd
.check_na_ratio <- function(x, threshold = 0.3, msg) {
  na_ratio <- mean(is.na(x))

  if (na_ratio > threshold) {
    cli::cli_warn(msg)
  }
}


#' Convert tidy SST data to zoo object
#'
#' Internal helper to construct a \code{zoo} time series object from
#' a tidy data frame containing SST observations. The function maps
#' the \code{Temp} column to values and the \code{date} column to the
#' time index.
#'
#' @param sst_tidy A data frame containing at least two columns:
#'   \code{date} (Date) and \code{Temp} (numeric).
#'
#' @return A \code{zoo} object with a single column named \code{Temp}
#'   indexed by \code{date}.
#'
#' @keywords internal
#' @noRd
.build_sst_zoo <- function(sst_tidy) {
  zoo::zoo(
    x = data.frame(Temp = sst_tidy$Temp),
    order.by = sst_tidy$date
  )
}


#' Retrieve daily mean sea-surface temperature as a zoo object from JMA
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the data as a \code{zoo} object indexed by date.
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @importFrom httr2 request req_user_agent req_error req_perform resp_body_raw resp_status
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
  raw <- .fetch_jma_raw(sea_area_id)

  sst_tidy <- .parse_jma_csv(raw)

  .check_na_ratio(
    sst_tidy$Temp,
    msg = "More than 30% of SST values are missing in retrieved data."
  )

  out <- .build_sst_zoo(sst_tidy)

  .tempssm_cli_inform(
    "Retrieved SST data with {nrow(sst_tidy)} observation{?s} for area {sea_area_id}"
  )

  out
}


#' Retrieve daily mean sea-surface temperature as a monthly ts object from JMA
#'
#' This function downloads publicly available daily mean
#' sea-surface temperature (SST) data for Japanese coastal waters
#' provided by the Japan Meteorological Agency (JMA),
#' and returns the monthly average data as a \code{ts} object.
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @importFrom httr2 request req_user_agent req_error req_perform resp_body_raw resp_status
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
#'
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
get_jma_sst_ts <- function(sea_area_id, na_prop_max = 1) {
  if (!is.numeric(na_prop_max) || na_prop_max < 0 || na_prop_max > 1) {
    cli::cli_abort("`na_prop_max` must be between 0 and 1.")
  }

  raw <- .fetch_jma_raw(sea_area_id)

  sst_tidy <- .parse_jma_csv(raw)

  .check_na_ratio(
    sst_tidy$Temp,
    msg = "More than 30% of daily SST values are missing."
  )

  sst_zoo <- .build_sst_zoo(sst_tidy)

  .tempssm_cli_debug("Aggregating daily data to monthly time series")

  monthly_sst_ts <- tempssm::aggregate_daily_zoo_to_monthly_ts(
    sst_zoo,
    na_prop_max = na_prop_max
  )

  .tempssm_cli_inform(
    "Retrieved and aggregated SST data to {length(monthly_sst_ts)} monthly observation{?s} (area {sea_area_id})"
  )

  monthly_sst_ts
}
