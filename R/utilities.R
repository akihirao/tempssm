
#' Function to convert from data.csv to monthly time sieries object
#' @import readr
#'
#' @param csv input csv file
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
  return(Temp_ts)
}



#' Function to convert from dairy zoo object to monthly ts object
#' @import zoo
#'
#' @param csv input csv file
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


