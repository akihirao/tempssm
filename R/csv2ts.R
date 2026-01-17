
#' Function to convert from data.csv to time sieries object
#' @import readr
#'
#' @param csv input csv file
#'
#' @encoding UTF-8
#'
#' @export


csv2ts <- function(csv) {

  raw_data <- readr::read_csv(csv)
  message("Columns in a loading file: ", paste(names(raw_data), collapse=", "))

  Temp_ts <- ts(
    as.matrix(raw_data["Temp"]),
    start = c(raw_data$Year[1], raw_data$Month[1]),
    frequency = 12 # 12か月1周期)
    )

  colnames(Temp_ts) <- "Temp"
  return(Temp_ts)
}





