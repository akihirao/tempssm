
#' Function to plot seasonal trend of monthly temperature
#' @import ggplot2
#' @param ts temperature time series object
#'
#' @encoding UTF-8
#'
#' @export
#'

plot_typical_seasonal_cycle <- function(ts){

  stopifnot(frequency(ts) == 12)  # check for monthly time series

  temp <- as.numeric(ts)
  mnum <- cycle(ts)
  mfac <- factor(mnum, levels = 1:12)
  monthly_ave_temp <- tapply(temp, mfac, function(v) mean(v, na.rm = TRUE))

  monthly_ave_temp_tidy = tibble(Month=seq(1:12),
                                 Temperature=monthly_ave_temp
                                 )

  plot_monthly_ave_temp <- ggplot2::ggplot(data=monthly_ave_temp_tidy,
                                  aes(x=Month,y=Temperature)
                                  ) +
    geom_point(size = 2) +
    geom_line(linetype= "dashed") +
    labs(y = paste0("Monthly average",expression(Temperature~(degree*C)))) +
    scale_x_discrete(
      labels = function(x) sprintf("%02d", as.integer(x))
    )


  plot(plot_monthly_ave_temp )
  #return(plot_monthly_ave_temp)

}



#' Function to plot time series of monthly temperature and temperature deviation
#' @import ggplot2
#' @param ts temperature time series object
#'
#' @encoding UTF-8
#'
#' @export
#'

plot_temp_dev <- function(ts){

  stopifnot(frequency(ts) == 12)  # check for monthly time series

  temp <- as.numeric(ts)
  mnum <- cycle(ts)
  mfac <- factor(mnum, levels = 1:12)
  monthly_ave_temp <- tapply(temp, mfac, function(v) mean(v, na.rm = TRUE))

  clim <- monthly_ave_temp[as.integer(mfac)]
  anom <- temp - clim

  ts_dev <- cbind(Temp=ts, Dev=anom)


  temp_plot <- forecast::autoplot(ts_dev[,"Temp"]) +
  labs(y = expression(Temperature~(degree*C)), x = "Time") +
  ggtitle("Temperature")

  dev_plot <- forecast::autoplot(ts_dev[,"Dev"]) +
  labs(y = expression(Temperature~(degree*C)), x = "Time") +
  ggtitle("Temperature anomalies")

  plot(dev_plot)
  #return(dev_plot)
  
}


