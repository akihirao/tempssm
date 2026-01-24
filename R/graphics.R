
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

  
}




#' Function to plot level component with confidence interval
#' @import ggplot2
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param ci_level Confidence level for the interval estimation.
#'   Must be a numeric value between 0 and 1.
#'   The default is \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @encoding UTF-8
#'
#' @export
#'
plot_level_ci <- function(res, ci_level=0.95){

  alpha_hat <- res$kfs$alphahat

  ci_lab <- paste0(round(ci_level*100,0),"% CI")
  level <- alpha_hat[,"level"]
  ci <- confint(res$kfs, level = ci_level)

  level_tidy <- cbind(
    data.frame(time=time(level),
               level=level),
  as.data.frame(ci$level)
  ) %>%
  as_tibble()

  level_plot <- ggplot(data=level_tidy,
                         aes(x=time,y=level)) +
    labs(title=paste0("Level component (grey area: ",ci_lab,")"),
       x="Year", y="Temperature") +
    geom_line(aes(y=level), size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)

  return(level_plot)

}





#' Function to plot drift component with confidence interval
#' @import ggplot2
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param ci_level Confidence level for the interval estimation.
#'   Must be a numeric value between 0 and 1.
#'   The default is \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @encoding UTF-8
#'
#' @export
#'
plot_drift_ci <- function(res, ci_level=0.95){

  alpha_hat <- res$kfs$alphahat

  drift <- alpha_hat[,"slope"]
  ci <- confint(res$kfs, level = ci_level)

  mean_drift_year <- round(mean(drift) * 12,4)

  drift_tidy <- cbind(
    data.frame(time=time(drift),
               drift=drift),
  as.data.frame(ci$slope)
  ) %>%
  as_tibble()

  ci_lab <- paste0("Drift component (grey area: ", round(ci_level*100,0),"% CI)")
  sub_lab <- paste0("Average drift rate per year = " ,mean_drift_year)

  drift_plot <- ggplot(data=drift_tidy,
                         aes(x=time,y=drift)) +
    labs(title=ci_lab,
      subtitle= sub_lab,
       x="Year", y="Drift") +
    geom_line(aes(y=drift), size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_hline(yintercept=0, linetype="dashed") 

  return(drift_plot)

}
