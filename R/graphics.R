
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
    labs(title="Seasonal trend of monthly average temperature",
      y = expression(Temperature~(degree*C))) +
    scale_x_discrete(
      labels = function(x) sprintf("%02d", as.integer(x))
    )


  plot(plot_monthly_ave_temp )

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

  
}


#' Function to plot estimated components: level, drift, seasonal, and auto-regression
#' @import ggplot2
#' @param res output of model with using lgssm()
#'
#' @encoding UTF-8
#'
#' @export
#'
plot_level_trend_season_ar <- function(res){

  alpha_hat <- res[[2]]$alphahat

  model_level_plot <- forecast::autoplot(alpha_hat[,"level"]) + 
    labs(y = "", x = "Time") + 
    ggtitle("Level component")

  model_slope_plot <- forecast::autoplot(alpha_hat[,"slope"]) +
    labs(y = "", x = "Time") +
    ggtitle("Drift component")

  model_season_plot <- forecast::autoplot(alpha_hat[,"sea_dummy1"]) +
    labs(y = "", x = "Time") +
    ggtitle("Seasonal component")

  model_arima1_plot <- forecast::autoplot(alpha_hat[,"arima1"]) +
    labs(y = "", x = "Time") +
    ggtitle("Auto-regression component")

   level_drift_seaon_ar_plots <- list(
    model_level_plot,
    model_slope_plot,
    model_season_plot,
    model_arima1_plot
    )

   return(level_drift_seaon_ar_plots)

}



#' Function to plot level component with confidence interval
#' @import ggplot2
#' @param res output of model with using lgssm()
#'
#' @encoding UTF-8
#'
#' @export
#'
plot_level_ci <- function(res,ci_range=0.95){

  smooth <- res[[2]] # smoothing
  alpha_hat <- smooth$alphahat

  ci_lab <- paste0(round(ci_range*100,0),"% CI")
  level <- alpha_hat[,"level"]
  ci <- confint(res[[2]], level = ci_range)

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
#' @param res output of model with using lgssm()
#'
#' @encoding UTF-8
#'
#' @export
#'
plot_drift_ci <- function(res,ci_range=0.95){

  ci <- confint(res[[2]], level = ci_range)
  smooth <- res[[2]] # smoothing
  alpha_hat <- smooth$alphahat

  drift <- alpha_hat[,"slope"]

  mean_drift_year <- mean(drift) * 12
  

  drift_tidy <- cbind(
    data.frame(time=time(drift),
               level=drift),
  as.data.frame(ci$slope)
  ) %>%
  as_tibble()

  ci_lab <- paste0("Drift component (grey area: ", round(ci_range*100,0),"% CI)")
  sub_lab <- paste0("Average drift rate per year = " ,mean_drift_year)

  drift_plot <- ggplot(data=drift_tidy,
                         aes(x=time,y=drift)) +
    labs(title=ci_lab,
      subtitle= sub_lab,
       x="Year", y="Drift") +
    geom_line(aes(y=level), size = 1.2) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_hline(yintercept=0, linetype="dashed") 

  return(drift_plot)

}
