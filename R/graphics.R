#' Plot the estimated level component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated level component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @param res
#' An object returned by \code{tempssm()} from the \pkg{tempssm} package.
#'
#' @param ci
#' Logical; if TRUE (default), pointwise confidence intervals are shown
#' as a shaded ribbon.
#' 
#' @param ylab
#' Character string giving label of y-axis. 
#' Defalut is "temperature".
#'
#' @param ci_level
#' Numeric confidence level between 0 and 1.
#' Defaults to \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @details
#' The confidence interval is computed using
#' \code{stats::confint()} applied to the Kalman filter and smoother
#' results stored in \code{res$kfs}. The shaded ribbon represents
#' pointwise confidence intervals for the level state.
#'
#' @return
#' A \code{ggplot} object, allowing further customization by adding
#' standard \pkg{ggplot2} layers.
#'
#' @seealso
#' \code{\link{tempssm}}, \code{\link[stats]{confint}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' # Default 95% confidence interval
#' autoplot_level(res)
#'
#' # Custom confidence level
#' autoplot_level(res, ci_level = 0.9)
#' }
autoplot_level <- function(res,
                           ci = TRUE,
                           ci_level = 0.95,
                           ylab = "Temperature") {
  
  # ---- Input validation ----
  if (!is.list(res) || is.null(res$kfs)) {
    stop("`res` must be an object returned by tempssm().",
         call. = FALSE)
  }
  
  if (!is.logical(ci) || length(ci) != 1) {
    stop("`ci` must be a single logical value.", call. = FALSE)
  }
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a numeric value between 0 and 1.",
         call. = FALSE)
  }
  
  if (is.null(res$kfs$alphahat)) {
    stop("State estimates not found in `res$kfs$alphahat`.",
         call. = FALSE)
  }
  
  # ---- Extract estimates ----
  alpha_hat <- res$kfs$alphahat
  level <- alpha_hat[, "level"]
  
  level_df <- data.frame(
    time  = time(level),
    level = as.numeric(level)
  )
  
  if (ci) {
    ci_res <- stats::confint(res$kfs, level = ci_level)
    level_df <- cbind(level_df, as.data.frame(ci_res$level))
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }
  
  level_tidy <- tibble::as_tibble(level_df)
  
  # ---- Plot ----
  p <- ggplot2::ggplot(
    level_tidy,
    ggplot2::aes(x = .data$time, y = .data$level)
  ) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::labs(
      title = if (ci) {
        paste0("Level component (", ci_lab, ")")
      } else {
        "Level component"
      },
      x = "Time",
      y = ylab
    )
  
  if (ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
        alpha = 0.3
      )
  }
  
  p
}



#' Plot the estimated drift (slope) component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @param res
#' An object returned by \code{tempssm()} from the \pkg{tempssm} package.
#'
#' @param ci
#' Logical; if TRUE (default), pointwise confidence intervals are shown
#' as a shaded ribbon.
#' 
#' @param ylab
#' Character string giving label of y-axis. 
#' Defalut is "temperature".
#'
#' @param ci_level
#' Numeric confidence level between 0 and 1.
#' Defaults to \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @details
#' The confidence interval is computed using
#' \code{stats::confint()} applied to the Kalman filter and smoother
#' results stored in \code{res$kfs}. The shaded ribbon represents
#' pointwise confidence intervals for the level state.
#'
#' @return
#' A \code{ggplot} object, allowing further customization by adding
#' standard \pkg{ggplot2} layers.
#'
#' @seealso
#' \code{\link{tempssm}}, \code{\link[stats]{confint}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' # Default 95% confidence interval
#' autoplot_drift(res)
#'
#' # Custom confidence level
#' autoplot_drift(res, ci_level = 0.9)
#' }
autoplot_drift <- function(res,
                           ci = TRUE,
                           ci_level = 0.95,
                           ylab = "Temperature") {
  
  # ---- Input validation ----
  if (!is.list(res) || is.null(res$kfs)) {
    stop("`res` must be an object returned by tempssm().",
         call. = FALSE)
  }
  
  if (!is.logical(ci) || length(ci) != 1) {
    stop("`ci` must be a single logical value.", call. = FALSE)
  }
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a numeric value between 0 and 1.",
         call. = FALSE)
  }
  
  if (is.null(res$kfs$alphahat)) {
    stop("State estimates not found in `res$kfs$alphahat`.",
         call. = FALSE)
  }
  
  # ---- Extract estimates ----
  alpha_hat <- res$kfs$alphahat
  slope <- alpha_hat[, "slope"]
  
  drift_df <- data.frame(
    time  = time(slope),
    slope = as.numeric(slope)
  )
  
  if (ci) {
    ci_res <- stats::confint(res$kfs, level = ci_level)
    slope_df <- cbind(drift_df, as.data.frame(ci_res$slope))
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }
  
  slope_tidy <- tibble::as_tibble(slope_df)
  
  # ---- Plot ----
  p <- ggplot2::ggplot(
    slope_tidy,
    ggplot2::aes(x = .data$time, y = .data$slope)
  ) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::labs(
      title = if (ci) {
        paste0("Drift component (", ci_lab, ")")
      } else {
        "Drift component"
      },
      x = "Time",
      y = ylab
    )
  
  if (ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
        alpha = 0.3
      )
  }
  
  p
}



#' Plot the estimated seasonal component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @param res
#' An object returned by \code{tempssm()} from the \pkg{tempssm} package.
#'
#' @param ci
#' Logical; if TRUE (default), pointwise confidence intervals are shown
#' as a shaded ribbon.
#' 
#' @param ylab
#' Character string giving label of y-axis. 
#' Defalut is "temperature".
#'
#' @param ci_level
#' Numeric confidence level between 0 and 1.
#' Defaults to \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @details
#' The confidence interval is computed using
#' \code{stats::confint()} applied to the Kalman filter and smoother
#' results stored in \code{res$kfs}. The shaded ribbon represents
#' pointwise confidence intervals for the level state.
#'
#' @return
#' A \code{ggplot} object, allowing further customization by adding
#' standard \pkg{ggplot2} layers.
#'
#' @seealso
#' \code{\link{tempssm}}, \code{\link[stats]{confint}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' # Default 95% confidence interval
#' autoplot_season(res)
#'
#' # Custom confidence level
#' autoplot_season(res, ci_level = 0.9)
#' }
autoplot_season <- function(res,
                            ci = TRUE,
                            ci_level = 0.95,
                            ylab = "Temperature") {
  
  # ---- Input validation ----
  if (!is.list(res) || is.null(res$kfs)) {
    stop("`res` must be an object returned by tempssm().",
         call. = FALSE)
  }
  
  if (!is.logical(ci) || length(ci) != 1) {
    stop("`ci` must be a single logical value.", call. = FALSE)
  }
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a numeric value between 0 and 1.",
         call. = FALSE)
  }
  
  if (is.null(res$kfs$alphahat)) {
    stop("State estimates not found in `res$kfs$alphahat`.",
         call. = FALSE)
  }
  
  # ---- Extract estimates ----
  alpha_hat <- res$kfs$alphahat
  season <- alpha_hat[, "sea_dummy1"]
  
  season_df <- data.frame(
    time  = time(season),
    season = as.numeric(season)
  )
  
  if (ci) {
    ci_res <- stats::confint(res$kfs, level = ci_level)
    season_df <- cbind(season_df, as.data.frame(ci_res$sea_dummy1))
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }
  
  season_tidy <- tibble::as_tibble(season_df)
  
  # ---- Plot ----
  p <- ggplot2::ggplot(
    season_tidy,
    ggplot2::aes(x = .data$time, y = .data$season)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = if (ci) {
        paste0("Seasonal component (", ci_lab, ")")
      } else {
        "Seasonal component"
      },
      x = "Time",
      y = ylab
    )
  
  if (ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
        alpha = 0.3
      )
  }
  
  p
}



#' Plot the estimated autoregressive component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @param res
#' An object returned by \code{tempssm()} from the \pkg{tempssm} package.
#'
#' @param ci
#' Logical; if TRUE (default), pointwise confidence intervals are shown
#' as a shaded ribbon.
#' 
#' @param ylab
#' Character string giving label of y-axis. 
#' Defalut is "temperature".
#'
#' @param ci_level
#' Numeric confidence level between 0 and 1.
#' Defaults to \code{0.95}, corresponding to a 95\% confidence interval.
#'
#' @details
#' The confidence interval is computed using
#' \code{stats::confint()} applied to the Kalman filter and smoother
#' results stored in \code{res$kfs}. The shaded ribbon represents
#' pointwise confidence intervals for the level state.
#'
#' @return
#' A \code{ggplot} object, allowing further customization by adding
#' standard \pkg{ggplot2} layers.
#'
#' @seealso
#' \code{\link{tempssm}}, \code{\link[stats]{confint}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' # Default 95% confidence interval
#' autoplot_ar1(res)
#'
#' # Custom confidence level
#' autoplot_ar1(res, ci_level = 0.9)
#' }
autoplot_ar1 <- function(res,
                        ci = TRUE,
                        ci_level = 0.95,
                        ylab = "Temperature") {
  
  # ---- Input validation ----
  if (!is.list(res) || is.null(res$kfs)) {
    stop("`res` must be an object returned by tempssm().",
         call. = FALSE)
  }
  
  if (!is.logical(ci) || length(ci) != 1) {
    stop("`ci` must be a single logical value.", call. = FALSE)
  }
  
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a numeric value between 0 and 1.",
         call. = FALSE)
  }
  
  if (is.null(res$kfs$alphahat)) {
    stop("State estimates not found in `res$kfs$alphahat`.",
         call. = FALSE)
  }
  
  # ---- Extract estimates ----
  alpha_hat <- res$kfs$alphahat
  ar1 <- alpha_hat[, "arima1"]
  
  ar1_df <- data.frame(
    time  = time(ar1),
    ar1 = as.numeric(ar1)
  )
  
  if (ci) {
    ci_res <- stats::confint(res$kfs, level = ci_level)
    ar1_df <- cbind(season_df, as.data.frame(ci_res$arima1))
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }
  
  ar1_tidy <- tibble::as_tibble(ar1_df)
  
  # ---- Plot ----
  p <- ggplot2::ggplot(
    ar1_tidy,
    ggplot2::aes(x = .data$time, y = .data$ar1)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = if (ci) {
        paste0("Autoregressive (1) component (", ci_lab, ")")
      } else {
        "Autoregressive (1)"
      },
      x = "Time",
      y = ylab
    )
  
  if (ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
        alpha = 0.3
      )
  }
  
  p
}



#' Function to plot time series of monthly temperature and temperature deviation
#'
#' @import ggplot2
#'
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
