
#' Plot estimated components of a ThermoSSM object
#'
#' @param x An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param components Components to plot. Any of \code{"level"}, \code{"drift"},
#'   \code{"seasonal"}, \code{"ar"}.
#' @param ci Logical; if TRUE, confidence intervals are shown for
#'   level and drift components.
#' @details
#' When \code{ci = TRUE}, shaded bands represent pointwise confidence intervals
#' for the smoothed state estimates of the level and drift (slope) components.
#'
#' These intervals are computed from the smoothed state covariance matrix
#' obtained by Kalman smoothing under the assumption of a linear Gaussian
#' state-space model. Specifically, confidence intervals are given by
#' the estimated state mean plus or minus a normal quantile multiplied by
#' the square root of the corresponding state variance.
#'
#' Note that these confidence intervals quantify uncertainty in the latent
#' state estimates conditional on the observed data, and should not be
#' interpreted as prediction intervals for future observations.
#' 
#' @param ci_level Confidence level for intervals (default 0.95).
#' @param layout Layout of plots: \code{"grid"} or \code{"list"}.
#' @param ... Further arguments (currently unused).
#'
#' @return A ggplot object or a patchwork object.
#'
#' @export
plot.ThermoSSM <- function(
  x,
  components = c("level", "drift", "seasonal", "ar"),
  ci = FALSE,
  ci_level = 0.95,
  layout = c("grid", "list"),
  ...
) {

  layout <- match.arg(layout)

  if (!inherits(x, "ThermoSSM")) {
    stop("Object must be of class 'ThermoSSM'.", call. = FALSE)
  }

  freq <- frequency(x$data_temp)
  if(freq==12){
    drift_plot_y_lab <- "Temperature (\u00B0C/month)"
  }else{
    drift_plot_y_lab <- "Temperature (\u00B0C/unit)"
  }

  alpha_hat <- x$kfs$alphahat
  time_index <- time(alpha_hat)

  if (ci) {
    ci_obj <- confint(x$kfs, level = ci_level)
  }

  plots <- list()

  ## ---- level ----
  if ("level" %in% components) {

    df <- tibble::tibble(
      time  = time_index,
      value = as.numeric(alpha_hat[, "level"])
    )

    if (ci) {
      df$lwr <- ci_obj$level[,"lwr"]
      df$upr <- ci_obj$level[,"upr"]
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::labs(title = "Level component", x = "Time", y = "Temperature (\u00B0C)")

    if (ci) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lwr, ymax = upr),
          alpha = 0.3#, fill = "grey"
        )
    }

    plots$level <- p
  }

  ## ---- drift (slope) ----
  if ("drift" %in% components) {

    df <- tibble::tibble(
      time  = time_index,
      value = as.numeric(alpha_hat[, "slope"])
    )

    if (ci) {
      df$lwr <- ci_obj$slope[,"lwr"]
      df$upr <- ci_obj$slope[,"upr"]
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::labs(title = "Drift (slope) component", x = "Time", y = drift_plot_y_lab)

    if (ci) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lwr, ymax = upr),
          alpha = 0.3#, fill = "grey70"
        )
    }

    plots$drift <- p
  }

  ## ---- seasonal ----
  if ("seasonal" %in% components) {
    plots$seasonal <-
      forecast::autoplot(alpha_hat[, "sea_dummy1"]) +
      ggplot2::labs(title = "Seasonal component", x = "Time", y = "Temperature (\u00B0C)")
  }

  ## ---- AR ----
  if ("ar" %in% components) {
    plots$ar <-
      forecast::autoplot(alpha_hat[, "arima1"]) +
      ggplot2::labs(title = "Autoregressive component", x = "Time", y = "Temperature (\u00B0C)")
  }

  if (layout == "list") {
    return(plots)
  }

  patchwork::wrap_plots(plots)
}
