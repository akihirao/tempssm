#' Plot the estimated level component from a tempssm model
#'
#' @importFrom rlang .data
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated level component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @inheritParams get_level_ts
#'
#' @param ylab
#' Character string giving label of y-axis.
#' Defalut is "temperature".
#'
#' @param show_ci_in_title Logical; should the confidence level be shown in
#'   the plot title when \code{ci = TRUE}? The default is \code{FALSE}.
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
                           ylab = "Temperature",
                           show_ci_in_title = FALSE) {
  ## ---- input checks ---------------------------------------------------
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  if (!is.logical(ci) || length(ci) != 1) {
    cli::cli_abort("`ci` must be a single logical value.")
  }

  if (!is.logical(show_ci_in_title) || length(show_ci_in_title) != 1) {
    cli::cli_abort("`show_ci_in_title` must be a single logical value.")
  }

  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
      cli::cli_abort(
        "`ci_level` must be a numeric value between 0 and 1."
      )
    }
  }

  .tempssm_cli_debug(
    "Creating level plot (ci={ci}, ci_level={ci_level})"
  )

  ## ---- extract data ---------------------------------------------------
  level_ts <- get_level_ts(res, ci = TRUE, ci_level = ci_level)

  level_df <- data.frame(
    time  = time(level_ts),
    level = as.numeric(level_ts[, "level"])
  )

  if (ci) {
    .tempssm_cli_debug("Including confidence intervals in plot")

    ci_df <- data.frame(
      lwr = as.numeric(level_ts[, "lwr"]),
      upr = as.numeric(level_ts[, "upr"])
    )

    level_df <- cbind(level_df, ci_df)
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }

  level_tidy <- tibble::as_tibble(level_df)

  ## ---- plotting -------------------------------------------------------
  p <- ggplot2::ggplot(
    level_tidy,
    ggplot2::aes(x = .data$time, y = .data$level)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::labs(
      title = if (ci && show_ci_in_title) {
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

  .tempssm_cli_debug("Level plot created successfully")

  return(p)
}


#' Plot the estimated drift (slope) component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @inheritParams get_level_ts
#' @inheritParams autoplot_level
#'
#' @param ylab
#' Character string giving label of y-axis.
#' Defalut is "temperature".
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
                           ylab = "Temperature",
                           show_ci_in_title = FALSE) {
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  # ---- check ci first ----
  if (!is.logical(ci) || length(ci) != 1) {
    cli::cli_abort("`ci` must be a single logical value.")
  }

  if (!is.logical(show_ci_in_title) || length(show_ci_in_title) != 1) {
    cli::cli_abort("`show_ci_in_title` must be a single logical value.")
  }

  # ---- check ci_level only if ci is TRUE ----
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
      cli::cli_abort(
        "`ci_level` must be a numeric value between 0 and 1."
      )
    }
  }

  .tempssm_cli_debug(
    "Creating drift plot (ci={ci}, ci_level={ci_level})"
  )

  ## ---- extract data ---------------------------------------------------
  drift_ts <- get_drift_ts(res, ci = TRUE, ci_level = ci_level)

  drift_df <- data.frame(
    time  = time(drift_ts),
    drift = as.numeric(drift_ts[, "drift"])
  )

  if (ci) {
    .tempssm_cli_debug("Including confidence intervals in plot")

    ci_df <- data.frame(
      lwr = as.numeric(drift_ts[, "lwr"]),
      upr = as.numeric(drift_ts[, "upr"])
    )
    drift_df <- cbind(drift_df, ci_df)
    ci_lab <- paste0(round(ci_level * 100), "% CI")
  }

  drift_tidy <- tibble::as_tibble(drift_df)

  # ---- Plot ----
  p <- ggplot2::ggplot(
    drift_tidy,
    ggplot2::aes(x = .data$time, y = .data$drift)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::labs(
      title = if (ci && show_ci_in_title) {
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

  .tempssm_cli_debug("Drift plot created successfully")

  return(p)
}


#' Plot the estimated seasonal component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @inheritParams get_level_ts
#' @inheritParams autoplot_level
#'
#' @param ylab
#' Character string giving label of y-axis.
#' Defalut is "temperature".
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
                            ylab = "Temperature",
                            show_ci_in_title = FALSE) {
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  # ---- check ci first ----
  if (!is.logical(ci) || length(ci) != 1) {
    cli::cli_abort("`ci` must be a single logical value.")
  }

  if (!is.logical(show_ci_in_title) || length(show_ci_in_title) != 1) {
    cli::cli_abort("`show_ci_in_title` must be a single logical value.")
  }

  # ---- check ci_level only if ci is TRUE ----
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
      cli::cli_abort(
        "`ci_level` must be a numeric value between 0 and 1."
      )
    }
  }

  .tempssm_cli_debug(
    "Creating season plot (ci={ci}, ci_level={ci_level})"
  )

  ## ---- extract data ---------------------------------------------------
  season_ts <- get_season_ts(res, ci = TRUE, ci_level = ci_level)

  season_df <- data.frame(
    time = time(season_ts),
    season = as.numeric(season_ts[, "season"])
  )

  if (ci) {
    .tempssm_cli_debug("Including confidence intervals in plot")

    ci_df <- data.frame(
      lwr = as.numeric(season_ts[, "lwr"]),
      upr = as.numeric(season_ts[, "upr"])
    )
    season_df <- cbind(season_df, ci_df)
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
      title = if (ci && show_ci_in_title) {
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

  .tempssm_cli_debug("Seasonal plot created successfully")

  return(p)
}


#' Plot the estimated autoregressive component from a tempssm model
#'
#' @description
#' Create a \pkg{ggplot2} visualization of the estimated drift component
#' obtained from a state space model fitted by \code{tempssm()}.
#' A pointwise confidence interval is shown as a shaded ribbon.
#'
#' @inheritParams get_level_ts
#' @inheritParams autoplot_level
#'
#' @param ylab
#' Character string giving label of y-axis.
#' Defalut is "temperature".
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
                         ylab = "Temperature",
                         show_ci_in_title = FALSE) {
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  # ---- check ci first ----
  if (!is.logical(ci) || length(ci) != 1) {
    cli::cli_abort("`ci` must be a single logical value.")
  }

  if (!is.logical(show_ci_in_title) || length(show_ci_in_title) != 1) {
    cli::cli_abort("`show_ci_in_title` must be a single logical value.")
  }

  # ---- check ci_level only if ci is TRUE ----
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
      cli::cli_abort(
        "`ci_level` must be a numeric value between 0 and 1."
      )
    }
  }

  .tempssm_cli_debug(
    "Creating ar1 plot (ci={ci}, ci_level={ci_level})"
  )

  ## ---- extract data ---------------------------------------------------
  ar1_ts <- get_ar1_ts(res, ci = TRUE, ci_level = ci_level)

  ar1_df <- data.frame(
    time = time(ar1_ts),
    ar1 = as.numeric(ar1_ts[, "ar1"])
  )

  if (ci) {
    .tempssm_cli_debug("Including confidence intervals in plot")

    ci_df <- data.frame(
      lwr = as.numeric(ar1_ts[, "lwr"]),
      upr = as.numeric(ar1_ts[, "upr"])
    )
    ar1_df <- cbind(ar1_df, ci_df)
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
      title = if (ci && show_ci_in_title) {
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

  .tempssm_cli_debug("AR1 plot created successfully")

  return(p)
}


#' Plot temperature anomalies
#'
#' @description
#' Plots a temperature time series together with its corresponding
#' temperature anomalies.
#'
#' The anomalies are computed by subtracting the long-term seasonal mean
#' for each period in the seasonal cycle from the observed temperature.
#'
#' @param ts
#' A univariate time series object of class \code{ts} with an integer
#' frequency greater than 1.
#'
#' @details
#' The function first computes the seasonal climatological mean across all
#' years, and then calculates anomalies as deviations from these seasonal
#' averages. For monthly data, this is equivalent to subtracting the
#' long-term mean for each calendar month.
#'
#' @return
#' A \code{ggplot2} plot object showing the temperature anomaly time series.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' p <- plot_temp_dev(niigata_sst)
#' print(p)
#' }
#'
#' @importFrom forecast autoplot
#' @importFrom ggplot2 labs ggtitle
#' @export
plot_temp_dev <- function(ts) {
  if (!inherits(ts, "ts")) {
    cli::cli_abort("`ts` must be an object of class {.cls ts}.")
  }

  anom <- tempssm::compute_temp_anomaly(ts)

  dev_plot <- forecast::autoplot(anom) +
    ggplot2::labs(
      y = expression(Temperature ~ (degree * C)),
      x = "Time"
    ) +
    ggplot2::ggtitle("Temperature anomalies")

  return(dev_plot)
}
