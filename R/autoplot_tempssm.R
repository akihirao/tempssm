
#' Arrange ggplot objects and return a gtable
#'
#' @param plots A list of ggplot objects
#' @param nrow Number of rows
#' @param ncol Number of columns
#'
#' @return A gtable object (grob)
#' @keywords internal
#' @noRd
.make_gtable_layout <- function(plots,
                               nrow = length(plots),
                               ncol = 1) {
  if (!is.list(plots)) {
    stop("`plots` must be a list.", call. = FALSE)
  }

  grobs <- lapply(plots, ggplot2::ggplotGrob)

  mat <- matrix(grobs,
                nrow = nrow,
                ncol = ncol,
                byrow = TRUE)

  gtable::gtable_matrix(
    name = "layout",
    grobs = mat,
    widths  = grid::unit(rep(1, ncol), "null"),
    heights = grid::unit(rep(1, nrow), "null")
  )
}


#' Autoplot method for tempssm objects
#'
#' @description
#' Create diagnostic and component plots for a model fitted by
#' \code{tempssm()}. By default, all available components are plotted.
#'
#' @param object
#' An object returned by \code{tempssm()}.
#'
#' @param component
#' Character string specifying which component to plot.
#' One of \code{"level"}, \code{"drift"}, \code{"season"}, or \code{"ar1"}.
#' If NULL (default), all components are plotted.
#'
#' @param ci
#' Logical; if TRUE (default), pointwise confidence intervals are shown.
#'
#' @param ci_level
#' Numeric confidence level between 0 and 1.
#' Defaults to \code{0.95}.
#'
#' @param ...
#' Additional arguments passed to the corresponding
#' \code{autoplot_*()} function.
#'
#' @return
#' A \code{ggplot} object if a single component is requested,
#' or a \code{gtable} object combining all component plots.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' # plot all components with 95% confidence interval
#' autoplot(res)
#'
#' # plot all components without 95% confidence interval
#' autoplot(res, ci = FALSE)
#'
#' # plot each of components
#' autoplot(res, component = "level", ci = FALSE)
#' }
#'
#' @method autoplot tempssm
#' @export
autoplot.tempssm <- function(object,
                             component = NULL,
                             ci = TRUE,
                             ci_level = 0.95,
                             ...) {

  plotters <- list(
    level  = autoplot_level,
    drift  = autoplot_drift,
    season = autoplot_season,
    ar1    = autoplot_ar1
  )

  if (!is.null(component)) {
    if (!is.character(component) || length(component) != 1) {
      stop("`component` must be a single character string.",
           call. = FALSE)
    }

    if (!component %in% names(plotters)) {
      stop(
        "`component` must be one of: ",
        paste0(names(plotters), collapse = ", "),
        call. = FALSE
      )
    }

    return(
      plotters[[component]](
        object,
        ci = ci,
        ci_level = ci_level,
        ...
      )
    )
  }

  plots <- lapply(
    plotters,
    function(f) f(object, ci = ci, ci_level = ci_level, ...)
  )

  g <- .make_gtable_layout(plots, ncol = 1)

  grid::grid.newpage()
  grid::grid.draw(g)

  invisible(g)
}
