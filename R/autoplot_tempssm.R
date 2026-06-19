
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

  if (!is.numeric(nrow) ||
      length(nrow) != 1 ||
      !is.finite(nrow) ||
      nrow < 1 ||
      nrow != as.integer(nrow)) {
    stop("`nrow` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(ncol) ||
      length(ncol) != 1 ||
      !is.finite(ncol) ||
      ncol < 1 ||
      ncol != as.integer(ncol)) {
    stop("`ncol` must be a positive integer.", call. = FALSE)
  }

  if (nrow * ncol < length(plots)) {
    stop("`nrow * ncol` must be at least the number of plots.", call. = FALSE)
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
#' @inheritParams get_level_ts
#'
#' @param nrow,ncol Positive integers specifying the layout used when all
#'   components are plotted. The default is a 2 by 2 layout. These arguments
#'   are ignored when \code{component} is specified.
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
#' # plot all components in one column
#' autoplot(res, nrow = 4, ncol = 1)
#'
#' # plot each of components
#' autoplot(res, component = "level", ci = FALSE)
#' }
#'
#' @method autoplot tempssm
#' @importFrom ggplot2 autoplot
#' @export
autoplot.tempssm <- function(object,
                             component = NULL,
                             ci = TRUE,
                             ci_level = 0.95,
                             nrow = 2,
                             ncol = 2,
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
        toString(names(plotters)),
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

  g <- .make_gtable_layout(plots, nrow = nrow, ncol = ncol)

  grid::grid.newpage()
  grid::grid.draw(g)

  invisible(g)
}


#' Plot method for tempssm objects
#'
#' @description
#' Default base R plot method for objects fitted by \code{tempssm()}.
#' This method delegates to \code{autoplot.tempssm()} so that \code{plot(res)}
#' displays the same component plots as \code{autoplot(res)}.
#'
#' @param x An object returned by \code{tempssm()}.
#' @param ... Additional arguments passed to \code{autoplot.tempssm()}.
#'
#' @return
#' Invisibly returns the plotted \code{ggplot} object for a single component,
#' or the \code{gtable} object combining all component plots.
#'
#' @method plot tempssm
#' @export
plot.tempssm <- function(x, ...) {
  p <- autoplot.tempssm(x, ...)

  if (inherits(p, "ggplot")) {
    print(p)
  }

  invisible(p)
}
