#' @keywords internal
.get_verbosity <- function() {
  v <- Sys.getenv("TEMPSSM_VERBOSITY", "")

  if (v == "") {
    v <- getOption("tempssm.verbosity", "inform")
  }

  v <- as.character(v)

  allowed <- c("none", "inform", "debug")
  if (!v %in% allowed) {
    v <- "inform"
  }

  v
}

#' @keywords internal
.tempssm_cli_inform <- function(..., .envir = parent.frame()) {
  if (.get_verbosity() %in% c("inform", "debug")) {
    cli::cli_inform(..., .envir = .envir)
  }
}

#' @keywords internal
.tempssm_cli_debug <- function(..., .envir = parent.frame()) {
  if (.get_verbosity() == "debug") {
    cli::cli_inform(..., .envir = .envir)
  }
}
