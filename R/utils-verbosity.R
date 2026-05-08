#' Get current verbosity level
#'
#' Internal helper to retrieve the current verbosity setting.
#'
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
