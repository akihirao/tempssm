#' Validate models supplied for AIC comparison
#'
#' @param models A named list expected to contain fitted \code{tempssm}
#'   objects.
#'
#' @return \code{models}, unchanged.
#'
#' @keywords internal
#' @noRd
.validate_aic_model_list <- function(models) {
  if (!is.list(models) || length(models) < 2L) {
    cli::cli_abort(
      "{.arg models} must be a named list containing at least two models."
    )
  }

  model_names <- names(models)
  valid_names <- !is.null(model_names) &&
    length(model_names) == length(models) &&
    !anyNA(model_names) &&
    all(nzchar(trimws(model_names)))
  if (!valid_names) {
    cli::cli_abort(
      "Every model in {.arg models} must have a non-empty name."
    )
  }

  if (anyDuplicated(model_names)) {
    cli::cli_abort(
      "Names in {.arg models} must be unique."
    )
  }

  for (i in seq_along(models)) {
    tryCatch(
      .validate_tempssm_for_ic(models[[i]]),
      error = function(e) {
        cli::cli_abort(
          paste0(
            "Model {.val {model_names[[i]]}} cannot be compared: ",
            conditionMessage(e)
          )
        )
      }
    )

    if (!inherits(models[[i]]$temp_data, "ts")) {
      cli::cli_abort(
        paste0(
          "Model {.val {model_names[[i]]}} does not contain valid ",
          "{.cls ts} response data."
        )
      )
    }
  }

  models
}


#' Check whether two models use the same response period
#'
#' @param reference Reference \code{tempssm} object.
#' @param candidate Candidate \code{tempssm} object.
#'
#' @return A logical scalar.
#'
#' @keywords internal
#' @noRd
.same_aic_response_period <- function(reference, candidate) {
  reference_data <- reference$temp_data
  candidate_data <- candidate$temp_data

  identical(NROW(reference_data), NROW(candidate_data)) &&
    identical(
      stats::frequency(reference_data),
      stats::frequency(candidate_data)
    ) &&
    identical(
      as.numeric(stats::time(reference_data)),
      as.numeric(stats::time(candidate_data))
    )
}


#' Check whether two models use identical response observations
#'
#' @inheritParams .same_aic_response_period
#'
#' @return A logical scalar.
#'
#' @keywords internal
#' @noRd
.same_aic_response_values <- function(reference, candidate) {
  reference_values <- as.numeric(reference$temp_data)
  candidate_values <- as.numeric(candidate$temp_data)

  identical(is.na(reference_values), is.na(candidate_values)) &&
    identical(reference_values, candidate_values)
}


#' Validate that fitted models can be compared by AIC
#'
#' @inheritParams .validate_aic_model_list
#'
#' @return A logical vector containing the likelihood setting for each model.
#'
#' @keywords internal
#' @noRd
.validate_aic_comparison_data <- function(models) {
  reference <- models[[1L]]
  model_names <- names(models)

  for (i in seq.int(2L, length(models))) {
    if (!.same_aic_response_period(reference, models[[i]])) {
      cli::cli_abort(
        paste0(
          "Model {.val {model_names[[i]]}} does not use the same response ",
          "period and time index as model {.val {model_names[[1L]]}}."
        )
      )
    }

    if (!.same_aic_response_values(reference, models[[i]])) {
      cli::cli_abort(
        paste0(
          "Model {.val {model_names[[i]]}} does not use the same response ",
          "values and missing-value pattern as model ",
          "{.val {model_names[[1L]]}}."
        )
      )
    }
  }

  marginal <- vapply(
    models,
    .resolve_tempssm_marginal,
    logical(1),
    marginal = NULL
  )
  if (length(unique(marginal)) != 1L) {
    cli::cli_abort(
      paste0(
        "Models fitted with diffuse and marginal likelihoods cannot be ",
        "compared in the same AIC table. Refit all models using the same ",
        "{.arg marginal} setting."
      )
    )
  }

  marginal
}


#' Build an AIC comparison table
#'
#' @inheritParams .validate_aic_model_list
#' @param marginal Logical vector containing model likelihood settings.
#'
#' @return A \code{tibble} ordered by increasing AIC.
#'
#' @keywords internal
#' @noRd
.build_aic_comparison_table <- function(models, marginal) {
  rows <- lapply(seq_along(models), function(i) {
    fitted_model <- models[[i]]
    ll <- logLik(fitted_model)
    aic <- -2 * as.numeric(ll) + 2 * attr(ll, "df")
    if (!is.finite(aic)) {
      cli::cli_abort(
        "AIC for model {.val {names(models)[[i]]}} is not finite."
      )
    }

    tibble::tibble(
      model = names(models)[[i]],
      logLik = as.numeric(ll),
      df = attr(ll, "df"),
      nobs = attr(ll, "nobs"),
      observed_n = sum(!is.na(fitted_model$temp_data)),
      start = paste(stats::start(fitted_model$temp_data), collapse = "-"),
      end = paste(stats::end(fitted_model$temp_data), collapse = "-"),
      frequency = stats::frequency(fitted_model$temp_data),
      likelihood = if (marginal[[i]]) "marginal" else "diffuse",
      AIC = aic
    )
  })

  result <- tibble::as_tibble(do.call(rbind, rows))
  result <- result[order(result$AIC, seq_len(nrow(result))), , drop = FALSE]
  result$delta_AIC <- result$AIC - min(result$AIC)
  relative_likelihood <- exp(-0.5 * result$delta_AIC)
  result$weight <- relative_likelihood / sum(relative_likelihood)
  rownames(result) <- NULL
  result
}


#' Compare AIC across fitted tempssm models
#'
#' @param models A named list containing at least two fitted and converged
#'   \code{tempssm} objects. Names must be non-empty and unique and are used as
#'   model labels in the returned table.
#'
#' @details
#' AIC comparisons require all models to use the same response observations.
#' This function verifies response length, frequency, complete time index,
#' values, and missing-value positions before calculating the table. Model
#' structure may otherwise differ, including autoregressive order, inclusion
#' of a seasonal component, and inclusion or number of exogenous variables.
#'
#' All models must also have been fitted using the same likelihood type. The
#' function uses the \code{marginal} setting stored in each fitted model and
#' does not provide an argument for changing the likelihood after fitting.
#' Models fitted with diffuse and marginal likelihoods are therefore rejected
#' when supplied together.
#'
#' In the current \pkg{tempssm} model, changing \code{use_season} changes the
#' diffuse state structure, and adding exogenous variables adds diffuse
#' regression states. Changing only the AR order adds non-diffuse AR states and
#' does not ordinarily change the diffuse rank. When comparing models whose
#' seasonal or exogenous structure differs, fitting every model with
#' \code{marginal = TRUE} is recommended because the marginal likelihood
#' includes the correction associated with diffuse initialization. Comparisons
#' using \code{marginal = FALSE} remain available when that setting is shared
#' by every model.
#'
#' Akaike weights are calculated as normalized relative likelihoods from
#' \code{delta_AIC}.
#'
#' @return A \code{tibble} ordered by increasing AIC with columns
#'   \code{model}, \code{logLik}, \code{df}, \code{nobs},
#'   \code{observed_n}, \code{start}, \code{end}, \code{frequency},
#'   \code{likelihood}, \code{AIC}, \code{delta_AIC}, and \code{weight}.
#'
#' @seealso \code{\link{tempssm}}, \code{\link{AIC.tempssm}},
#'   \code{\link{get_aic}}
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#'
#' models <- list(
#'   seasonal = tempssm(
#'     niigata_sst,
#'     use_season = TRUE,
#'     marginal = TRUE
#'   ),
#'   no_season = tempssm(
#'     niigata_sst,
#'     use_season = FALSE,
#'     marginal = TRUE
#'   )
#' )
#'
#' compare_tempssm_aic(models)
#' }
#'
#' @export
compare_tempssm_aic <- function(models) {
  models <- .validate_aic_model_list(models)
  marginal <- .validate_aic_comparison_data(models)
  .build_aic_comparison_table(models, marginal)
}
