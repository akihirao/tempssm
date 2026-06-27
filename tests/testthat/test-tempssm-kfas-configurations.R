test_that("KFAS model configurations retain their state-name mappings", {
  y <- ts(rnorm(16), frequency = 4)
  exo_matrix <- cbind(
    pdo = rnorm(16),
    enso = rnorm(16)
  )
  configurations <- list(
    list(use_season = FALSE, exo_matrix = NULL, ar_order = 1),
    list(use_season = TRUE, exo_matrix = NULL, ar_order = 1),
    list(use_season = FALSE, exo_matrix = exo_matrix, ar_order = 1),
    list(use_season = TRUE, exo_matrix = exo_matrix, ar_order = 2)
  )

  for (config in configurations) {
    param_idx <- .get_param_index(
      ar_order = config$ar_order,
      use_season = config$use_season
    )
    initial_model <- .define_build_model(
      y = y,
      freq = 4,
      use_season = config$use_season,
      exo_mat = config$exo_matrix,
      ar_order = config$ar_order
    )
    update_model <- .define_update_func(
      y = y,
      freq = 4,
      use_season = config$use_season,
      exo_mat = config$exo_matrix,
      ar_order = config$ar_order,
      ar_idx = param_idx$ar,
      var_idx = param_idx$var,
      H_idx = param_idx$H
    )
    updated_model <- update_model(
      .default_tempssm_inits(config$ar_order, config$use_season),
      initial_model
    )

    exo_names <- colnames(config$exo_matrix)
    expected_exo_names <- if (is.null(exo_names)) character(0) else exo_names
    n_states <- dim(updated_model$T)[1]
    state_names <- .make_tempssm_state_names(
      exo_names = exo_names,
      freq = 4,
      use_season = config$use_season,
      ar_order = config$ar_order,
      n_states = n_states
    )

    expect_s3_class(initial_model, "SSModel")
    expect_s3_class(updated_model, "SSModel")
    expect_identical(dim(initial_model$T), dim(updated_model$T))
    expect_length(state_names, n_states)
    expect_identical(
      head(state_names, length(expected_exo_names)),
      expected_exo_names
    )
    expect_true(all(c("level", "slope") %in% state_names))
  }
})
