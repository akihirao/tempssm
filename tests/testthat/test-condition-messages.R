extract_literal_condition_messages <- function(lines) {
  pattern <- paste0(
    "(cli::cli_abort|cli::cli_warn|cli::cli_inform|",
    "stop|warning|message)\\(\"([^\"]+)\""
  )

  matches <- regmatches(lines, gregexpr(pattern, lines, perl = TRUE))
  matches <- unlist(matches, use.names = FALSE)
  if (length(matches) == 0L) {
    return(character())
  }

  sub(pattern, "\\2", matches, perl = TRUE)
}


test_that("literal condition messages in R code are unique", {
  r_files <- list.files(
    test_path("../../R"),
    pattern = "[.]R$",
    full.names = TRUE
  )

  source_lines <- unlist(lapply(r_files, readLines, warn = FALSE))
  messages <- extract_literal_condition_messages(source_lines)
  duplicated_messages <- unique(messages[duplicated(messages)])

  expect_identical(duplicated_messages, character())
})


find_condition_expectations <- function(expr, file) {
  calls <- list()

  walk <- function(node) {
    if (!is.call(node)) {
      return(invisible(NULL))
    }

    call_name <- as.character(node[[1L]])[[1L]]
    if (identical(call_name, "function")) {
      return(invisible(NULL))
    }

    if (call_name %in% c("expect_error", "expect_warning", "expect_message")) {
      calls[[length(calls) + 1L]] <<- list(
        file = file,
        call = paste(deparse(node, width.cutoff = 80), collapse = " "),
        n_args = length(as.list(node)) - 1L
      )
    }

    for (idx in seq_along(node)[-1L]) {
      child <- tryCatch(node[[idx]], error = function(e) NULL)
      if (is.null(child)) {
        next
      }
      tryCatch(
        walk(child),
        missingArgError = function(e) NULL
      )
    }
  }

  walk(expr)
  calls
}


test_that("condition expectations compare expected messages", {
  test_files <- list.files(
    test_path("."),
    pattern = "[.]R$",
    full.names = TRUE
  )
  test_files <- test_files[basename(test_files) != "test-condition-messages.R"]

  condition_calls <- unlist(
    lapply(test_files, function(file) {
      expressions <- parse(file, keep.source = FALSE)
      unlist(
        lapply(as.list(expressions), find_condition_expectations, file = file),
        recursive = FALSE
      )
    }),
    recursive = FALSE
  )

  missing_expected <- vapply(
    condition_calls,
    function(call) {
      call$n_args < 2L
    },
    logical(1)
  )

  missing_calls <- vapply(
    condition_calls[missing_expected],
    function(call) {
      paste(basename(call$file), call$call, sep = ": ")
    },
    character(1)
  )

  expect_identical(missing_calls, character())
})
