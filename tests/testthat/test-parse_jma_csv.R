# tests/testthat/test-parse_jma_csv.R

test_that(".parse_jma_csv parses correctly", {
  raw <- charToRaw(
    "Year,Month,Day,areaNo,flag,Temp
2000,1,1,138,0,15.2
2000,1,2,138,0,16.0
END"
  )

  res <- .parse_jma_csv(raw)

  expect_identical(nrow(res), 1L)
  expect_named(res, c("date", "Temp", "flag"))
})


test_that("date column is correct", {
  raw <- charToRaw(
    "Year,Month,Day,areaNo,flag,Temp
2000,1,2,138,0,15
2000,1,3,138,0,16
END"
  )

  res <- .parse_jma_csv(raw)

  expect_identical(res$date[1], as.Date("2000-01-02"))
})


test_that("-999 is converted to NA", {
  raw <- charToRaw(
    "Year,Month,Day,areaNo,flag,Temp
2000,1,1,138,0,-999
2000,1,2,138,0,16
END"
  )

  res <- .parse_jma_csv(raw)

  expect_true(is.na(res$Temp[1]))
})


test_that("last row is removed", {
  raw <- charToRaw(
    "Year,Month,Day,areaNo,flag,Temp
2000,1,1,138,0,15
2000,1,2,138,0,16
END"
  )

  res <- .parse_jma_csv(raw)

  expect_identical(nrow(res), 1L)
})


test_that("fails on too small dataset", {
  raw <- charToRaw(
    "Year,Month,Day,areaNo,flag,Temp
2000,1,1,138,0,15"
  )

  expect_error(
    .parse_jma_csv(raw),
    "unexpected or empty"
  )
})
