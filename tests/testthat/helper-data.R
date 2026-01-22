# helper-data.R

set.seed(123)

temp_ts_test <- ts(
  10 +
    seq(0, 1, length.out = 120) +      # weak trend
    rep(rnorm(12, 0, 0.5), 10) +       # seasonal cycle
    rnorm(120, 0, 0.3),                # noise
  start = c(2000, 1),
  frequency = 12
)
