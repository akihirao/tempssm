# tests/testthat/setup-tempssm.R

set.seed(123)

# ---- base temperature series ------------------------------------------

temp_ts_test <- ts(
  10 +
    seq(0, 1, length.out = 120) +  # weak trend
    rep(rnorm(12, 0, 0.5), 10) +  # seasonal cycle
    rnorm(120, 0, 0.3),
  start = c(2000, 1),
  frequency = 12
)

# ---- baseline model (no exogenous) ------------------------------------

res_tempssm <- tempssm(temp_ts_test)

# ---- exogenous variables ----------------------------------------------

exo_ts_test <- ts(
  matrix(
    rnorm(120),
    ncol = 1
  ),
  start = c(2000, 1),
  frequency = 12
)


tempssm::set_ts_name(exo_ts_test, label=c("var1"))

# ---- model with exogenous variables -----------------------------------

res_tempssm_exo <- tempssm(
  temp_data = temp_ts_test,
  exo_data = exo_ts_test
)



# ---- small ts dataset -----------------------------------
temp_ts_small <- ts(rnorm(48), frequency = 12)

exo_ts_small <- ts(
    matrix(rnorm(48), ncol = 1),
    frequency = 12
  )
tempssm::set_ts_name(exo_ts_small, label=c("var1"))
