data(LondonHP)
m <- NULL

test_that("GW Correlation: pair works", {
  m <<- expect_no_error(gwcorrelation(
    formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
    data = LondonHP
  ))
})

test_that("GW Correlation: works parallel", {
  m <<- expect_no_error(gwcorrelation(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP,
    parallel_method = "omp",
    parallel_arg = 4
  ))
})

test_that("GW Correlation: specific config", {
  #### all by name
  expect_no_error({
    gwcorrelation(
      formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        PURCHASE_UNEMPLOY = gwcorr_config(kernel="bisquare", adaptive = TRUE),
        PURCHASE_PROF = gwcorr_config(kernel="bisquare", adaptive = FALSE),
        FLOORSZ_UNEMPLOY = gwcorr_config(kernel="gaussian", adaptive = TRUE),
        FLOORSZ_PROF = gwcorr_config(kernel="gaussian", adaptive = FALSE)
      )
    )
  })
  #### use default value
  expect_no_error({
    gwcorrelation(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        PURCHASE_FLOORSZ = gwcorr_config(bw = 30, optim_bw = "no", adaptive = TRUE),
        .default = gwcorr_config(adaptive = FALSE)
      )
    )
  })
  #### default
  gwcorrelation(
    formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
    data = LondonHP,
    config = list(gwcorr_config(adaptive = FALSE, kernel = "bisquare"))
    )
})

test_that("GW Correlation: expect error", {
  #### error when default is missing and name is wrong
  expect_error({
    gwcorrelation(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        FLOORSZ = gwcorr_config(bw = 30, optim_bw = "no", adaptive = TRUE)
      )
    )
  }, "Either provide configs for all variable combinations using '_', or supply a default config named '.default'.")

  #### error when default is missing and variable numbers wrong
  expect_error({
    gwcorrelation(
      formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
        # gwcorr_config(adaptive = TRUE, kernel = "gaussian"),
        # gwcorr_config(adaptive = TRUE, kernel = "bisquare"),
        gwcorr_config(adaptive = TRUE, kernel = "bisquare"),
        gwcorr_config(adaptive = TRUE, kernel = "gaussian")
      )
    )
  }, "The length of config mush be equal to the number of variables' product.")

})

test_that("GW Correlation: verbose", {
  skip_on_ci()
  expect_no_error(gwcorrelation(
    formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 1
  ))
  expect_no_error(gwcorrelation(
    formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 2
  ))
  # expect_no_error(gwcorrelation(
  #   formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
  #   data = LondonHP,
  #   verbose = 3
  # ))
})
