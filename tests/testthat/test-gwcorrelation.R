data(LondonHP)
m <- NULL

test_that("GW Correlation: works", {
  m <<- expect_no_error(gw_correlation(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP
  ))
})

test_that("GW Correlation: pair works", {
  m <<- expect_no_error(gw_correlation(
    formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
    data = LondonHP
  ))
})

test_that("GW Correlation: specific config by names", {
  #### all by name
  expect_no_error({
    gw_correlation(
      formula = PURCHASE + FLOORSZ ~ UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        PURCHASE_UNEMPLOY = mgwr_config(adaptive = TRUE),
        PURCHASE_PROF = mgwr_config(adaptive = TRUE),
        FLOORSZ_UNEMPLOY = mgwr_config(kernel="gaussian", adaptive = TRUE),
        FLOORSZ_PROF = mgwr_config(kernel="gaussian", adaptive = TRUE)
      )
    )
  })
  #### use default value
  expect_no_error({
    gw_correlation(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        PURCHASE_FLOORSZ = mgwr_config(bw = 30, optim_bw = "no", adaptive = TRUE),
        .default = mgwr_config(adaptive = TRUE)
      )
    )
  })
  #### error when default is missing and config is missing
  expect_error({
    gw_correlation(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        FLOORSZ = mgwr_config(bw = 30, optim_bw = "no", adaptive = TRUE)
      )
    )
  }, "Either provide configs for all variable combinations using '_', or supply a default config named '.default'.")
})

# test_that("GW Correlation: verbose", {
#   skip_on_ci()
#   expect_no_error(gw_correlation(
#     formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#     data = LondonHP,
#     verbose = 1
#   ))
#   expect_no_error(gw_correlation(
#     formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#     data = LondonHP,
#     verbose = 2
#   ))
#   expect_no_error(gw_correlation(
#     formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#     data = LondonHP,
#     verbose = 3
#   ))
#   expect_no_error(gw_correlation(
#     formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
#     data = LondonHP,
#     verbose = 4
#   ))
# })
