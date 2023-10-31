data(LondonHP)
m <- NULL

test_that("Multiscale GWR: works", {
  m <<- expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP
  ))
})

test_that("Multiscale GWR: specific config by names", {
  #### all by name
  expect_no_error({
    gwr_multiscale(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        Intercept = mgwr_config(adaptive = TRUE),
        FLOORSZ = mgwr_config(adaptive = TRUE),
        UNEMPLOY = mgwr_config(adaptive = TRUE),
        PROF = mgwr_config(adaptive = TRUE)
      )
    )
  })
  #### use default value
  expect_no_error({
    gwr_multiscale(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        FLOORSZ = mgwr_config(bw = 30, optim_bw = "no", adaptive = TRUE),
        .default = mgwr_config(adaptive = TRUE)
      )
    )
  })
  #### error when default is missing and config is missing
  expect_error({
    gwr_multiscale(
      formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
      data = LondonHP,
      config = list(
        FLOORSZ = mgwr_config(bw = 30, optim_bw = "no", adaptive = TRUE)
      )
    )
  }, "Either specific configs for all variables, or provide a default config naming '.default'!")
})

test_that("Multiscale GWR: verbose", {
  skip_on_ci()
  expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 1
  ))
  expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 2
  ))
  expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 3
  ))
  expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP,
    verbose = 4
  ))
})
