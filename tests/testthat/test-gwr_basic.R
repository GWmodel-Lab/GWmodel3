data(LondonHP)
m <- NULL

test_that("Basic GWR: works", {
  m <<- expect_no_error(
    gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})

test_that("Basic GWR: Bandwidth selection", {
  m <<- expect_no_error(
    gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE)
  )
})

test_that("Basic GWR: Variable selection", {
  m1 <<- expect_no_error(
    model_sel(m)
  )
})

test_that("Basic GWR: predict", {
  predict(m, LondonHP)
})

test_that("Basic GWR: verbose", {
  expect_no_error(
    model_sel(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE, verbose = 1))
  )
  expect_no_error(
    model_sel(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE, verbose = 2))
  )
  expect_no_error(
    predict(m, LondonHP, verbose = 1)
  )
})
