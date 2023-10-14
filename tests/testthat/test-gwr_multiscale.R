data(LondonHP)
m <- NULL

test_that("Multiscale GWR: works", {
  m <<- expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP
  ))
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
