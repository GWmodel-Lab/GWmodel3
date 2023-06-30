data(LondonHP)
m <- NULL

test_that("Multiscale GWR: works", {
  m <<- expect_no_error(gwr_multiscale(
    formula = PURCHASE ~ FLOORSZ + UNEMPLOY + PROF,
    data = LondonHP
  ))
})
