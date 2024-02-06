data(LondonHP)
m <- NULL

test_that("SDR: works", {
  m <<- expect_no_error(sdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP))
})

test_that("SDR: verbose", {
  skip_on_ci()
  expect_no_error(sdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 1))
  expect_no_error(sdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 2))
})
