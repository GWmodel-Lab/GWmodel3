data(LondonHP)
m <- NULL

test_that("GWDR: works", {
  m <<- expect_no_error(gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP))
})

test_that("GWDR: verbose", {
  skip_on_ci()
  expect_no_error(gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 1))
  expect_no_error(gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 2))
})
