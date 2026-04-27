data(LondonHP)
m <- NULL

test_that("GTDR: works", {
  m <<- expect_no_error(gtdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP))
})

test_that("GTDR: verbose", {
  skip_on_ci()
  expect_no_error({
    gtdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 1)
  })
  expect_no_error({
    gtdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, optim_bw = "AIC", verbose = 2)
  })
})

test_that("GTDR: local periodical", {
  skip_on_ci()
  expect_no_error({
    gtdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, config = list(
      gtdr_config(),
      gtdr_config(kernel = "localperiodical", kernel_params = c(1))
    ), verbose = 1)
  })
})
