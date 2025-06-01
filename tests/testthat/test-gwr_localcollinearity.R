data(LondonHP)
m <- NULL

test_that("GWR LocalCollinearity: works", {
  m1 <<- expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})

test_that("GWR LocalCollinearity: bw selection", {
  expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP)
  )
})

test_that("GWR LocalCollinearity: omp", {
  expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, parallel_method = "omp", parallel_arg = 4)
  )
})

test_that("GWR LocalCollinearity: lambda & cnthresh", {
  m2 <<- expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, lambda = 0.3, cn_thresh = 30)
  )
})

test_that("GWR LocalCollinearity: lambda & cnthresh", {
  m3 <<- expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, lambda_adjust = TRUE,cn_thresh = 30)
  )
})

test_that("GWR LocalCollinearity: lambda error", {
  expect_error({
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, lambda = 10, cn_thresh = 30)
  }, "Error: lambda must in \\[0,1\\]")
})
