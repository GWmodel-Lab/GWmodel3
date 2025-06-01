data(LondonHP)
m <- NULL

test_that("Basic GWR LocalCollinearity: works", {
  m <<- expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})

test_that("Basic GWR LocalCollinearity: bw selection", {
  m <<- expect_no_error(
    gwr_lcr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP)
  )
})