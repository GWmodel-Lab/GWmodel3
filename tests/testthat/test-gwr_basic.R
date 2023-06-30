data(LondonHP)
m <- NULL

test_that("Basic GWR: works", {
  m <<- expect_no_error(
    gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})
