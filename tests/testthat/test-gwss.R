data(LondonHP)
m <- NULL

test_that("gw.average: works", {
  m <<- expect_no_error({
    gw.average(LondonHP, c("PURCHASE","FLOORSZ","UNEMPLOY"), 64, TRUE)
  })
})

test_that("gw.correlation: works", {
  m <<- expect_no_error({
    gw.correlation(LondonHP, c("PURCHASE"),c("FLOORSZ","UNEMPLOY"), adaptive=TRUE)
  })
})

test_that("gw.correlation: works", {
  m <<- expect_no_error({
    gw.correlation(LondonHP, c("PURCHASE"),c("FLOORSZ","UNEMPLOY"),c(50,50), adaptive=TRUE)
  })
})