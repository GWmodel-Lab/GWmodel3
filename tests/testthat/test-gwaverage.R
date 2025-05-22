data(LondonHP)
m <- NULL

test_that("GWAverage: works", {
  m <<- expect_no_error({
    gwaverage(~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  })
})