data(LondonHP)
LondonHP$time <- as.integer(round(runif(nrow(LondonHP), 1, 365)))
m <- NULL

test_that("GTWR: set timestamps", {
  m <<- expect_no_error({
    gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "time", 64, 0.05, TRUE)
  })
  expect_no_error({
    gtwr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, LondonHP$time, 64, 0.05, TRUE)
  })
})

test_that("GTWR: helper functions", {
  expect_no_error({
    coef(m)
    fitted(m)
    residuals(m)
  })
  expect_no_error({
    plot(m)
  })
})
