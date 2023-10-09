data(LondonHP)
m <- NULL

test_that("Basic GWR: works", {
  m <<- expect_no_error(
    gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})

test_that("Basic GWR: formula with dot", {
  expect_no_error({
    sub_data <- LondonHP[, c("PURCHASE", "FLOORSZ", "UNEMPLOY")]
    gwr_basic(PURCHASE~., data = sub_data)
  })
})

test_that("Basic GWR: Bandwidth selection", {
  m <<- expect_no_error(
    gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE)
  )
})

test_that("Basic GWR: Variable selection", {
  m1 <<- expect_no_error(
    step(m)
  )
})

test_that("Basic GWR: predict", {
  predict(m, LondonHP)
})

test_that("Basic GWR: verbose", {
  expect_no_error(
    step(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE, verbose = 1))
  )
  expect_no_error(
    step(gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE, verbose = 2))
  )
  expect_no_error(
    predict(m, LondonHP, verbose = 1)
  )
})
