data(LondonHP)
m <- NULL

test_that("GWR Robust: works", {
  m <<- expect_no_error(
    gwr_robust(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  )
})

test_that("GWR Robust: filtered", {
  m <<- expect_no_error(
    gwr_robust(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, filter=TRUE)
  )
})

test_that("GWR Robust: formula with dot", {
  expect_no_error({
    sub_data <- LondonHP[, c("PURCHASE", "FLOORSZ", "UNEMPLOY")]
    gwr_robust(PURCHASE~., data = sub_data)
  })
})

test_that("GWR Robust: missing values", {
  expect_no_error({
    sub_data <- LondonHP[, c("PURCHASE", "FLOORSZ", "UNEMPLOY")]
    sub_data$UNEMPLOY[1] <- NA
    gwr_robust(PURCHASE~., data = sub_data)
  })
})

test_that("GWR Robust: Bandwidth selection", {
  m <<- expect_no_error(
    gwr_robust(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "AIC", TRUE)
  )
})

# test_that("GWR Robust: Variable selection", {
#   m1 <<- expect_no_error(
#     step(m)
#   )
# })