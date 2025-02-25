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

test_that("Basic GWR: missing values", {
  expect_no_error({
    sub_data <- LondonHP[, c("PURCHASE", "FLOORSZ", "UNEMPLOY")]
    sub_data$UNEMPLOY[1] <- NA
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
  skip_on_ci()
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

test_that("Basic GWR: CUDA", {
  if (Sys.getenv("CUDA_PATH") == "") {
    expect_error({
      m_cuda <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, parallel_method = "cuda", parallel_arg = c(0, 64))
      predict(m_cuda, LondonHP)
    })
  } else {
    expect_no_error({
      m_cuda <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, parallel_method = "cuda", parallel_arg = c(0, 64))
      predict(m_cuda, LondonHP)
    })
  }
})

test_that("Basic GWR: CUDA for Windows development", {
  skip("This test case is only used for development on Windows")
  library(GWmodel3, lib.loc = "../Rlibrary.tmp")
  data(LondonHP)
  expect_no_error({
    m_cuda <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, parallel_method = "cuda", parallel_arg = c(0, 64), verbose = 1)
    predict(m_cuda, LondonHP, verbose = 1)
  })
})
