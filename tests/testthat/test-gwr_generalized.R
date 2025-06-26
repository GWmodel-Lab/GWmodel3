data(LondonHP)
# m <- NULL

test_that("GWR Generalized: poisson family", {
  m1 <<- expect_no_error(
    gwr_generalized(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "poisson", 64, TRUE)
  )
})

test_that("GWR Generalized: binomial family", {
  expect_no_error(
    gwr_generalized(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, "binomial", 64, TRUE)
  )
})

test_that("GWR Generalized: Bandwidth selection", {
  m2 <<- expect_no_error(
    gwr_generalized(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, "binomial", bw="CV", TRUE)
  )
})

test_that("GWR Generalized: Bandwidth selection, omp", {
  expect_no_error(
    gwr_generalized(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, "binomial", bw="AIC", TRUE, parallel_method = "omp", parallel_arg = 4)
  )
})

# 测试问题可能在内核库的predict函数中，应该先设置一下mHasRegressionData，否则导致现在出现矩阵列数计算的错误。
# 在内核库中测试通过predict后，再完善这部分内容。
# test_that("GWR Generalized: predict", {
#   predict(m1, LondonHP)
# })