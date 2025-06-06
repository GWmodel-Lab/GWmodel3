data(LondonHP)

test_that("GWDA: wqda", {
  m1 <<- expect_no_error({
    gwda(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  })
})

test_that("GWDA: wlda", {
  m2 <<- expect_no_error({
    gwda(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, method = "wlda")
  })
})

test_that("GWDA: omp", {
  expect_no_error({
    gwda(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, parallel_method = "omp", parallel_arg = 2)
  })
})