data(LondonHP)
m <- NULL

test_that("GWPCA: works", {
  m <<- expect_no_error(gwpca(~PURCHASE + FLOORSZ + UNEMPLOY + PROF, LondonHP, bw = 50, adaptive = TRUE, components = 3))
})

test_that("GWPCA: error test", {
  expect_error(gwpca(~FLOORSZ+UNEMPLOY, LondonHP, bw = 50, adaptive = TRUE, components = 3),
  "Components to keep must be lower than variable counts!")
})

test_that("GWPCA: error components", {
  expect_error(gwpca(~FLOORSZ+UNEMPLOY, LondonHP, bw = 50, adaptive = TRUE, components = 0),
  "Components must be an interger greater than 0 !!")
})