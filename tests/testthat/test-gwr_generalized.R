data(LondonHP)
m <- NULL

test_that("GWR Generalized: poisson family", {
  m <<- expect_no_error(
    gwr_generalized(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, "poisson", 64, TRUE)
  )
})

test_that("GWR Generalized: binomial family", {
  m <<- expect_no_error(
    gwr_generalized(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, "binomial", 64, TRUE)
  )
})

test_that("GWR Generalized: Bandwidth selection", {
  m <<- expect_no_error(
    gwr_generalized(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, "binomial", bw="CV", TRUE)
  )
})
