data(LondonHP)
m <- NULL

test_that("GWDA: works", {
  m <<- expect_no_error({
    gwda(TYPEDETCH~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  })
})