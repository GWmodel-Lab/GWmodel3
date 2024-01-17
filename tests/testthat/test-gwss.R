data(LondonHP)
m <- NULL

test_that("GWSS: works", {
  m <<- expect_no_error({
    gwss(~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE)
  })
})