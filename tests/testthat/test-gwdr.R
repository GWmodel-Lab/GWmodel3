data(LondonHP)
m <- NULL

test_that("GWDR: works", {
  m <<- expect_no_error(gwdr(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP))
})
