data(LondonHP)
m <- NULL

test_that("GWRScalable: works", {
  m <<- expect_no_error({
    gwr_scalable(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64)
  })
})

# test_that("GWRScalable: bw", {
#   m <<- expect_no_error({
#     gwr_scalable(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP)
#   })
# })