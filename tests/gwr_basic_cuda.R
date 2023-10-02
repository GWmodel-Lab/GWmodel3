library(GWmodel3, lib.loc = "../Rlibrary.tmp")
data(LondonHP)
m_cuda <- gwr_basic(PURCHASE~FLOORSZ+UNEMPLOY, LondonHP, 64, TRUE, parallel_method = "cuda", parallel_arg = c(0, 64), verbose = 1)
predict(m_cuda, LondonHP, verbose = 1)

# whhp <- with(wuhan.hp, sf::st_as_sf(data, coords = c("lon", "lat")))
# gwr_basic(Price ~ d.PrimarySchool + BuildingArea + Fee, whhp, bw = 64, parallel_method = "cuda", parallel_arg = c(0, 64))