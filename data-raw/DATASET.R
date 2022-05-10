## code to prepare `FedYieldCurve` dataset
FedYieldCurve = read.csv("FedYieldCurve.csv", sep=";")
FedYieldCurve$Month = paste0(FedYieldCurve$Month, "-01")
FedYieldCurve$Month = base::as.Date(FedYieldCurve$Month)
usethis::use_data(FedYieldCurve, overwrite = TRUE)
