usethis::use_package("data.table")
usethis::use_package("sf")
usethis::use_package("sp")
usethis::use_package("gstat")
usethis::use_package("dplyr")
usethis::use_package("magic")
usethis::use_package("methods")
usethis::use_package("raster")
usethis::use_package("spdep")
usethis::use_package("tidyr")
usethis::use_pipe() # can use %>% after this

usethis::use_mit_license("Taro Mieno")

devtools::document() # regenerate documents reflecting the changes and apply load_all()
