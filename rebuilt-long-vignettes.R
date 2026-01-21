old_wd <- getwd()

setwd("vignettes/")
knitr::knit("ALDEx3-Quickstart.Rmd.orig", output = "ALDEx3-Quickstart.Rmd")

setwd(old_wd)
