# precompile vignettes that depend on internet access

setwd("./vignettes")
knitr::knit("./individual-steps.Rmd.orig", "./individual-steps.Rmd")
knitr::knit("./single-call.Rmd.orig", "./single-call.Rmd")
setwd("../")
