Sys.setenv("R_TESTS" = "")
library(testthat)
library(spatzie)
test_check("spatzie")
