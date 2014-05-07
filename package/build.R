#Instructions to build the package

#require(devtools)
#build("multiseq/package/multiseq",binary=FALSE)
#instead of building the package
#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("multiseq", path=".", code <- files=c("multiseq.R", "PoissonBinomial.funcs.R", "deltamethod.R"), cpp <- files="multiseq.cpp", example <- code = FALSE, attributes = TRUE)

#/data/tools/R-3.0.3/bin/R
#devtools::install_github("klutometis/roxygen")
require(roxygen2)
require(devtools)
roxygenize("~/Downloads/multiseq")
document("~/Downloads/multiseq")




#in the terminal
#cd ..
#tar -pczf multiseq.tar.gz multiseq
#in R
install.packages("~/src/multiseq/package/multiseq.tar.gz",repos=NULL,type="source")
