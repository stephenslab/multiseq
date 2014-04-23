#Instructions to build the package

#require(devtools)
#build("multiseq/package/multiseq",binary=FALSE)
#instead of building the package
#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("multiseq", path=".", code <- files=c("multiseq.R", "PoissonBinomial.funcs.R", "deltamethod.R"), cpp <- files="multiseq.cpp", example <- code = FALSE, attributes = TRUE)
#in the terminal
#cd ..
#tar -pczf multiseq.tar.gz multiseq
#in R
install.packages("~/src/multiseq/package/multiseq.tar.gz",repos=NULL,type="source")
