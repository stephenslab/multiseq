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

######################
# Instruction to build the package using Rstudio
######################
1. open Rstudio
2. specify project : file -> new project -> existing directory ->  specify "multiseq" directory under "package" directory 
3. install all necessary packages 
4. require(roxygen2)
5. require(devtools)
6. roxygenize()
7. document()
8. build -> build source package 
9. before push, be careful! You don't want to push some unnecessary files created during building package.



#in the terminal
#cd ..
#tar -pczf multiseq.tar.gz multiseq
#in R
install.packages("~/src/multiseq/package/multiseq.tar.gz",repos=NULL,type="source")
