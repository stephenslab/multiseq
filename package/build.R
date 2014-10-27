#Instructions to build the package

#require(devtools)
#build("multiseq/package/multiseq",binary=FALSE)
#instead of building the package
#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("multiseq", path=".", code_files=c("multiseq.R", "PoissonBinomial.funcs.R", "deltamethod.R"), cpp_files="multiseq.cpp", example_code = FALSE, attributes = TRUE)

#/data/tools/R-3.0.3/bin/R
#devtools::install_github("klutometis/roxygen")
require(roxygen2)
require(devtools)
roxygenize("~/Downloads/multiseq")
document("~/Downloads/multiseq")

######################
# Instruction to document the package using R
######################
install.packages("~/src/ash/package/ashr.no.cxx.tar.gz",repos=NULL,type="source")
require(roxygen2)
require(devtools)
roxygenize()
document()
build_vignettes()
#in bash
#generate and archive from folder multiseq
#tar -pczf multiseq.tar.gz multiseq
install.packages("~/src/multiseq/package/multiseq_0.1.tar.gz",repos=NULL,type="source")
help(package="multiseq")

#in Rstudio
Sys.setenv(MOUNTPOINT_PATH="/data/internal/solexa_mountpoint/epantaleo")
Sys.setenv(MOUNTPOINT_HTTP_ADDRESS="https://*******:pipeline@solexa-compute1.uchicago.edu/stephenslab/epantaleo")
Sys.setenv(PATH=paste0(Sys.getenv('PATH'),":/usr/local/bin/",":/data/tools/ucsctools/"))
