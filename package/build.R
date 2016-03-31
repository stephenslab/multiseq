#Instructions to build the package

#require(devtools)
#build("multiseq/package/multiseq",binary=FALSE)
#instead of building the package
#put all files (R and cpp files) in a directory and cd into that directory
#in R:
library(Rcpp)
Rcpp.package.skeleton("multiseq", path=".", code_files=c("multiseq.R", "PoissonBinomial.funcs.R", "deltamethod.R"), cpp_files="multiseq.cpp", example_code = FALSE, attributes = TRUE)

######################
# Instruction to document the package using R
######################
#from within the package/multiseq directory:
#/data/tools/R-3.1.1/bin/R
#devtools::install_github("klutometis/roxygen")
library(ashr)
require(roxygen2)
require(devtools)
roxygenize()
document()
#log out of R
#generate an archive from folder multiseq
#tar -pczf multiseq_0.1.tar.gz multiseq
#log into R from within the package/multiseq directory
#/data/tools/R-3.1.1/bin/R
install.packages("~/src/multiseq/package/multiseq_0.1.tar.gz",repos=NULL,type="source")
#make sure that path.package("multiseq") is
#[1] "/mnt/lustre/home/epantaleo/R/x86_64-unknown-linux-gnu-library/3.1/multiseq"
#before you build the vignette
require(devtools)
library(multiseq)
build_vignettes()
#in bash
#generate an archive from folder multiseq
#tar -pczf multiseq_0.1.tar.gz multiseq
install.packages("~/src/multiseq/package/multiseq_0.1.tar.gz",repos=NULL,type="source")
help(package="multiseq")

#in Rstudio
Sys.setenv(MOUNTPOINT_PATH="/data/internal/solexa_mountpoint/epantaleo")
Sys.setenv(MOUNTPOINT_HTTP_ADDRESS="https://*******:pipeline@solexa-compute1.uchicago.edu/stephenslab/epantaleo")
Sys.setenv(PATH=paste0(Sys.getenv('PATH'),":/usr/local/bin/",":/data/tools/ucsctools/"))
