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
# Instruction to document the package using Rstudio
######################
1. open Rstudio
2. specify project : file -> new project -> existing directory ->  specify "multiseq" directory under "package" directory 
3. install all necessary packages 
4. require(roxygen2)
5. require(devtools)
6. roxygenize()
7. document()
8. build -> build source package 
9. before push, be careful! You dont want to push some unnecessary files created during building package.



#in the terminal
#cd ..
#tar -pczf multiseq.tar.gz multiseq
#IMPORTANT!!! /data/tools/R-3.1.1/bin/R
install.packages("~/src/multiseq/package/multiseq_0.1.tar.gz",repos=NULL,type="source")


#instructions to build the vignette using Rstudio (check version of R)
path.to.multiseq.gz="~/src/multiseq/package/multiseq_0.1.tar.gz"
path.to.ashr.gz="~/src/ash/package/ashr.no.cxx.tar.gz"
Sys.setenv(MOUNTPOINT_PATH="/data/internal/solexa_mountpoint/epantaleo")
Sys.setenv(MOUNTPOINT_HTTP_ADDRESS="https://*******:pipeline@solexa-compute1.uchicago.edu/stephenslab/epantaleo")
Sys.setenv(PATH=paste0(Sys.getenv('PATH'),":/usr/local/bin/",":/data/tools/ucsctools/"))
install.packages(path.to.ashr.gz,repos=NULL,type="source")
install.packages(path.to.multiseq.gz,repos=NULL,type="source")



#change data
#save(x, region, samples, file="~/src/multiseq/package/multiseq/data/example2.RData")
