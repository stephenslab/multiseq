% Generated by roxygen2 (4.0.2): do not edit by hand
\name{split_region}
\alias{split_region}
\title{Split a region string into sequence name, region start position, and region end position.}
\usage{
split_region(region)
}
\arguments{
\item{region:}{a string specifying a genomic region: reference sequence name, start position, end position}
}
\value{
a list with elements \code{chr}, \code{start}, \code{end}.
}
\description{
Split a region string into sequence name, region start position, and region end position.
}
\examples{
region=split_region("chr1:11740409-11756792")
print(region)
print(region$chr)
print(region$start)
}
\keyword{internal}

