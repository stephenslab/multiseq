\name{get.counts}
\alias{get.counts}
\title{Get sample counts in a genomic region from bam, hdf5, or bigWig file}
\usage{
get.counts(samples, chr, locus.start, locus.end)
}
\arguments{
  \item{samples:}{a data frame extracted from the
  samplesheet with N rows}

  \item{chr:}{a string representing a reference sequence
  name}

  \item{locus.start:}{an integer representing a position on
  chr}

  \item{locus.end:}{an integer representing a position on
  chr (must be locus.end>locus.start)}
}
\value{
a matrix with N rows and locus.end-locus.start+1 columns
containing the number of reads that start at each base in
the specified region in each sample
}
\description{
If samples$bamReads is specified then this function
extracts reads from the bam files in samples$bamReads using
`samtools` (which must be in the USER's path). Else if
samples$h5FullPath is specified this function extracts
reads from the hdf5 files in samples$h5FullPath using the R
package rhdf5. Else if samples$bigWigPath is specified this
function extracts reads from bigWig files using the
executable `bigWigToWig` (which must be in the USER's
path).
}

