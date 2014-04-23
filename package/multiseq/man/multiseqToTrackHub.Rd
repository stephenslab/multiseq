\name{multiseqToTrackHub}
\alias{multiseqToTrackHub}
\title{Create a UCSC genome browser "Track Hub" from multiseq output.

This funxtion requires the executables
"wigToBigWig" and "bedToBigBed" to be in the user's PATH}
\usage{
multiseqToTrackHub(region, hub_name = "multiseq",
  multiseq_folder = "./results_run_multiseq/",
  chrom_file = "~/src/multiseq/data/chromosome.lengths.hg19.txt",
  assembly = "hg19", mountpoint = MOUNTPOINT_PATH,
  http_address = MOUNTPOINT_HTTP_ADDRESS)
}
\arguments{
  \item{region:}{a region (e.g. chr1:2345-234567)}

  \item{hub_name:}{name of the track hub; it could contain
  a path, in which case the path will be relative to the
  mountpoint. This string can be set to any value;
  default="multiseq"}

  \item{multiseq_folder:}{path to the folder containing
  results from multiseq; this script requires output from
  multiseq to be in the format effect_mean_var.txt.gz where
  first column is effect and second column variance}

  \item{chrom_file:}{path to the file containing chromosome
  names and lengths;
  default="~/src/multiseq/data/chromosome.lengths.hg19.txt"}

  \item{assembly:}{genome assembly that reads were mapped
  to; default="hg19"}

  \item{mountpoint:}{path to the directory where the track
  hub folder will be saved in. This directory should be
  associated with an http address or an ftp address;
  default=Sys.getenv("MOUNTPOINT_PATH")}

  \item{http_address:}{http or ftp address associated with
  the mountpoint; default=Sys.getenv("MOUNTPOINT_PATH")}
}
\description{
Create a UCSC genome browser "Track Hub" from multiseq
output.

This funxtion requires the executables "wigToBigWig" and
"bedToBigBed" to be in the user's PATH
}

