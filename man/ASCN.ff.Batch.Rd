\name{ASCN.ff.Batch}
\alias{ASCN.ff.Batch}
\title{Allele-Specific Copy Number estimation, from RDS files in batch mode, with multithreading.}
\usage{
  ASCN.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = "\\\\.SEG\\\\..*\\\\.RDS$", full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with one of the segmentation functions (\code{Segment.ASCAT()}, \code{Segment.FACETS()}, \code{Segment.SEQUENZA()}), or their corresponding ".ff" or ".ff.Batch" wrappers.}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter to the copy number estimation function, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{ASCN.ff()} with multithreading).
}
\author{
  Bastien Job
}

