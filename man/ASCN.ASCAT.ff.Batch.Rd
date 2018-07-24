\name{ASCN.ASCAT.ff.Batch}
\alias{ASCN.ASCAT.ff.Batch}
\title{Allele-Specific Copy Number estimation using ASCAT, from a RDS file in batch mode, with multithreading.}
\usage{
  ASCN.ASCAT.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = ".EaCoN.ASPCF.RDS$", full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE),  nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{Segment.ASCAT()} (or corresponding wrappers \code{Segment.ASCAT.ff()} and \code{Segment.ASCAT.ff.Batch()}.}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{ASCN.ASCAT()}, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{ASCN.ASCAT.ff()} with multithreading).
}
\author{
  Bastien Job
}

