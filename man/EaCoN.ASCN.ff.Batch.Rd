\name{EaCoN.ASCN.ff.Batch}
\alias{EaCoN.ASCN.ff.Batch}
\title{Allele-Specific Copy Number estimation using ASCAT, from a RDS file in batch mode, with multithreading.}
\usage{
  EaCoN.ASCN.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = ".EaCoN.ASPCF.RDS$", full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE),  nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{EaCoN.Segment()} (or corresponding wrappers \code{EaCoN.Segment.ff()} and \code{EaCoN.Segment.ff.Batch()}.}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{EaCoN.ASCN()}, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{EaCoN.WES.ASCN.ff()} with multithreading).
}
\author{
  Bastien Job
}

