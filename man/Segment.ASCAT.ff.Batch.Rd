\name{Segment.ASCAT.ff.Batch}
\alias{Segment.ASCAT.ff.Batch}
\title{L2R and BAF joint segmentation using ASCAT, from a RDS file in batch mode, with multithreading.}
\usage{
  Segment.ASCAT.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = "_processed.RDS$",  full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{WES.Normalize()} (for WES data), \code{CS.Process()} (for Affymetrix CytoScan data), \code{OS.Process()} (for Affymetrix OncoScan data), or \code{SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{Segment.ASCAT()}, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{Segment.ASCAT.ff()} with multithreading).
}
\author{
  Bastien Job
}

