\name{EaCoN.Segment.ff.Batch}
\alias{EaCoN.Segment.ff.Batch}
\title{L2R and BAF joint segmentation using ASCAT, from a RDS file in batch mode, with multithreading.}
\usage{
  EaCoN.Segment.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = "_processed.RDS$",  full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{EaCoN.WES.Normalize()} (for WES data), \code{EaCoN.CS.Process()} (for Affymetrix CytoScan data), \code{EaCoN.OS.Process()} (for Affymetrix OncoScan data), or \code{EaCoN.SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{EaCoN.Segment()}, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{EaCoN.WES.Segment.ff()} with multithreading).
}
\author{
  Bastien Job
}

