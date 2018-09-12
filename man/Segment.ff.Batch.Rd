\name{Segment.ff.Batch}
\alias{Segment.ff.Batch}
\title{L2R and BAF joint segmentation, from a RDS file in batch mode, with multithreading.}
\usage{
  Segment.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = "_processed.RDS$",  full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), segmenter = "ASCAT",
  nthread = 1, cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with a normalization function : \code{WES.Normalize()} (for WES data), \code{CS.Process()} (for Affymetrix CytoScan data), \code{OS.Process()} (for Affymetrix OncoScan data), or \code{SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{segmenter}{The (case sensitive!) name of the segmenter ("ASCAT" or "FACETS")}
  \item{nthread}{Number of threads.}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter to the segmentation function refered by the \code{segmenter} name, except \code{data}.}
}
\description{
  This function is just a wrapper to \code{Segment.ff()} with multithreading).
}
\author{
  Bastien Job
}

