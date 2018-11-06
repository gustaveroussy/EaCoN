\name{ASCN.ff}
\alias{ASCN.ff}
\title{Allele-Specific Copy Number estimation, from a RDS file.}
\usage{
  ASCN.ff(RDS.file = NULL, ...)
}
\arguments{
  \item{RDS.file}{A RDS file generated with one of the segmentation functions (\code{Segment.ASCAT()}, \code{Segment.FACETS()}), or their corresponding ".ff" or ".ff.Batch" wrappers.}
  \item{...}{Any parameter to the copy number estimation method.}
}
\description{
  This function is just a wrapper to any of the copy number estimation functions (\code{ASCN.ASCAT()}, \code{ASCN.FACETS()} and \code{ASCN.SEQUENZA()}), starting from a RDS file rather than a data object.
}
\author{
  Bastien Job
}

