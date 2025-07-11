\name{Segment.ff}
\alias{Segment.ff}
\title{L2R and BAF joint segmentation, from a RDS file.}
\usage{
  Segment.ff(RDS.file = NULL, segmenter = "ASCAT", ...)
}
\arguments{
  \item{RDS.file}{A RDS file generated with a normalization function : \code{WES.Normalize()} (for WES data), \code{CS.Process()} (for Affymetrix CytoScan data), \code{OS.Process()} (for Affymetrix OncoScan data), or \code{SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{segmenter}{The (case sensitive!) name of the segmenter ("ASCAT" or "FACETS")}
  \item{...}{Any parameter to the segmentation function refered by the \code{segmenter} name, except \code{data}.}
}
\description{
  This function is just a wrapper to any of the segmentation functions (\code{Segment.ASCAT()}, \code{Segment.FACETS()}), starting from a RDS file rather than a data object.
}
\author{
  Bastien Job
}

