\name{Segment.ASCAT.ff}
\alias{Segment.ASCAT.ff}
\title{L2R and BAF joint segmentation using ASCAT, from a RDS file.}
\usage{
  Segment.ASCAT.ff(RDS.file = NULL, ...)
}
\arguments{
  \item{RDS.file}{a RDS file generated with \code{WES.Normalize()} (for WES data), \code{CS.Process()} (for Affymetrix CytoScan data), \code{OS.Process()} (for Affymetrix OncoScan data), or \code{SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
\item{...}{Any parameter for the function \code{Segment.ASCAT()}, except \code{data}.}
}
\description{
  This function is just a wrapper to the \code{Segment.ASCAT()}, starting from a RDS file rather than a data object.).
}
\author{
  Bastien Job
}

