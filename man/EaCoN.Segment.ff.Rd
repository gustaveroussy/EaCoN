\name{EaCoN.Segment.ff}
\alias{EaCoN.Segment.ff}
\title{L2R and BAF joint segmentation using ASCAT, from a RDS file.}
\usage{
  EaCoN.Segment.ff(RDS.file = NULL, ...)
}
\arguments{
  \item{RDS.file}{a RDS file generated with \code{EaCoN.WES.Normalize()} (for WES data), \code{EaCoN.CS.Process()} (for Affymetrix CytoScan data), \code{EaCoN.OS.Process()} (for Affymetrix OncoScan data), or \code{EaCoN.SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
\item{...}{Any parameter for the function \code{EaCoN.Segment()}, except \code{data}.}
}
\description{
  This function is just a wrapper to the \code{EaCoN.Segment()}, starting from a RDS file rather than a data object.).
}
\author{
  Bastien Job
}

