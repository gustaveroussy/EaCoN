\name{EaCoN.ASCN.ff}
\alias{EaCoN.ASCN.ff}
\title{Allele-Specific Copy Number estimation using ASCAT, from a RDS file.}
\usage{
  EaCoN.ASCN.ff(RDS.file = NULL, ...)
}
\arguments{
  \item{RDS.file}{a RDS file generated with \code{EaCoN.Segment()} (or corresponding wrappers \code{EaCoN.Segment.ff()} and \code{EaCoN.Segment.ff.Batch()}.}
  \item{...}{Any parameter for the function \code{EaCoN.ASCN()}, except \code{data}.}
}
\description{
  This function is just a wrapper to the \code{EaCoN.ASCN()}, starting from a RDS file rather than a data object.).
}
\author{
  Bastien Job
}

