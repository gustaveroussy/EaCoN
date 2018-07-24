\name{ASCN.ASCAT.ff}
\alias{ASCN.ASCAT.ff}
\title{Allele-Specific Copy Number estimation using ASCAT, from a RDS file.}
\usage{
  ASCN.ASCAT.ff(RDS.file = NULL, ...)
}
\arguments{
  \item{RDS.file}{a RDS file generated with \code{Segment.ASCAT()} (or corresponding wrappers \code{Segment.ASCAT.ff()} and \code{Segment.ASCAT.ff.Batch()}.}
  \item{...}{Any parameter for the function \code{EaCoN.ASCN()}, except \code{data}.}
}
\description{
  This function is just a wrapper to the \code{ASCN.ASCAT()}, starting from a RDS file rather than a data object.).
}
\author{
  Bastien Job
}

