\name{EaCoN.ASCN}
\alias{EaCoN.ASCN}
\title{Allele-Specific Copy Number estimation using ASCAT.}
\usage{
  EaCoN.ASCN(data = NULL, gammaRange = c(.35,.95), out.dir = getwd())
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.Segment()} (or corresponding wrappers \code{EaCoN.Segment.ff()} and \code{EaCoN.Segment.ff.Batch()}.}
  \item{gammaRange}{The range of gamma values to assess the modelling of ASCN (see \code{?ASCAT::ascat.runAscat}).}
  \item{out.dir}{Output directory.}
}
\description{
  This function performs the estimation of allele-specific copy number using ASCAT.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables in several folders containing the different ASCN (allele-specific copy number) analyses performed for each value in \code{gammaRange}. The directory corresponding to the best evaluated gamma will also have an "_optimal" suffix.
}
\details{
  Please note that depending on noise and complexity, ASCAT will not always be able to perform ASCN.
}
\author{
  Bastien Job
}

