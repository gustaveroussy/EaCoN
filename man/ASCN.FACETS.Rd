\name{ASCN.FACETS}
\alias{ASCN.FACETS}
\title{Allele-Specific Copy Number estimation using FACETS}
\usage{
  ASCN.FACETS(data = NULL, out.dir = getwd(), force = FALSE, ...)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{Segment.FACETS()} (or corresponding wrappers \code{Segment.ff()} and \code{Segment.ff.Batch()}, using the "FACETS" segmenter.}
  \item{out.dir}{Output directory.}
  \item{force}{If some result files already exist, delete them before runing.}
  \item{...}{Any additional parameter to give to \code{facets::emcncf()}.}
}
\description{
  This function performs the estimation of allele-specific copy number using FACETS.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables.
}
\author{
  Bastien Job
}

