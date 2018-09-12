\name{ASCN.SEQUENZA}
\alias{ASCN.SEQUENZA}
\title{Allele-Specific Copy Number estimation using SEQUENZA}
\usage{
  ASCN.SEQUENZA(data = NULL, max.ploidy = 4, ploidy.step = .1,
  seg.min.size = 1E+06, out.dir = getwd(), force = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{Segment.FACETS()} (or corresponding wrappers \code{Segment.ff()} and \code{Segment.ff.Batch()}, using the "FACETS" segmenter.}
  \item{max.ploidy}{Maximal median ploidy value to evaluate.}
  \item{ploidy.step}{Step for the ploidy value to evaluate in loop.}
  \item{seg.min.size}{Minimal size (in nt) of a segment to be used in the model fit.}
  \item{out.dir}{Output directory.}
  \item{force}{If some result files already exist, delete them before runing.}
}
\description{
  This function performs the estimation of allele-specific copy number using SEQUENZA.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables.
}
\author{
  Bastien Job
}

