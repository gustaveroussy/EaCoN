\name{EaCoN.WES.Normalize}
\alias{EaCoN.WES.Normalize}
\title{Normalization of binned WES data.}
\usage{
  EaCoN.WES.Normalize(data = NULL, BINpack = NULL, L2R.RD.min.Ref = 20,
  L2R.RD.min.Test = 20, BAF.RD.min = 25, out.dir = getwd(), return.data = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.WES.Bin()}.}
  \item{BINpack}{A data pack generated with \code{EaCoN.BINpack.Maker()} containing to a binned version of the capture BED used to generate the BAM files, and corresponding precomputed GC\% tracks.}
  \item{L2R.RD.min.Ref}{Minimum read depth in the reference to keep a log2(ratio) bin.}
  \item{L2R.RD.min.Test}{Minimum read depth in the test to keep a log2(ratio) bin.}
  \item{BAF.RD.min}{Minimum read depth in the test to keep a BAF position.}
  \item{out.dir}{Output directory.}
  \item{return.data}{Return the results.}
}
\description{
  This function performs the normalization of the L2R and BAF signals from a data object obtained with \code{EaCoN.WES.Bin()} (or \code{EaCoN.WES.Bin.Batch()}).
}
\value{
  A RDS file that can be used with \code{EaCoN.Segment.ff()}.
}
\author{
  Bastien Job
}

