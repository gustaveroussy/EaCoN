\name{EaCoN.WES.Process}
\alias{EaCoN.WES.Process}
\title{Performs the binning, then the normalization of WES read depth data from BAM files.}
\usage{
  EaCoN.WES.Process(testBAM = NULL, refBAM = NULL, BINpack = NULL,
  samplename = "SAMPLE", Q = 20, L2R.RD.min.Ref = 10,
  L2R.RD.min.Test = 10, BAF.RD.min = 25,  out.dir = getwd(),
  return.data = FALSE)
}
\arguments{
  \item{testBAM}{BAM file corresponding to the tested sample (typically a patient tumor sample).}
  \item{refBAM}{BAM file corresponding to the reference sample (typically a patient normal sample).}
  \item{BINpack}{A data pack generated with \code{EaCoN.BINpack.Maker()} containing to a binned version of the capture BED used to generate the BAM files, and corresponding precomputed GC\% tracks.}
  \item{samplename}{The name of the sample.}
  \item{Q}{Minimum baseQ/BAQ quality score as defined by \code{samtools mpileup}.}
  \item{L2R.RD.min.Ref}{Minimum read depth in the reference to keep a log2(ratio) bin.}
  \item{L2R.RD.min.Test}{Minimum read depth in the test to keep a log2(ratio) bin.}
  \item{BAF.RD.min}{Minimum read depth in the test to keep a BAF position.}
  \item{out.dir}{Output directory.}
  \item{return.data}{Return the results.}
}
\description{
  This is just a wrapper function that chains \code{EaCoN.WES.Bin.ff()} then \code{EaCoN.WES.Normalize()}.
}
\author{
  Bastien Job
}

