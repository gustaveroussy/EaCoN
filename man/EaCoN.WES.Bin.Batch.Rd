\name{EaCoN.WES.Bin.Batch}
\alias{EaCoN.WES.Bin.Batch}
\title{Run EaCoN.WES.Bin() in batch mode with multithreading}
\description{
  This is just a wrapper to EaCoN.WES.Bin() with multithreading.
}
\usage{
  EaCoN.WES.Bin.Batch(BAM.list.file = NULL, BINpack = NULL, nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{BAM.list.file}{A tab-separated file giving the file names of the test and ref BAM files and corresponding sample name (see Details section).}
  \item{BINpack}{A data pack generated with \code{EaCoN.BINpack.Maker()} containing to a binned version of the capture BED used to generate the BAM files, and corresponding precomputed GC\% tracks.}
  \item{nthread}{Number of threads for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{EaCoN.WES.Bin()}, except \code{testBAM}, \code{refBAM} and \code{samplename}.}
}
\details{
  \code{BAM.list.file} is a tab-separated text file containing 3 columns (header and specified column names are mandatory) :
    \itemize{
      \item{testBAM : Name (and path) of the test BAM file(s)}
      \item{refBAM : Name (and path) of the reference BAM file(s)}
      \item{SampleName : The output sample name(s)}
    }
}
\author{
  Bastien Job
}

