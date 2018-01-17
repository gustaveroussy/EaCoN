\name{EaCoN.OS.Process.Batch}
\alias{EaCoN.OS.Process.Batch}
\title{Runs EaCoN.OS.Process in batch mode with multithreading}
\description{
  This is just a wrapper to EaCoN.OS.Process() with multithreading.
}
\usage{
  EaCoN.OS.Process.Batch(pairs.file = NULL, nthread = 1, cluster.type = "PSOCK", ...)
}
\arguments{
  \item{pairs.file}{A tab-separated file giving the file names of the pair(s) of OncoScan CEL files and sample name.}
  \item{nthread}{Number of thread(s) for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any other parameter for \code{EaCoN.OS.Process())} except \code{ATChannelCel}, \code{GCChannelCel} and \code{samplename}.}
}
\details{
  \code{pairs.file} is a tab-separated text file containing 3 columns (header and specified column names are mandatory) :
  \itemize{
    \item{ATChannelCel : Name (and path) of the A-channel CEL file(s)}
    \item{GCChannelCel : Name (and path) of the C-channel CEL file(s)}
    \item{SampleName : The output sample name(s)}
  }
}
\author{
Bastien Job
}

