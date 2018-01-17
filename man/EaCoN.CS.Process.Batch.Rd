\name{EaCoN.CS.Process.Batch}
\alias{EaCoN.CS.Process.Batch}
\title{Run EaCoN.CS.Process() in batch mode with multithreading}
\description{
  This is just a wrapper to EaCoN.CS.Process() with multithreading.
}
\usage{
  EaCoN.CS.Process.Batch(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...)
}
\arguments{
  \item{CEL.list.file}{A tab-separated file giving the file names of the CEL files and corresponding sample name (see Details section).}
  \item{nthread}{Number of threads for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{EaCoN.CS.Process()}, except \code{CEL} and \code{samplename}.}
}
\details{
  \code{CEL.list.file} is a tab-separated text file containing 2 columns (header and specified column names are mandatory) :
    \itemize{
      \item{cel_files : Name (and path) of the CEL file(s)}
      \item{SampleName : The output sample name(s)}
    }
}
\author{
  Bastien Job
}

