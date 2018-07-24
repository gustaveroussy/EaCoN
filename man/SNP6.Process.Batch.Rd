\name{SNP6.Process.Batch}
\alias{SNP6.Process.Batch}
\title{Run SNP6.Process() in batch mode with multithreading}
\description{
  This is just a wrapper to SNP.Process() with multithreading.
}
\usage{
  SNP6.Process.Batch(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...)
}
\arguments{
  \item{CEL.list.file}{A tab-separated file giving the file names of the CEL files and corresponding sample name (see Details section).}
  \item{nthread}{Number of threads for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{SNP.Process()}, except \code{CEL}, \code{samplename} and \code{return.data}.}
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

