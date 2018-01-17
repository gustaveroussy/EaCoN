\name{EaCoN.Annotate.ff.Batch}
\alias{EaCoN.Annotate.ff.Batch}
\title{Runs EaCoN.Annotate.ff in batch mode with multithreading}
\description{
  This is just a wrapper to EaCoN.Annotate.ff() with multithreading.
}
\usage{
  EaCoN.Annotate.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = ".EaCoN.ASPCF.RDS$", full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{EaCoN.Segment()}.}
  \item{nthread}{Number of thread(s) for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any other parameter for \code{EaCoN.Annotate())} except \code{data}.}
}
\author{
  Bastien Job
}

