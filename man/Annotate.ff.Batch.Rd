\name{Annotate.ff.Batch}
\alias{Annotate.ff.Batch}
\title{Runs Annotate.ff in batch mode with multithreading}
\description{
  This is just a wrapper to Annotate.ff() with multithreading.
}
\usage{
  Annotate.ff.Batch(RDS.files = list.files(path = getwd(),
  pattern = "\\.SEG\\.ASCAT\\.RDS$", full.names = TRUE, recursive = TRUE,
  ignore.case = TRUE, include.dirs = FALSE), nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{RDS.files}{A list of RDS files generated with \code{Segment.ASCAT()}.}
  \item{nthread}{Number of thread(s) for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any other parameter for \code{Annotate()} except \code{data}.}
}
\author{
  Bastien Job
}

