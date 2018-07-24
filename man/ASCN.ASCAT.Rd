\name{ASCN.ASCAT}
\alias{ASCN.ASCAT}
\title{Allele-Specific Copy Number estimation using ASCAT.}
\usage{
  ASCN.ASCAT(data = NULL, gammaRange = c(.35,.95), nsubthread = 1,
  cluster.type = "PSOCK", out.dir = getwd(), force = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{Segment.ASCAT()} (or corresponding wrappers \code{Segment.ASCAT.ff()} and \code{Segment.ASCAT.ff.Batch()}.}
  \item{gammaRange}{The range of gamma values to assess the modelling of ASCN (see \code{?ASCAT::ascat.runAscat}).}
  \item{nsubthread}{Number of subthreads for multithreading.}
  \item{cluster.type}{Type of subthreading cluster (see \code{?parallel::makeCluster()}).}
  \item{out.dir}{Output directory.}
  \item{force}{If some result files already exist, delete them before runing.}
}
\description{
  This function performs the estimation of allele-specific copy number using ASCAT.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables in several folders containing the different ASCN (allele-specific copy number) analyses performed for each value in \code{gammaRange}. The optimal fit among the generated ones will be output in a separate text file.
}
\details{
  Please note that depending on noise and complexity of the L2R and BAF profiles, ASCAT will not always be able to perform ASCN.
}
\author{
  Bastien Job
}

