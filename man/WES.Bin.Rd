\name{WES.Bin}
\alias{WES.Bin}
\title{Computes the read depth of a pair of BAM files and bins it.}
\usage{
WES.Bin(testBAM = NULL, refBAM = NULL, BINpack = NULL, samplename = "SAMPLE",
Q = 20, nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(),
return.data  = FALSE, write.data = TRUE, plot = TRUE, force = FALSE)
}
\arguments{
\item{testBAM}{BAM file corresponding to the tested sample (typically a patient tumor sample).}
\item{refBAM}{BAM file corresponding to the reference sample (typically a patient normal sample).}
\item{BINpack}{A data pack generated with \code{EaCoN.BINpack.Maker()} containing to a binned version of the capture BED used to generate the BAM files, and corresponding precomputed GC\% tracks.}
\item{samplename}{The name of the sample.}
\item{Q}{Minimum baseQ/BAQ quality score as defined by \code{samtools mpileup}.}
\item{nsubthread}{Number of subthreads for multithreading.}
\item{cluster.type}{Type of subthreading cluster (see \code{?parallel::makeCluster()}).}
\item{out.dir}{Output directory.}
\item{return.data}{Return the results.}
\item{write.data}{Write the results as a RDS file.}
\item{plot}{Perform plots.}
\item{force}{If some result files already exist, delete them before runing.}
}
\description{
This function computes the read depth of a pair of BAM files (typicaly tumor and normal) and bins it thanks to a "BINpack" generated with \code{BINpack.Maker()}.
}
\value{
A RDS object containing binned read depths and SNP depths.
If \code{return.data} is set to \code{TRUE}, its content is also returned.
}
\author{
Bastien Job
}

