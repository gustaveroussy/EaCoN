\name{BINpack.Maker}
\alias{BINpack.Maker}
\title{Generates a "BINpack" containing binned data and GC tracks from a capture BED file for WES data.}
\usage{
BINpack.Maker(bed.file = NULL, bin.size = 50, genome.pkg =
"BSgenome.Hsapiens.UCSC.hg19", extend.multi = c(0, 50, 100, 200, 400, 800,
1600, 3200, 6400), blocksize = 1E+04, nthread = 1, out.dir = getwd(),
return.data = FALSE)
}
\arguments{
\item{bed.file}{BED File corresponding to the captured genomic regions in the WES experiment (usually given by the manufacturer.}
\item{bin.size}{Size of the bins (in nucleotides) to reduce WES data.}
\item{genome.pkg}{Name of the BSgenome package corresponding to the genome used to map the BAM files.}
\item{extend.multi}{Vector of integers corresponding to lengths to extend (on both sides) the BED regions to compute GC\%.}
\item{blocksize}{Size (in lines) of the block of BED regions to compute GC\% in multithreading mode (to limit RAM consumption).}
\item{nthread}{Number of thread(s).}
\item{out.dir}{Output directory.}
\item{return.data}{Return the results.}
}
\description{
This function performs the binning of WES data from a BED of the captured regions and compute corresponding GC\% with optional elongation.
}
\value{
A data pack (.rda) corresponding to the binned BED of captured regions with its precomputed GC\%.
If \code{return.data} is set to \code{TRUE}, its content is also returned.
}
\details{
Before BED slicing in bins of specified size and GC\% computations, a cleaning of the BED is automaticaly performed, to fix some possible BED errors (inversed coordinates, duplicated or overlapping intervals, etc).
}
\author{
Bastien Job
}

