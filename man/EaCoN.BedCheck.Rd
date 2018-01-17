\name{EaCoN.BedCheck}
\alias{EaCoN.BedCheck}
\title{Performs checks un cleanups on a WES capture bed file}
\usage{
EaCoN.BedCheck(bed.file = NULL, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19",
out.dir = getwd(), return.data = FALSE)
}
\arguments{
\item{bed.file}{BED File corresponding to the captured genomic regions in the WES experiment (usually given by the manufacturer.}
\item{genome.pkg}{Name of the BSgenome package corresponding to the genome used to map the BAM files.}
  \item{out.dir}{Output directory.}
\item{return.data}{Return the checked / cleaned BED as a data.frame (and still write it as a BED file).}
}
\description{
This function cleans (when necessary) a BED file by detecting and merging overlapped regions and sorting by chr > start > end.
}
\value{
Checked / cleaned BED is written in \code{out.dir}, and returned as a \code{data.frame} if \code{return.data} is set to \code{TRUE}.
}
\details{
Only genomes available in the \code{chomosomes} data package are compatible with this function.
}
\author{
Bastien Job
}

