\name{WES.Normalize}
\alias{WES.Normalize}
\title{Normalization of binned WES data.}
\usage{
WES.Normalize(data = NULL, BINpack = NULL, gc.renorm = TRUE,
wave.renorm = FALSE, wave.rda = NULL, RD.tot.min = 20, RD.alt.min = 3,
BAF.hetmin = .33, sex.chr = c("chrX", "chrY"), TumorBoost = FALSE,
out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE)
}
\arguments{
\item{data}{Data object contained in a RDS file generated with \code{EaCoN.WES.Bin()}.}
\item{BINpack}{A data pack generated with \code{EaCoN.BINpack.Maker()} containing to a binned version of the capture BED used to generate the BAM files, and corresponding precomputed GC\% tracks.}
\item{gc.renorm}{Perform a 'GC-content' normalization of read depth.}
\item{wave.renorm}{Perform a 'Wave-effect' normalization of read depth.}
\item{wave.rda}{If \code{wave.norm} set to \code{TRUE}, specifies the data package containing Wave tracks.}
\item{RD.tot.min}{Minimum total read depth in the test and REF BAMs to keep a log2(ratio) bin.}
\item{RD.alt.min}{Minimum alternate nucleotide read depth in the test BAM to keep a SNP position.}
\item{BAF.hetmin}{BAF threshold for reference to consider a SNP position as heterozygous.}
\item{sex.chr}{Name of the gonosomes. Please verify that provided name(s) is (are) compatible with chromosome names used in the corresponding BSgenome package (if none, use \code{NULL}).}
\item{TumorBoost}{Perform a normalization of the Test BAF profile using the Ref BAF profile, thanks to the TumorBoost methods (see \code{?aroma.light::normalizeTumorBoost())}).}
\item{out.dir}{Output directory.}
\item{return.data}{Return the results.}
\item{write.data}{Write the results as a RDS file.}
\item{plot}{Perform plots.}
}
\description{
This function performs the normalization of the L2R and BAF signals from a data object obtained with \code{WES.Bin()} (or \code{WES.Bin.Batch()}).
}
\value{
A RDS file that can be used with \code{EaCoN.Segment.ff()}.
}
\author{
Bastien Job
}

