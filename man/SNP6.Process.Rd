\name{SNP6.Process}
\alias{SNP6.Process}
\title{Normalize GenomeWideSNP_6 arrays using Affymetrix Power Tools (APT)}
\description{
This function performs the normalization of a CEL file of Affymetrix GenomeWideSNP_6 microarray using the APT binary (apt-copynumber-workflow).
}
\usage{
SNP6.Process(CEL = NULL, samplename = NULL, l2r.level = "normal",
gc.renorm = TRUE, gc.rda = NULL, wave.renorm = TRUE, wave.rda = NULL,
mingap = 1E+06, out.dir = getwd(), oschp.keep = FALSE,force.OS = NULL,
apt.version = "1.20.0", apt.build = "na35.r1", genome.pkg =
"BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE, write.data = TRUE,
plot = TRUE, force = FALSE)
}
\arguments{
\item{CEL}{A CEL file.}
\item{samplename}{Name of the analyzed sample.}
\item{l2r.level}{Type of log2(ratio) to extract : \code{'weighted'} or \code{'normal'}.}
\item{gc.renorm}{Perform a 'GC-content' re-normalization after APT processing.}
\item{gc.rda}{If \code{gc.renorm} set to \code{TRUE}, specifies the data package containing GC-content tracks. If \code{NULL}, a default one (included in the \code{affy.CN.norm.data} package) will be used.}
\item{wave.renorm}{Perform a 'Wave-effect' re-normalization after APT processing.}
\item{wave.rda}{If \code{wave.renorm} set to \code{TRUE}, specifies the data package containing Wave tracks. If \code{NULL}, a default one (included in the \code{affy.CN.norm.data} package) will be used.}
\item{mingap}{Minimum genomic size to be considered as an uncovered gap in the design (allows to discard centromeres).}
\item{out.dir}{The path to output results.}
\item{oschp.keep}{Keep the intermediate OSCHP (HDF5) file.}
\item{force.OS}{Specifiy the OS type of the execution machine (see Details).}
\item{apt.version}{Version of the APT tool \code{apt-copynumber-cyto}.}
\item{apt.build}{Affymetrix genome build version.}
\item{genome.pkg}{The name of an installed \code{BSgenome} package. Please verify that provided package is compatible with the provided \code{apt.build} genome version.}
\item{return.data}{Return data when finished.}
\item{write.data}{Write the results as a RDS file.}
\item{plot}{Perform plots.}
\item{force}{If some result files already exist, delete them before runing.}
}
\details{
  \code{force.OS} is \code{NULL} by default, but can take any of these values :
  \itemize{
    \item{linux}
    \item{windows}
    \item{osx}
  }
}
\author{
Bastien Job
}

