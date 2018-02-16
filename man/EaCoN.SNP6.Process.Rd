\name{EaCoN.SNP6.Process}
\alias{EaCoN.SNP6.Process}
\title{Normalize GenomeWideSNP_6 arrays using Affymetrix Power Tools (APT)}
\description{
This function performs the normalization of a CEL file of Affymetrix GenomeWideSNP_6 microarray using the APT binary (apt-copynumber-workflow).
}
\usage{
EaCoN.SNP6.Process(CEL = NULL, samplename = NULL, l2r.level = "normal",
renormalize = TRUE, renorm.rda = NULL, out.dir = getwd(),
oschp.keep = TRUE, force.OS = NULL, apt.version = "1.20.0",
apt.build = "na35.r1", genome.pkg = 
"BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE)
}
\arguments{
\item{CEL}{A CEL file.}
\item{samplename}{Name of the analyzed sample.}
\item{l2r.level}{Type of log2(ratio) to extract : \code{'weighted'} or \code{'normal'}.}
\item{renormalize}{Perform a GC-content re-normalization after APT processing.}
\item{renorm.rda}{If \code{renormaize} set to \code{TRUE}, specifies the data package containing GC-content tracks. If \code{NULL}, a default one (included in the package) will be used.}
\item{out.dir}{The path to output results.}
\item{oschp.keep}{Keep the intermediate OSCHP (HDF5) file.}
\item{force.OS}{Specifiy the OS type of the execution machine (see Details).}
\item{apt.version}{Version of the APT tool \code{apt-copynumber-cyto}.}
\item{apt.build}{Affymetrix genome build version.}
\item{genome.pkg}{The name of an installed \code{BSgenome} package. Please verify that provided package is compatible with the provided \code{apt.build} genome version.}
\item{return.data}{Return data when finished.}
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

