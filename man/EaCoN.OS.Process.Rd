\name{EaCoN.OS.Process}
\alias{EaCoN.OS.Process}
\title{Normalize OncoScan / OncoScan_CNV arrays using Affymetrix Power Tools (APT)}
\description{
This function performs the normalization of a CEL file pair of Affymetrix OncoScan or OncoScan_CNV microarray using the APT binary (apt-copynumber-onco-ssa).
}
\usage{
EaCoN.OS.Process(ATChannelCel = NULL, GCChannelCel = NULL,
samplename = NULL, dual.norm = TRUE, l2r.level = "weighted",
wave.renorm = TRUE, wave.renorm.rda = NULL, gc.renorm = TRUE,
gc.renorm.rda = NULL, out.dir = getwd(), oschp.keep = TRUE,
force.OS = NULL, apt.version = "2.4.0", apt.build = "na33.r2",
genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE)
}
\arguments{
\item{ATChannelCel}{The first CEL file (version A).}
\item{GCChannelCel}{The second CEL file (version C).}
\item{samplename}{Name of the analyzed sample.}
\item{dual.norm}{Perform a dual-channel normalization (recommended for a better wave-effect bias normalization).}
\item{l2r.level}{Type of log2(ratio) to extract : \code{'weighted'} or \code{'normal'}.}
\item{wave.renorm}{Perform a 'Wave-effect' re-normalization after APT processing.}
\item{wave.renorm.rda}{If \code{wave.renorm} set to \code{TRUE}, specifies the data package containing GC-content tracks. If \code{NULL}, a default one (included in the \code{affy.CN.norm.data} package) will be used.}
\item{gc.renorm}{Perform a 'GC-content' re-normalization after APT processing.}
\item{gc.renorm.rda}{If \code{gc.renorm} set to \code{TRUE}, specifies the data package containing GC-content tracks. If \code{NULL}, a default one (included in the \code{affy.CN.norm.data} package) will be used.}
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

