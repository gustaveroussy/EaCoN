\name{EaCoN.Segment}
\alias{EaCoN.Segment}
\title{L2R and BAF joint segmentation using ASCAT.}
\usage{
  EaCoN.Segment(data = NULL, mingap = 5E+06, segmentLength = 10, homoCut = .05,
  mc.G = 3:4, nfactor = 2, bafbin.size = 1E+07, BAF.rescale = FALSE,
  BAF.filter = .90,  BAF.cutter = 0, smooth.k = NULL, ASPCF.pen = 100,
  recenter = "l2r.centeredpeak", calling.method = "mad",  nrf = .5,
  SER.pen = 1-(1E-04),  out.dir = getwd(), return.data = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.WES.Normalize()} (for WES data), \code{EaCoN.CS.Process()} (for Affymetrix CytoScan data), \code{EaCoN.OS.Process()} (for Affymetrix OncoScan data), or \code{EaCoN.SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{mingap}{Minimum size (in nt) of a genomic gap (mainly, centromere), to allow ASCAT to treat flanking parts independently.}
  \item{segmentLength}{Minimal segment size for germline prediction and (optional) PELT segmentation.}
  \item{homoCut}{Threshold value to consider a BAF segment as homozygote.}
  \item{mc.G}{Number of BAF clusters to evaluate for germline prediction (see parameter \code{G} in \code{?mclust::Mclust()}.}
  \item{nfactor}{A multiplier to identify homozygous regions (the lower, the more sensitive).}
  \item{bafbin.size}{Size (in nucleotide) of the window for the BAF germline prediction.}
  \item{BAF.rescale}{Perform a re-scaling of BAF profile to have homozygous regions centered at O and 1. Recommended for CytoScan arrays.}
  \item{BAF.filter}{Perform a filtering of BAF noise, keeping a fraction of heterozygous BAF markers, represented by this parameter. Can be avoided if set to \code{NULL}.}
  \item{BAF.cutter}{A threshold to discard homozygous BAF probes when performing germline prediction. increase this value (slightly, up to 0.05) in case of noisy (smeary / blurry) BAF profile.}
  \item{smooth.k}{The half-window size for data smoothing by winsorization. If \code{NULL}, no winsorization performed.}
  \item{ASPCF.pen}{The penalty for \code{ASPCF} segmentation from ASCAT (see \code{?ASCAT::ascat.aspcf()}).}
  \item{recenter}{Method to perform re-centering of the L2R profile. Can be : \code{"l2r.mainpeak"} (use the highest peak of the L2R density), \code{"l2r.centeredpeak"} (use the most-centered peak), \code{"l2r.median"} (just center to the median of the L2R distribution. Not recommended), any signed numeric (to force centering to a desired level), or \code{NULL} (do not recenter).}
  \item{calling.method}{Method to perform the calling of aberrations (attribution of the Loss / Normal / Gain status). Can be \code{"mad"} (build threshold based on the L2R spread noise, altered by the \code{nrf} coefficient) or \code{"density"} (automatically build thresholds based on the L2R density).}
  \item{nrf}{Coefficient to multiply the L2R spread noise, to build a calling threshold when \code{calling.method} is set to \code{"mad"}. Not used when \code{calling.method} set to \code{"density"}}
  \item{SER.pen}{Penalty for the second-step PELT segmentation of the L2R profile, to retrieve small events ignored by ASCAT. Can be a numeric in ]0..1[ (but typicaly very, very close to 1), or any of : \code{"MBIC"}, \code{"BIC"},\code{"BIC0"}, \code{"AIC"}, \code{"AIC0"}, \code{"Hannan-Quinn"}, or \code{NULL} (do not perform a second-step PELT segmentation at all).}
  \item{out.dir}{Output directory.}
\item{return.data}{Return the results.}
}
\description{
  This function performs the BAF germline prediction, joint L2R and BAF segmentation, centering and calling of WES, CytoScan, OncoScan or SNP6 data.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables. If \code{ASCN} set to \code{TRUE}, several folders containing the different ASCN (allele-specific copy number) analyses performed for each value in \code{gammaRange}. The directory corresponding to the best evaluated gamma will also have an "_optimal" suffix.
}
\details{
  Depending on the data source, different parameter values are recommended :
  \itemize{
    \item{for WES : TO DEFINE.}
    \item{for CytoScan array : TO DEFINE.}
    \item{for OncoScan array : TO DEFINE.}
    \item{for OncoScan array : TO DEFINE.}
  }
  \code{SER.pen} : The ASPCF algorithm from ASCAT::ascat.aspcf commonly ignores (very) small events (especially in loss) when not supported by the BAF profile). This option allows to detect such events using a PELT segmentation. See \code{?changepoint::cpt.mean()} to see which penalties are available. Please set it to NULL to deactivate this step.
  Please note :
  \itemize{
    \item{1) Depending on noise and complexity, ASCAT will not always be able to perform ASCN.}
    \item{2) Default parameters correspond to recommended ones for CytoScanHD microarray.}
  }
}
\author{
Bastien Job
}

