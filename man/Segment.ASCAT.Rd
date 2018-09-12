\name{Segment.ASCAT}
\alias{Segment.ASCAT}
\title{L2R and BAF joint segmentation using ASCAT.}
\usage{
  Segment.ASCAT(data = NULL, mingap = 5E+06, smooth.k = NULL, BAF.filter = .75,
  homoCut = .05, ASPCF.pen = 50, recenter = "l2r.centeredpeak",
  calling.method = "mad", nrf = 0.5, SER.pen = 40, out.dir = getwd(),
  return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.WES.Normalize()} (for WES data), \code{EaCoN.CS.Process()} (for Affymetrix CytoScan data), \code{EaCoN.OS.Process()} (for Affymetrix OncoScan data), or \code{EaCoN.SNP6.Process()} (for Affymetrix GenomeWideSNP_6 data).}
  \item{mingap}{Minimum size (in nt) of a genomic gap (mainly, centromere), to allow ASCAT to treat flanking parts independently.}
  \item{smooth.k}{The half-window size for data smoothing by winsorization. If \code{NULL}, no winsorization performed.}
  \item{BAF.filter}{Perform a filtering of BAF noise, keeping a fraction of heterozygous BAF markers, represented by this parameter.}
  \item{homoCut}{Threshold value to consider a BAF segment as homozygote.}
  \item{ASPCF.pen}{The penalty for \code{ASPCF} segmentation from ASCAT (see \code{?ASCAT::ascat.aspcf()}).}
  \item{recenter}{Method to perform re-centering of the L2R profile. Can be : \code{"l2r.mainpeak"} (use the highest peak of the L2R density), \code{"l2r.centeredpeak"} (use the most-centered peak), \code{"l2r.median"} (just center to the median of the L2R distribution. Not recommended), any signed numeric (to force centering to a desired level), or \code{NULL} (do not recenter).}
  \item{calling.method}{Method to perform the calling of aberrations (attribution of the Loss / Normal / Gain status). Can be \code{"mad"} (build threshold based on the L2R spread noise, altered by the \code{nrf} coefficient) or \code{"density"} (automatically build thresholds based on the L2R density).}
  \item{nrf}{Coefficient to multiply the L2R spread noise, to build a calling threshold when \code{calling.method} is set to \code{"mad"}. Not used when \code{calling.method} set to \code{"density"}}
  \item{SER.pen}{Penalty for the small events rescue steop (a second-step segmentation of the L2R profile using PELT, to retrieve small events ignored by ASCAT). Can be a positive integer, or \code{NULL} (do not perform a second-step PELT segmentation at all).}
  \item{out.dir}{Output directory.}
  \item{return.data}{Return the results.}
  \item{write.data}{Write the results as a RDS file.}
  \item{plot}{Perform plots.}
  \item{force}{If some result files already exist, delete them before runing.}
}
\description{
  This function performs the joint L2R and BAF segmentation using ASCAT PCF-based segmenter, centering and calling of WES, CytoScan, OncoScan or SNP6 data.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables.
}
\details{
  Depending on the data source, different parameter values are recommended :
  \itemize{
    \item{for WES :             SER.pen = 2}
    \item{for CytoScan array :  SER.pen = 20, smooth.k = 5}
    \item{for OncoScan array :  SER.pen = 40}
    \item{for SNP6 array :      SER.pen = 60, smooth.k = 5, nrf = .25}
  }
  \code{SER.pen} : The ASPCF algorithm from ASCAT::ascat.aspcf commonly ignores (very) small events (especially in loss) when not supported by the BAF profile). This option allows to detect such events using a PELT segmentation. See \code{?changepoint::cpt.mean()}.
}
\author{
Bastien Job
}

