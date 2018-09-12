\name{Segment.FACETS}
\alias{Segment.FACETS}
\title{L2R and BAF joint segmentation using FACETS.}
\usage{
  Segment.FACETS(data = NULL, smooth.k = NULL, BAF.filter = .9, homoCut = .05,
  FACETS.pen = 150, recenter = "l2r.centeredpeak", calling.method = "mad",
  nrf = .5, SER.pen = 2, out.dir = getwd(), return.data = FALSE,
  write.data = TRUE, plot = TRUE, force = FALSE)
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.WES.Normalize()}, as this segmenter is only compatible with HTS results.}
  \item{smooth.k}{The half-window size for data smoothing by winsorization. If \code{NULL}, no winsorization performed.}
  \item{BAF.filter}{Perform a filtering of BAF noise, keeping a fraction of heterozygous BAF markers, represented by this parameter.}
  \item{homoCut}{Threshold value to consider a BAF segment as homozygote.}
  \item{FACETS.pen}{The penalty for segmentation from FACETS (see \code{?facets::preProcSample()} and \code{?facets::procSample()}).}
  \item{recenter}{Method to perform re-centering of the L2R profile. Can be : \code{"l2r.mainpeak"} (use the highest peak of the L2R density), \code{"l2r.centeredpeak"} (use the most-centered peak), \code{"l2r.median"} (just center to the median of the L2R distribution. Not recommended), any signed numeric (to force centering to a desired level), or \code{NULL} (do not recenter).}
  \item{calling.method}{Method to perform the calling of aberrations (attribution of the Loss / Normal / Gain status). Can be \code{"mad"} (build threshold based on the L2R spread noise, altered by the \code{nrf} coefficient) or \code{"density"} (automatically build thresholds based on the L2R density).}
  \item{nrf}{Coefficient to multiply the L2R spread noise, to build a calling threshold when \code{calling.method} is set to \code{"mad"}. Not used when \code{calling.method} set to \code{"density"}}
  \item{SER.pen}{Penalty for the small events rescue steop (a second-step segmentation of the L2R profile using PELT, to retrieve small events ignored by FACETS). Can be a positive integer, or \code{NULL} (does not perform a second-step PELT segmentation at all).}
  \item{out.dir}{Output directory.}
  \item{return.data}{Return the results.}
  \item{write.data}{Write the results as a RDS file.}
  \item{plot}{Perform plots.}
  \item{force}{If some result files already exist, delete them before runing.}
}
\description{
  This function performs the joint L2R and BAF segmentation using FACETS CBS-based bivariate segmenter, centering and calling of WES data.
}
\value{
  Several RDS objects, PNG plots, CBS files and data tables.
}
\author{
Bastien Job
}

