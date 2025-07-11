\name{GC.data}
\docType{data}
\alias{GC.data}
\title{Dataset containing GC-content data}
\description{
  This data package exists in two versions (for the OncoScan and CytoScan HD designs), and contains several tracks corresponding to the GC percent of regions surrounding the design probes position.
  The different tracks correspond to different sizes of elongation around the probes position, namely 2^(0:13)*100 bp.
}
\usage{GC.data}
\format{A dataframe with first four column corresponding to the probe information (chr, start, end, name), then columns corresponding to the GC percent values.}
\details{This dataset is available to perform a GC-content re-normalization of APT-prenormalized data.}
\source{Depends on data type. For Affymetrix microarrays, annotation files downloaded from Affymetrix website. For WES data, the (binned) capture BED.}
\author{Bastien Job (bastien.job@inserm.fr)}
\keyword{datasets, object}
