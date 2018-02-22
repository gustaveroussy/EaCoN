\name{EaCoN.Annotate}
\alias{EaCoN.Annotate}
\title{Annotation of segmented regions and generation of an HTML QC report}
\usage{
  EaCoN.Annotate(data = NULL, refGene.table = NULL, targets.table = NULL,
  report = TRUE, author.name, out.dir = getwd())
}
\arguments{
  \item{data}{Data object contained in a RDS file generated with \code{EaCoN.Segment()}..}
  \item{refGene.table}{A table containing gene symbols. If \code{NULL} and \code{genome} set to "hg18", "hg19" or "hg38", a default list will be used.}
  \item{targets.table}{A table containing a selected list of gene symbols (typically, a list of drug-actionable genes). If \code{NULL} and \code{genome} set to "hg19", a default list of 475 cancer genes will be used.}
  \item{report}{Generate the HTML report.}
  \item{author.name}{Name of the author/lab/facility to display in the HTML report.}
  \item{out.dir}{The output directory.}
}
\description{
  This function performs the genomic annotation of segmented regions, and builds an HTML interactive QC report.
}
\author{
  Bastien Job
}

