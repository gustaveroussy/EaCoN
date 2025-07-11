\name{EaCoN.set.bitmapType}
\alias{EaCoN.set.bitmapType}
\title{Changes the environment variable bitmapType}
\description{
  This functions allows to change the value of the environment variable \code{bitmapType} to perform thePNG plots.
}
\usage{
  EaCoN.set.bitmapType(type = "cairo")
}
\arguments{
  \item{type}{The type of render (should be one of \code{Xlib}, \code{cairo} or \code{quartz}.}
}
\details{
  This function allows to perform the PNG plots on platforms which do not have the default window ibrary (Xlib) but an alternate one.
}
\author{
  Bastien Job
}

