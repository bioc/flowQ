\name{qaProcess.timeline}
\alias{qaProcess.timeline}

\title{  Create QA process of type 'timeline' }
\description{
   This function takes a \code{\link[flowCore]{flowSet}} as input and
  creates all necessary output for a 'timeline' type QA
  process. Objects created by this function can be laid out as HTML
  using \code{\link{writeQAReport}}.
}
\usage{
qaProcess.timeline(set, channel, outdir, cutoff = 0.1)
}

\arguments{
  \item{set}{ A \code{\link[flowCore]{flowSet}} }
  \item{channel}{ A character vector of channel names for which margin
  events are to be recorded }
  \item{outdir}{ The directory to which the graphical output is to be
  saved. If multiple QA processes are to be combined, make sure to use
  the same directory every time. }
  \item{cutoff}{ The threshold at which the QA process is considered to
    be failed. An absolute value in the timeline deviation score as
  computed by the \code{\link[flowViz]{timeLinePlot}} function}
}
\details{
  QA processes of type 'timeline' detect unusal patterns in the data
  acquisition over time.
  
  For more details on how to layout \code{\link{qaProcess}} objects to
  HTML, see \code{\link{writeQAReport}} and \code{\link{qaReport}}.
}
\value{
  An object of class \code{\link{qaProcess}}.
}

\author{ Florian Hahne }
\note{ This function is still experimental }


\seealso{ \code{\link{writeQAReport}}, \code{\link{qaReport}},
  \code{\link{qaProcess}},  \code{\link{qaProcess.marginevents}}}


\examples{
\dontrun{
data(GvHD)
dest <- tempdir()
qp1 <- qaProcess.timeline(GvHD, channel="FL1-H", outdir=dest)
qp1
}

}

\keyword{ dynamic }


