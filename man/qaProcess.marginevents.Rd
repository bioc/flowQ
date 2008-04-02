\name{qaProcess.marginevents}
\alias{qaProcess.marginevents}

\title{ Create QA process of type 'marginevents' }
\description{
  This function takes a \code{\link[flowCore]{flowSet}} as input and
  creates all necessary output for a 'marginevents' type QA
  process. Objects created by this function can be laid out as HTML
  using \code{\link{writeQAReport}}.
}
\usage{
qaProcess.marginevents(set, channels = NULL, outdir, cFactor = 3)
}

\arguments{
  \item{set}{ A \code{\link[flowCore]{flowSet}} }
  \item{channels}{ A character vector of channel names for which margin
  events are to be recorded }
  \item{outdir}{ The directory to which the graphical output is to be
  saved. If multiple QA processes are to be combined, make sure to use
  the same directory every time. }
  \item{cFactor}{ The threshold at which the QA process is considered to
    have failed. This is the fold change of margin event percentages
  compared to the median percentage of events for the respective channel.}
}
\details{
  QA processes of type 'marginevents' record the number of events that
  fall on the margins of the measurement range for each
  channel. Unproportionally high numbers of such events can indicate
  problems with the instrument settings.
  
  For more details on how to layout \code{\link{qaProcess}} objects to
  HTML, see \code{\link{writeQAReport}} and \code{\link{qaReport}}.
}
\value{
  An object of class \code{\link{qaProcess}}.
}

\author{ Florian Hahne }
\note{ This function is still experimental }


\seealso{ \code{\link{writeQAReport}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}}}

 \examples{
\dontrun{
data(GvHD)
dest <- tempdir()
qp <- qaProcess.marginevents(GvHD, channels=c("FL1-H", "FL2-H"),
  outdir=dest)
qp
}

 }

\keyword{ dynamic }

