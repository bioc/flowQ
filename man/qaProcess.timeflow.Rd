\name{qaProcess.timeflow}
\alias{qaProcess.timeflow}

\title{  Create QA process of type 'timeflow' }
\description{
   This function takes a \code{\link[flowCore]{flowSet}} as input and
  creates all necessary output for a 'timeflow' type QA
  process. Objects created by this function can be laid out as HTML
  using \code{\link{writeQAReport}}.
}
\usage{
qaProcess.timeflow(set, outdir, cutoff=2, name="time flow",
sum.dimensions=c(7,7), det.dimensions=c(7,7))
}

\arguments{
  
  \item{set}{ A \code{\link[flowCore]{flowSet}} }
  
  \item{outdir}{ The directory to which the graphical output is to be
    saved. If multiple QA processes are to be combined, make sure to use
    the same directory every time. }

  \item{cutoff}{ The threshold at which the QA process is considered to
    be failed. An absolute value in the timeline deviation score as
    computed by the \code{\link[flowViz]{timeLinePlot}} function}
  
  \item{name}{The name of the process used for the headings in the HTML
    output}
  
   \item{sum.dimensions, det.dimensions}{The pdf dimensions used for the
     summary and the detailed plots. }

 }
 \details{
   
  QA processes of type 'timeflow' detect disturbances in the flow of
  cells over time.
  
  For more details on how to layout \code{\link{qaProcess}} objects to
  HTML, see \code{\link{writeQAReport}} and \code{\link{qaReport}}.
  
}
\value{
  
  An object of class \code{\link{qaProcess}}.
  
}

\author{ Florian Hahne }



\seealso{

  \code{\link{writeQAReport}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.marginevents}},
  \code{\link{qaProcess.timeline}}, \code{\link{qaProcess.cellnumber}}

}


\examples{

\dontrun{
data(GvHD)
dest <- tempdir()
qp <- qaProcess.timeflow(GvHD, outdir=dest, cutoff=1)
qp
}

}

\keyword{ dynamic }


