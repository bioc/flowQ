\name{qaProcess.cellnumber}
\alias{qaProcess.cellnumber}

\title{  Create QA process of type 'cellnumber' }

\description{
  
   This function takes a \code{\link[flowCore]{flowSet}} as input and
  creates all necessary output for a 'cellnumber' type QA process. Objects
  created by this function can be laid out as HTML using
  \code{\link{writeQAReport}}.
  
}

\usage{

qaProcess.cellnumber(set, grouping=NULL, outdir, cFactor=0.5,
   name="cell number", sum.dimensions=c(7,7), pdf=TRUE, \dots)

}

\arguments{
  
  \item{set}{ A \code{\link[flowCore]{flowSet}} }

  \item{grouping}{ A character vector defining one of the variables in
  the phenoData of \code{set} used as a grouping variable. If this
  argument is uses, comparisons will be made within groups rather than
  across all samples. }
  
  \item{outdir}{ The directory to which the graphical output is to be
    saved. If multiple QA processes are to be combined, make sure to use
    the same directory every time. }

  \item{cFactor}{ The threshold at which the QA process is considered to
    be failed. The factor of standard deviations away from the average
    number of cells per sample }
  
  \item{name}{The name of the process used for the headings in the HTML
    output}
  
   \item{sum.dimensions}{The pdf dimensions used for the
     summary. }
   
   \item{pdf}{Create vectorized versions of images for this quality
     process.}
   
   \item{\dots}{Further arguments.}
   
 }
 \details{
   
  QA processes of type 'cellnumber' detect aberations in the number of
  cells analyzed per frame.
  
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
  \code{\link{qaProcess.timeflow}}, \code{\link{qaProcess.timeline}}

}


\examples{
\dontrun{
data(GvHD)
dest <- tempdir()
qp <- qaProcess.cellnumber(GvHD, outdir=dest, cFactor=1)
qp
}

}

\keyword{ dynamic }


