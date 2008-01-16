\name{discreteAggregator-class}
\docType{class}
\alias{discreteAggregator-class}
\alias{show,discreteAggregator-method}
\alias{writeLines,discreteAggregator,file,missing-method}

\title{Class "discreteAggregator" }
\description{ Abstraction of a binary type of aggregator with possible
  states "passes" and "not passed" and "warning"}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("discreteAggregator", ...)}.
}

\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"factor"} One in \code{0}
	 (not passes),
	 \code{1} (passed) or \code{2} (warning)}
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements (not eveluated) }
  }
}
\section{Extends}{
Class \code{"\linkS4class{qaAggregator}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "discreteAggregator")}: print
  object details }
    \item{writeLines}{\code{signature(text = "discreteAggregator", con = "file", sep = "missing")}:  write to HTML file connection }
	 }
 }

 
\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{binaryAggregator}},
  \code{\linkS4class{factorAggregator}},
  \code{\linkS4class{numericAggregator}},
  \code{\linkS4class{stringAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("discreteAggregator")
}
\keyword{classes}
 