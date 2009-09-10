\name{binaryAggregator-class}
\docType{class}
\alias{binaryAggregator-class}
\alias{binaryAggregator}
\alias{show,binaryAggregator-method}
\alias{writeLines,binaryAggregator,file,missing,missing-method}

\title{Class "binaryAggregator" }
\description{ Abstraction of a binary type of aggregator with possible
  states "passed" and "not passed"}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("binaryAggregator", ...)},
or using the constructor \code{binaryAggregator(passed)}, where
  \code{passed} is a logical scalar.
}
\section{Slots}{
	 \describe{
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{"\linkS4class{qaAggregator}"}, directly.
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "binaryAggregator")}: print
  object details }
    \item{writeLines}{\code{signature(text = "binaryAggregator", con =
  "file", sep = "missing", useBytes="missing")}: write to HTML file connection }
	 }
}

\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{discreteAggregator}},
  \code{\linkS4class{factorAggregator}},
  \code{\linkS4class{numericAggregator}},
  \code{\linkS4class{stringAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("binaryAggregator")
}
\keyword{classes}

