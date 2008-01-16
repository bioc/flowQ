\name{rangeAggregator-class}
\docType{class}
\alias{rangeAggregator-class}
\alias{rangeAggregator}
\alias{show,rangeAggregator-method}
\alias{writeLines,rangeAggregator,file,missing-method}

\title{Class "rangeAggregator"  }
\description{ Abstraction of a range type of aggregator where possible
  states are within certain ranges (e.g. percentages) }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rangeAggregator", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{min}:}{Object of class \code{"numeric"}, the range minimum }
    \item{\code{max}:}{Object of class \code{"numeric"}, the range maximum }
    \item{\code{x}:}{Object of class \code{"numeric"}, the value within
      the range }
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{"\linkS4class{numericAggregator}"}, directly.
Class \code{"\linkS4class{qaAggregator}"}, by class "numericAggregator", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "rangeAggregator")}: print
  object details }
    \item{writeLines}{\code{signature(text = "rangeAggregator", con = "file", sep = "missing")}:  write to HTML file connection }
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
  \code{\linkS4class{binaryAggregator}},
}
\examples{
showClass("rangeAggregator")
}
\keyword{classes}


 