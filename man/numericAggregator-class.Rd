\name{numericAggregator-class}
\docType{class}
\alias{numericAggregator-class}
\alias{numericAggregator}
\alias{show,numericAggregator-method}
\alias{writeLines,numericAggregator,file,missing-method}

\title{Class "numericAggregator" }
\description{ Abstraction of a numeric type of aggregator for which possible
  states are coded by a numeric value}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("numericAggregator", ...)}.

}
\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} coding the outcome
    state}
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{"\linkS4class{qaAggregator}"}, directly.
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "numericAggregator")}: print
  object details}
    \item{writeLines}{\code{signature(text = "numericAggregator", con = "file", sep = "missing")}: write to HTML file connection }
	 }
       }

  
\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{discreteAggregator}},
  \code{\linkS4class{factorAggregator}},
  \code{\linkS4class{binaryAggregator}},
  \code{\linkS4class{stringAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("numericAggregator")
}
\keyword{classes}

    