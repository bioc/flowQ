\name{stringAggregator-class}
\docType{class}
\alias{stringAggregator-class}
\alias{stringAggregator}
\alias{show,stringAggregator-method}
\alias{writeLines,stringAggregator,file,missing-method}

\title{Class "stringAggregator"  }

\description{ Abstraction of a string type of aggregator for which possible
  states are indicated through a textual description }

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("stringAggregator", ...)} or using the constructor
\code{stringAggregator(x, passed)}, where \code{x} is a character
scalar, and \code{passed} is a logical scalar.  }

\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"character"} which is a
      textual description of the outcome }
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{"\linkS4class{qaAggregator}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "stringAggregator")}:  print
  object details}
    \item{writeLines}{\code{signature(text = "stringAggregator", con = "file", sep = "missing")}:  write to HTML file connection}
	 }
}

\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{discreteAggregator}},
  \code{\linkS4class{factorAggregator}},
  \code{\linkS4class{numericAggregator}},
  \code{\linkS4class{binaryAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("stringAggregator")
}
\keyword{classes}

