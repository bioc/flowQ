\name{qaAggregator-class}
\docType{class}
\alias{qaAggregator-class}
\alias{qaAggregator}

\title{ Abstraction of the possible outcomes of a QA process }
\description{ Virtual parent class for different types of QA aggregators}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Methods}{
No methods defined with class "qaAggregator" in the signature.
}

\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{binaryAggregator}},
  \code{\linkS4class{discreteAggregator}},
  \code{\linkS4class{factorAggregator}},
  \code{\linkS4class{numericAggregator}},
  \code{\linkS4class{stringAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("qaAggregator")
}
\keyword{classes}
