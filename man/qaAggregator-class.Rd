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

\details{ In the context of this package, \code{qaAggregators} are
  objects that hold the outcome of a QA process. Each subclass
  implements its own \code{writeLine} method, which creates the appropriate
  HTML code for a graphical representation of the object.}

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
