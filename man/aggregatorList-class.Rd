\name{aggregatorList-class}
\docType{class}
\alias{aggregatorList-class}
\alias{aggregatorList}
\alias{initialize,aggregatorList-method}
\alias{show,aggregatorList-method}

\title{Class "aggregatorList" }
\description{ A list of qaAggregators }
\section{Objects from the Class}{
  Objects should be created by the constructor:

  \code{aggregatorList(\dots), where \code{\dots}} are objects
  inheriting form class \code{\linkS4class{qaAggregator}} or a list of
  such objects.
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, the list data }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "aggregatorList")}: constructor }
    \item{show}{\code{signature(object = "aggregatorList")}: print
	 object details }
	 }
}
    
\author{ Florian Hahne }


\seealso{ \code{\linkS4class{qaGraph}}, \code{\link{writeQAReport}},
    \code{\linkS4class{qaProcess}}, \code{\linkS4class{qaAggregator}}
}
\examples{
showClass("aggregatorList")
}
\keyword{classes}
  