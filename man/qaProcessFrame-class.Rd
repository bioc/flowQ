\name{qaProcessFrame-class}
\docType{class}
\alias{qaProcessFrame-class}
\alias{qaProcessFrame}
\alias{initialize,qaProcessFrame-method}
\alias{show,qaProcessFrame-method}

\title{Class "qaProcessFrame" }
\description{ Abstraction of subitems within a qQA process}
\section{Objects from the Class}{
  Objects should be created using the constructor:

  \code{qaProcessFrame(frameID, summaryAggregator, summaryGraph,
    frameAggregators, frameGraphs)} where \code{frameID} is the ID of the
    flowFrame the process is linked to, \code{summaryAggregator} is an
    object inheriting from class  \code{\linkS4class{qaAggregator}}
    which summarizes the outcome, \code{summaryGraph} is an object of
    class \code{\linkS4class{qaGraph}}, \code{frameAggregators} is an
    object of class \code{\linkS4class{aggregatorList}} and
    \code{frameGraphs} is an object of class \code{\linkS4class{qaGraphList}}
}
\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"character"}, a unique ID of
    the object}
    \item{\code{frameID}:}{Object of class \code{"character"}, ID of the
    flowFrame the process is linked to }
    \item{\code{summaryAggregator}:}{Object of class
    \code{"qaAggregator"}, an aggregator summarizing the output of the process }
    \item{\code{summaryGraph}:}{Object of class \code{"qaGraph"}, a
    graphical summary of the process }
    \item{\code{frameAggregators}:}{Object of class
    \code{"aggregatorList"}  a list of aggregators for the subprocesses}
    \item{\code{frameGraphs}:}{Object of class \code{"qaGraphList"} a
      list of graphical summaries for the subprocesses }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "qaProcessFrame")}: constructor }
	 }
       }
       

\author{ Florian Hahne }


\seealso{ \code{\linkS4class{qaGraphList}}, \code{\link{writeQAReport}}, \code{\linkS4class{qaProcess}}
}
\examples{
showClass("qaProcessFrame")
}
\keyword{classes}
