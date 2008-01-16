\name{qaReport}
\alias{qaReport}

\title{ Create HTML report using one or several QA process function(s) }
\description{
  
  This function combines all graphical output of multiple QA process
  functions for one \code{\link[flowCore]{flowSet}} in a single
  hyperlinked HTML document.

}
\usage{
qaReport(set, qaFunctions, outdir = "./qaReport", argLists, grouping = NULL)
}

\arguments{
  \item{set}{ A \code{\link[flowCore]{flowSet}} }
  \item{qaFunctions}{ A character vector of the names of QA process
  functions to be used }
  \item{outdir}{ The directory to which the HTML report is to be
  saved. }
  \item{argLists}{ lists of argument lists for each of the QA process
  functions specified via \code{qaFunctions} }
  \item{grouping}{ A character scalar indicating a variable in the
  \code{\link[flowCore]{flowSet}}'s \code{phenoData} that is used as a
  grouping factor in the output. }
}
\details{
  This is a simple convenience function to produce HTML QA reports
  for a single \code{\link[flowCore]{flowSet}} given a list of QA
  process functions. For more fine-grained control use function 
  \code{\link{writeQAReport}} directly.

  An entry point to the output of this function can be found at
  \code{outdir/index.html}.
}
\value{
 The function is called for it's side effects
}

\author{ Florian Hahne }

\seealso{ \code{\link{qaProcess.marginevents}}, \code{\link{writeQAReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}}}
% \examples{

% }

\keyword{ IO }
\keyword{ dynamic }
