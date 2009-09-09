\name{writeQATextReport}
\alias{writeQATextReport}


\title{ Create tab delimited text files  for (lists of) qaProcess objects }

\description{
  
  This function combines all the results of multiple QA processes
  for one or several \code{\link[flowCore]{flowSet}}s in a single
  tab delimited text document.
  
}
\usage{

writeQATextReport(set, processes=NULL, globalProcess=NULL, 
              fileName="textReport.txt")

}

\arguments{
  
  \item{set}{ A \code{\link[flowCore]{flowSet}} or a list of several
    \code{\link[flowCore]{flowSet}}s.}
  
  \item{processes}{ A list of  \code{\link{qaProcess}} objects or, in
    the case of multiple \code{\link[flowCore]{flowSet}}s, a list of
    lists of \code{\link{qaProcess}} objects. See below for further
    details.}

 \item{globalProcess}{ A list of  \code{\link{qaProcess}} objects generated
 by multipanel comparison qaProcess functions, when set is a list of 
  \code{\link[flowCore]{flowSet}}s .}

 
   \item{fileName}{ A character vector for the name of the file. }

}

\details{
  
  Qunatitative or qualitative results are stored  in objects of class 
  \code{\link{qaProcess}}. \code{\link{writeQATextReport}} takes
  lists of such objects and combines their information in a tab delimited
  text document. In the case of multiple panels, a list of
  \code{\link[flowCore]{flowSet}}s can be given to
  \code{\link{writeQAReport}}, and the function expects a list of lists
  of processes, where each process list is specific to one panel.
  
}
\value{
  
  The function is mostly called for it's side effects, that is, the
  generation of a tab delimited text file. .
  
}


\author{ Nishant Gopalakrishnan}

\seealso{ \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}}}
\examples{
\dontrun{
  data(qData)
  p1 <- qaProcess.marginevents(qData[[1]], channels=c("FSC-A", "SSC-A"),outdir=tempdir(),name=paste("Panel",1,sep="_"))
  writeQATextReport(qData[[1]],processes=list(p1),fileName="test.txt")
}

}

\keyword{ IO }
\keyword{ dynamic }
