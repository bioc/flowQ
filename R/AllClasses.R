## ===========================================================================
## virtual outlier test
## ---------------------------------------------------------------------------
setClass("outlier", 
         representation("VIRTUAL",test="character",
                        parameters="ANY"))


## ===========================================================================
## outlierResult
## ---------------------------------------------------------------------------
## A container for the results of outliers test 
## ---------------------------------------------------------------------------
setClass("outlierResult",
         representation(frameId="character", filterDetails="list"),
         contains="outlier",
         prototype=list(frameId=character(0), filterDetails=list()))




## ===========================================================================
## Virtual qaAggregator
## ---------------------------------------------------------------------------
## A class describing an aggregated QA value for a single flowFrame. Derived
## subclasses describe the various subtypes of aggregators. Slot 'frameID'
## stores the reference to the flowFrame (i.e., the sample name in the
## respective flowSet) and slot 'passed' contains a logical
## indicating whether the QA requirements have been met. Dedicated write
## methods of these subclasses produce the appropriate HTML output.
## ---------------------------------------------------------------------------
setClass("qaAggregator",
         representation("VIRTUAL", passed="logical"),
         prototype=list(passed=TRUE))



## ===========================================================================
## binaryAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the most simple binary type,
## i.e. indicating whether a QA requirement has been passed or not 
## ---------------------------------------------------------------------------
setClass("binaryAggregator",
         contains="qaAggregator")



## ===========================================================================
## discreteAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the most discrete type,
## i.e. it can have three states: passed(1), warning(2) and failed(0) 
## ---------------------------------------------------------------------------
setClass("discreteAggregator",
         representation(x="factor"),
         contains="qaAggregator",
         prototype=list(x=factor(1)))



## ===========================================================================
## factorAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the factor type, i.e.
## indicating the states of the QA results from a selection of
## different outcomes
## ---------------------------------------------------------------------------
setClass("factorAggregator",
         representation(x="factor"),
         contains="qaAggregator")



## ===========================================================================
## stringAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the string type, i.e. a
## character vector that was created by the QA process with a textual
## description of the result
## ---------------------------------------------------------------------------
setClass("stringAggregator",
         representation(x="character"),
         contains="qaAggregator")



## ===========================================================================
## numericAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the numeric type, i.e. a
## numeric skalar that was created by the QA process
## ---------------------------------------------------------------------------
setClass("numericAggregator",
         representation(x="numeric"),
         contains="qaAggregator")



## ===========================================================================
## rangeAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the range type, i.e. a
## numeric value within a defined range of values (e.g. a percentage)
## ---------------------------------------------------------------------------
setClass("rangeAggregator",
         representation(min="numeric", max="numeric"),
         contains="numericAggregator")
         


## ===========================================================================
## aggregatorList
## ---------------------------------------------------------------------------
## A list of aggregators for a whole flow set. All elements of the list must
## be aggregators. The class should be initiated via it's
## constructor
## ---------------------------------------------------------------------------
setClass("aggregatorList",
         contains="list")

setMethod("initialize", "aggregatorList",
          function(.Object, ...) {
              if(length(list(...))>0){
                  if(is.list(..1))
                      input <- ..1
                  else
                      input <- list(...)
                  if(!all(sapply(input, is, "qaAggregator")))
                      stop("All items of an aggregator list must ",
                           "inherit from class 'qaAggregator'")
                  .Object@.Data=input
              }
              return(.Object)
          })

aggregatorList <- function(...)
    new("aggregatorList", ...)



## ===========================================================================
## qaGraph
## ---------------------------------------------------------------------------
## A class that describes graphical output of a QA operation for a single
## flowFrame. This contains information about the type, the dimensions
## and the name of the image file. During initiation of the class we make
## sure that both bitmap and vector versions of the image are present and
## that the files are copied to the correct location. Object of the class
## should be created using the constructor
## ---------------------------------------------------------------------------
setClass("qaGraph",
         representation(fileNames="character",
                        dimensions="matrix",
                        types="character",
                        id="character"))
         
setMethod("initialize", "qaGraph",
          function(.Object, fileName, imageDir, width=NULL, empty=FALSE){
              if(!empty){
                  ## check arguments
                  if(!is.null(width) && (length(width)!=1 ||
                                         !is.numeric(width)))
                      stop("'width' must be numeric scalar")
                  if(!file.exists(fileName))
                      stop("Unable to find file '", fileName, "'")
                  
                  ## get file information
                  if(!file.exists(imageDir))
                      dir.create(imageDir, recursive=TRUE)
                  imageInfo <- strsplit(system(paste("identify", fileName),
                                               intern=TRUE), " ")[[1]][1:3]
                  names(imageInfo) <- c("file", "type", "dimensions")
                  bname <- basename(gsub("\\..*$", "", fileName))
                  dims <- as.numeric(strsplit(imageInfo["dimensions"],"x")[[1]])
                  newDims <- dims
                  if(!is.null(width)){
                      scaleFac <- width/dims[1]
                      if(scaleFac!=1)
                          newDims <-  dims*scaleFac
                  }
                  ft <- c("vectorized", "bitmap")
                  cf <- file.path(imageDir, basename(fileName))
                  
                  ## convert image to vectorized or bitmap version
                  if(tolower(imageInfo["type"])=="pdf"){
                      ## original image is vectorized
                      convType <- "jpg"
                      newFileName <- file.path(imageDir, paste(bname, convType,
                                                               sep="."))
                      if(!file.exists(newFileName))
                          system(paste("convert -resize", paste(newDims,
                                                                collapse="x"),
                                       fileName, newFileName))
                      type <- c("pdf", "jpg")
                      files <- c(cf, newFileName)
                      if(!file.exists(cf))
                          file.copy(fileName, cf)
                  }else{
                      ## original image is bitmap
                      convType <- "pdf"
                      newFileName <- file.path(imageDir, paste(bname, convType,
                                                               sep="."))
                      if(!file.exists(newFileName))
                          system(paste("convert", fileName, newFileName))
                      if(!file.exists(cf))
                          system(paste("convert -resize", paste(newDims,
                                                                collapse="x"),
                                       fileName, cf))
                      type <- c("pdf",tolower(imageInfo["type"]))
                      files <- c(newFileName, cf)
                  }
                  ## fill qaGraph object
                  .Object@dimensions <-  matrix(c(dims, newDims), ncol=2,
                                                byrow=TRUE,
                                                dimnames=list(ft, c("width",
                                                "height")))
                  names(type) <- names(files) <- ft
                  .Object@types <- type
                  .Object@fileNames <- files
                  .Object@id=guid()
              }
              return(.Object)
          })


qaGraph <- function(...)
    new("qaGraph", ...)



## ===========================================================================
## qaGraphList
## ---------------------------------------------------------------------------
## A list of qaGraphs. Elements of the list are individual qaGraphs. This
## mainly exists for the sake of a constructor to facilitate object creation
## ---------------------------------------------------------------------------
setClass("qaGraphList",
         contains="list")

setMethod("initialize", "qaGraphList",
          function(.Object, imageFiles, imageDir, width=NULL) {
              if(!all(file.exists(imageFiles)))
                  stop("'imageFiles' must be character vector of ",
                       "paths to image files")
              
              input <- lapply(imageFiles, qaGraph, imageDir=imageDir,
                              width=width)
              .Object@.Data=input
              return(.Object)
          })

qaGraphList <- function(...)
    new("qaGraphList", ...)



## ===========================================================================
## qaProcessFrame
## ---------------------------------------------------------------------------
## A class that bundles all information about the QA output for a single
## flowFrame. Slots of this class are:
##   id: a unique identifier for the whole block
##   frameID: the id of the respective flowFrame
##   summaryAggregator: a binaryAggregator indicating status
##   summaryGraph: an optional image linked to the summaryAggregator
##   frameAggregators: a list of aggregators for this frame
##   frameGraphs: a list of optional images linked to the respective
##                frameaggregators
## ---------------------------------------------------------------------------
setClass("qaProcessFrame",
         representation(id="character",
                        frameID="character",
                        summaryAggregator="qaAggregator",
                        summaryGraph="qaGraph",
                        frameAggregators="aggregatorList",
                        frameGraphs="qaGraphList"))



setMethod("initialize", "qaProcessFrame",
          function(.Object, frameID, summaryAggregator, summaryGraph,
                   frameAggregators, frameGraphs){
              .Object@id <- guid()
              .Object@frameID <- frameID
              .Object@summaryAggregator <- summaryAggregator
              if(missing(summaryGraph))
                  summaryGraph <- new("qaGraph", empty=TRUE)
              .Object@summaryGraph <- summaryGraph
              if(xor(missing(frameAggregators), missing(frameGraphs)))
                 stop("Both 'frameAggregators' and 'frameGraphs' must be ",
                       "specified")
              else if(!missing(frameAggregators)){
                 .Object@frameAggregators <- frameAggregators
                 .Object@frameGraphs <- frameGraphs
             }
              return(.Object)
          })

qaProcessFrame <- function(...)
    new("qaProcessFrame", ...)



## ===========================================================================
## qaProcess
## ---------------------------------------------------------------------------
## A class that describes a QA process and the associated graphical output
## that is generated during this process. Objects of this class encapsulate
## all the necessary information to generate the HTML output including all
## links that is bundled in a QA report as generated by function
## 'qaReport'
## ---------------------------------------------------------------------------
setClass("qaProcess",
         representation(id="character",
                        name="character",
                        type="character",
                        frameIDs="character",
                        summaryGraph="qaGraph",
                        frameProcesses="list"))
                        


setMethod("initialize", "qaProcess",
          function(.Object, id, name, type, frameIDs, summaryGraph,
                   frameProcesses) {

              
              fids <- names(frameProcesses)
              if(is.null(fids) || (!fids %in% frameIDs))
                  stop("frameIDs don't match")
              if(missing(name))
                  name <- "anonymous QA Process"
              .Object@id <- id
              .Object@name <- name
              .Object@type <- type
              .Object@frameIDs <- frameIDs
              .Object@summaryGraph <- summaryGraph
              .Object@frameProcesses <- frameProcesses
              validProcess(.Object)
              return(.Object)
          })

qaProcess <- function(...)
    new("qaProcess", ...)



