## Classes and methods to implement agregated values of a quality
## assessment process on one individual flowFrame




## ===========================================================================
## Virtual agregator
## ---------------------------------------------------------------------------
## A class describing an agregated QA value for a single flowFrame. Derived
## subclasses describe the various subtypes of agregators. Slot 'frameID'
## stores the reference to the flowFrame and slot 'passed' contains a logical
## indicating whether the QA requirements have been met. Dedicated write
## methods of these subclasses produce the appropriate HTML output.
## ---------------------------------------------------------------------------
setClass("agregator",
         representation("VIRTUAL", frameID="character", passed="logical"),
         prototype=list(passed=TRUE))



## ===========================================================================
## binaryAgregator
## ---------------------------------------------------------------------------
## A class describing agregated QA value of the most simple binary type,
## i.e. indicating whether a QA requirement has been passed or not 
## ---------------------------------------------------------------------------
setClass("binaryAgregator",
         contains="agregator")


## write method to create HTML output
setMethod("write", signature("binaryAgregator", "file", "missing",
                             "missing", "missing"),
          function(x, file){
              if(x@passed)
                  writeLines(paste("<img class=\"agregator\"",
                                   "src=\"images/bulbGreen.png\">"), file)
              else
                  writeLines(paste("<img class=\"agregator\"",
                                   "src=\"images/bulbRed.png\">"), file) 
          })

## display details about agregator
setMethod("show", signature("binaryAgregator"),
          function(object){
              cat("Binary quality score ", ifelse(object@passed, "", "not"),
                  "passing the requirements\n(linked to flowFrame '",
                  object@frameID, "')\n", sep="") 
          })

## ===========================================================================
## factorAgregator
## ---------------------------------------------------------------------------
## A class describing agregated QA value of the factor type, i.e.
## indicating the states of the QA results from a selection of
## different outcomes
## ---------------------------------------------------------------------------
setClass("factorAgregator",
         representation(x="factor"),
         contains="agregator")


## write method to create HTML output
setMethod("write", signature("factorAgregator", "file", "missing",
                             "missing", "missing"),
          function(x, file){
              col <- ifelse(x@passed, "green", "red")
              lx <- levels(x@x)
              fcol <- rep("lightgray", length(lx))
              fcol[match(x@x, lx)] <- col
              writeLines("<b>", con)
              writeLines(paste("<span style=\"margin:0 3px; ",
                               "color:", fcol, ";\">", lx,
                               "</span>", sep=""), file)
              writeLines("</b><br>", con)
          })


## display details about agregator
setMethod("show", signature("factorAgregator"),
          function(object){
              cat("Factorized quality score ", ifelse(object@passed, "", "not"),
                  "passing the requirements of value=", object@x,
                  "\n(linked to flowFrame '",
                  object@frameID, "')\n", sep="") 
          })


## ===========================================================================
## stringAgregator
## ---------------------------------------------------------------------------
## A class describing agregated QA value of the string type, i.e. a
## character vector that was created by the QA process with a textual
## description of the result
## ---------------------------------------------------------------------------
setClass("stringAgregator",
         representation(x="character"),
         contains="agregator")

## write method to create HTML output
setMethod("write", signature("stringAgregator", "file", "missing",
                             "missing", "missing"),
          function(x, file){
              col <- ifelse(x@passed, "green", "red")
              writeLines(paste("<b><p style=\"margin:0 3px; ",
                               "color:", col, ";\">", x@x,
                               "</p></b>", sep=""), file)
          })

## display details about agregator
setMethod("show", signature("stringAgregator"),
          function(object){
              cat("Textual quality score ", ifelse(object@passed, "", "not"),
                  "passing the requirements of value=", object@x,
                  "\n(linked to flowFrame '",
                  object@frameID, "')\n", sep="") 
          })




## ===========================================================================
## numericAgregator
## ---------------------------------------------------------------------------
## A class describing agregated QA value of the numeric type, i.e. a
## numeric skalar that was created by the QA process
## ---------------------------------------------------------------------------
setClass("numericAgregator",
         representation(x="numeric"),
         contains="agregator")

## write method to create HTML output
setMethod("write", signature("numericAgregator", "file", "missing",
                             "missing", "missing"),
          function(x, file){
              col <- ifelse(x@passed, "green", "red")
              writeLines(paste("<b><span style=\"margin:0 3px; ",
                               "color:", col, ";\">", signif(x@x,2),
                               "</span></b><br>", sep=""), file)
          })

## display details about agregator
setMethod("show", signature("numericAgregator"),
          function(object){
              cat("Numeric quality score ", ifelse(object@passed, "", "not"),
                  "passing the requirements of value=", object@x,
                  "\n(linked to flowFrame '",
                  object@frameID, "')\n", sep="") 
          })


## ===========================================================================
## rangeAgregator
## ---------------------------------------------------------------------------
## A class describing agregated QA value of the range type, i.e. a
## numeric value within a defined range of values (e.g. a percentage)
## ---------------------------------------------------------------------------
setClass("rangeAgregator",
         representation(min="numeric", max="numeric"),
         contains="numericAgregator")
         
## write method to create HTML output
setMethod("write", signature("rangeAgregator", "file", "missing",
                             "missing", "missing"),
          function(x, file){
              perc <- (x@x-x@min)/diff(c(x@min, x@max))*100
              col <- ifelse(x@passed, "green", "red")
              writeLines(paste("<table class=\"agregator\" width=\"40px\" ",
                               "height=\"10px\" style=\"border:1px solid ",
                               "black; border-spacing:0px; padding:1px; ",
                               "margin:2px\">",
                               "<tr><td width=\"", perc, "%\" ",
                               "height=\"10px\" style=\"background-color:",
                               col, "; padding:0px; margin:0px;\"></td>",
                               "<td width=\"", 100-perc, "%\" ",
                               "style=\"padding:0px; margin:0px;\"></td>",
                               "</tr></table>", sep=""), file)
          })

## display details about agregator
setMethod("show", signature("rangeAgregator"),
          function(object){
              cat("Range quality score ", ifelse(object@passed, "", "not"),
                  "passing the requirements of value=", object@x,
                  "\n(linked to flowFrame '",
                  object@frameID, "')\n", sep="") 
          })


## ===========================================================================
## agregatorList
## ---------------------------------------------------------------------------
## A list of agregators for a whole flow set. All elements of the list must
## be agregators of equal type. The class should be initiated via it's
## constructor
## ---------------------------------------------------------------------------
setClass("agregatorList",
         contains="list")

setMethod("initialize", "agregatorList",
          function(.Object, ...) {
              if(is.list(..1))
              input <- ..1
              else
                  input <- list(...)
              if(!all(sapply(input, is, "agregator")))
                  stop("All items of an agregator list must ",
                       "inherit from class 'agregator'")
              if(length(unique(sapply(input, class))) !=1)
                  stop("All items of an agregator list must be ",
                       "of the same class")
              .Object@.Data=input
              return(.Object)
          })

agregatorList <- function(...)
    new("agregatorList", ...)

## display details about list
setMethod("show", signature("agregatorList"),
          function(object)
              cat("List of", length(object), "agregators\n")
          )


## ===========================================================================
## qaSummaryGraph
## ---------------------------------------------------------------------------
## A class that describes graphical output of a QA operation for a whole.
## flowSet, which is a summary over the outcome.
## This contains information about the type, the dimensions
## and the name of the image file. During initiation of the class we make
## sure that both bitmap and vector versions of the image are present and
## that the files are copied to the correct location. Object of the class
## should be created using the constructor
## ---------------------------------------------------------------------------
setClass("qaSummaryGraph",
         representation(fileName="character",
                        dimensions="numeric",
                        types="character"))

## copy and convert images during object instantiation
initiateGraphs <-  function(.Object, fileName, frameID, imageDir)
{
    ## get file information
    if(!file.exists(fileName))
        stop("Unable to find file '", fileName, "'")
    if(!file.exists(imageDir))
        dir.create(imageDir, recursive=TRUE)
    imageInfo <- strsplit(system(paste("identify", fileName),
                                 intern=TRUE), " ")[[1]][1:3]
    names(imageInfo) <- c("file", "type", "dimensions")
    bname <- basename(gsub("\\..*$", "", fileName))
    dims <- strsplit(imageInfo["dimensions"], "x")[[1]]
    .Object@dimensions <- c(width=as.numeric(dims[1]),
                            height=as.numeric(dims[2]))
    ## copy file to image directory
    cf <- file.path(imageDir, basename(fileName))
    if(!file.exists(cf))
        file.copy(fileName, cf)
    .Object@fileName <- cf
    ## convert image to pdf or jpg if missing
    convType <- ifelse(tolower(imageInfo["type"])=="pdf", "jpg", "pdf")
    newFileName <- file.path(imageDir, paste(bname, convType,
                                             sep="."))
    if(!file.exists(newFileName))
        system(paste("convert", fileName, newFileName))
    .Object@types <- c(tolower(imageInfo["type"]), convType)
    if(!missing(frameID))
        .Object@frameID <- as.character(frameID)
    return(.Object)
}

setMethod("initialize", "qaSummaryGraph",
          initiateGraphs
         )    

qaSummaryGraph <- function(fileName, imageDir="./images")
    new("qaSummaryGraph", fileName, imageDir=imageDir)


## display details about image files
setMethod("show", signature("qaSummaryGraph"),
          function(object)
              cat(object@types[1], " and ", object@types[2], " images '",
                 object@fileName, "' (",object@dimensions["width"],
                  "x", object@dimensions["height"],  ")\n", sep="")
          )

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
         representation(frameID="character"),
         contains="qaSummaryGraph")

setMethod("initialize", "qaGraph",
          initiateGraphs
         )      

qaGraph <- function(fileName, frameID, imageDir="./images")
    new("qaGraph", fileName, frameID, imageDir)

## display details about image files
setMethod("show", signature("qaGraph"),
          function(object)
              cat(object@types[1], " and ", object@types[2], " images '",
                 object@fileName, "' (",object@dimensions["width"],
                  "x", object@dimensions["height"],  ")\n",
                  "(linked to flowFrame '", object@frameID, "')\n",sep="")
          )

## ===========================================================================
## qaGraphList
## ---------------------------------------------------------------------------
## A list of qaGraphs for a whole flow set. Elements of the list are
## qaGraphs for the individual flowFrames. Objects of the class should
## be created using the constructor
## ---------------------------------------------------------------------------
setClass("qaGraphList",
         contains="list")

setMethod("initialize", "qaGraphList",
          function(.Object, imageFiles, frameIDs, imageDir) {
              if(!all(file.exists(imageFiles)))
                  stop("'imageFiles' must be character vector of ",
                       "paths to image files")
              if(missing(frameIDs))
                  frameIDs <- names(imageFiles)
              if(is.null(frameIDs))
                  stop("Need frameIDs to match files to frames")
              input <- mapply(qaGraph, imageFiles, as.character(frameIDs),
                              MoreArgs=list(imageDir=imageDir))
              .Object@.Data=input
              return(.Object)
          })

qaGraphList <- function(imageFiles, frameIDs, imageDir="./images")
    new("qaGraphList", imageFiles, frameIDs, imageDir)


## display details about list
setMethod("show", signature("qaGraphList"),
          function(object)
              cat("List of", length(object), "images\n")
          )



## ===========================================================================
## qaProcess
## ---------------------------------------------------------------------------
## A class that describes a QA process and the associated graphical output
## that is generated during this process. Objects of this class encapsulate
## all the necessary information to generate the HTML output including all
## links that is bundled in a QA report as generated by function
## 'createQAReport'
## ---------------------------------------------------------------------------
setClass("qaProcess",
         representation(name="character",
                        type="character",
                        frameGraphs="list",
                        summaryGraphs="list",
                        agregators="list",
                        frameIDs="character",
                        path="character"))


setMethod("initialize", "qaProcess",
          function(.Object, name, type, summaryGraphs,
                   frameGraphs, agregators, frameIDs, path) {
              if(!file.exists(path))
                  stop(path, " is not a valid directory")
              imageFiles <- c(sapply(frameGraphs, function(x)
                                     sapply(x, slot, "fileName")),
                              sapply(summaryGraphs, slot, "fileName"))
              missing <- !file.exists(imageFiles)
              if(any(missing))
                  stop(paste("Unable to find file", imagesFiles[missing])) 
              fids <- sapply(frameGraphs, function(x)
                                     sapply(x, slot, "frameID"))
              if(any(!fids %in% frameIDs))
                  stop("frameIDs don't match")
              .Object@name <- name
              .Object@type <- type
              .Object@summaryGraphs <- summaryGraphs
              .Object@frameGraphs <- frameGraphs
              .Object@agregators <- agregators
              .Object@frameIDs <- frameIDs
              .Object@path <- path
              return(.Object)
          })




## display details about process
setMethod("show", signature("qaProcess"),
          function(object)
              cat("Quality process '", object@name, "' of type '",
                  object@type, "'\n", sep="")
          )





library(geneplotter)


con <- openHtmlPage("test")

y <- new("binaryAgregator")
write(y, con)
y@passed <- FALSE
write(y, con)


x <- new("rangeAgregator", x=11, min=3, max=13)
write(x, con)
x@passed <- FALSE
x@x <- 6
write(x, con)


f <- factor("b", levels=c("a", "b", "c"))
z <- new("factorAgregator", x=f)
write(z, con)
z@passed <- FALSE
write(z, con)


n <- new("numericAgregator", x=0.54212123)
write(n, con)
n@passed <- FALSE
write(n, con)


s <- new("stringAgregator", x="This is sample output")
write(s, con)
s@passed <- FALSE
write(s, con)




closeHtmlPage(con)
