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
