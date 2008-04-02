## re-evaluate a process for new thresholds or cutoff values 
evaluateProcess <- function(process, thresh, ...)
{
    switch(process@type,
           "time line"=
       {
           efun <- function(x, c){
               qaScore <- x@details$qaScore
               x@summaryAggregator@passed <- qaScore <= c
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh)    
       },
           "margin events"=
       {
           efun <- function(x, c){
               sums <- x@details$events
               msums <- x@details$mevents
               for(i in 1:length(x@frameAggregators)){
                   x@frameAggregators[[i]]@passed <- sums[i] < (msums[i]*c)
                   x@frameAggregators[[i]]@x <- sums[i]
               }
               nfail <- !sapply(x@frameAggregators, slot, "passed")
               val <- if(sum(nfail)==1) factor(2) else factor(0)
               if(sum(nfail)==0)
                   val <- factor(1)
               x@summaryAggregator <- discreteAggregator(val)
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh)    
       },           
           stop("Don't know how to deal with process of type '",  process@type,
                "'", call.=FALSE)
           )
    return(process)
}
