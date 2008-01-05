qaProcess.timeline <- function(set, channel, outdir)
{
    id <- paste("timeLine", channel, sep="_")
    tmp <- tempdir()
    sfile <- file.path(tmp, paste(id, "summary.jpg", sep="_"))
    jpeg(file=sfile)
    summary <- timeLinePlot(set, channel)
    dev.off()
    idir <- file.path(outdir, "images")
    sgraph <- qaSummaryGraph(sfile, idir)
    aggr <- list()
    fileNames <- NULL
    for(i in 1:length(set)){
        fid <- sampleNames(set)[i]
        tfile <- file.path(tmp, paste(id, "_", sprintf("%0.2d", i), ".jpg",
                                      sep=""))
        jpeg(file=tfile)
        timeLinePlot(set[i], channel)
        dev.off()
        aggr[[i]] <- new("binaryAgregator", frameID=fid, passed=TRUE)
        fileNames <- c(fileNames, tfile)
    }
    glist <- qaGraphList(fileNames, sampleNames(set), idir)
    alist <- agregatorList(aggr)
    return(new("qaProcess", name=paste("timeLine", channel), type="time line",
               frameGraphs=list(glist), summaryGraphs=list(sgraph),
               agregators=list(alist), frameIDs=sampleNames(set),
               path=outdir))
    
}    
