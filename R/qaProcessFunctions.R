## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address.
guid <- function()
    format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"])


## QA process indicating too many events on the margins
## FIXME: include grouping
qaProcess.marginevents <- function(set, channels=NULL, outdir, cFactor=3)
{
    ## count events on the margins
    frameIDs <- sampleNames(set)
    if(is.null(channels))
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
    else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    ranges <- range(set[[1]], parms)
    sums <- matrix(ncol=length(set), nrow=length(parms),
               dimnames=list(parms, frameIDs))
    cat("computing margin events...")
    for(p in parms){
        exp <- parse(text=paste("`", p, "`==", ranges[,p], sep="",
                     collapse="|"))
        attributes(exp) <- NULL
        ef <- expressionFilter(exp, filterId=p)
        sums[p,] <- sapply(filter(set, ef), function(x) summary(x)$p)
        cat(".")
    }
    cat("\n")

    ## create summary plot
    require("lattice")
    gid <- guid()
    tmp <- tempdir()
    sfile <- file.path(tmp, "summary.pdf")
    pdf(file=sfile)
    col.regions=colorRampPalette(c("white",  "darkblue"))(256)
    print(levelplot(t(sums)*100, scales = list(x = list(rot = 90)),
                    xlab="", ylab="", main="% margin events",
                    col.regions=col.regions))
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350)
    
    ## create graphs and aggregators for each frame (and each channel)
    cat("creating frame plots...")
    frameProcesses <- list()
    for(i in 1:length(set)){
        fnames <- NULL
        agTmp <- aggregatorList()
        for(j in 1:length(parms)){
            tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), "_",
                                          parms[j], ".pdf", sep=""))
            pdf(file=tfile, height=3)
            plot(density(exprs(set[[i]][,j])), ann=FALSE, axes=FALSE,
                 col="darkred", lwd=2)
            axis(1)
            dev.off()
            fnames <- c(fnames, tfile)
            agTmp[[j]] <- new("rangeAggregator",
                              passed=sums[j,i]<(mean(sums[j,])*cFactor),
                              x=sums[j,i], min=0, max=1)
            cat(".")
        }
        names(agTmp) <- parms
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir, width=150)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                summaryAggregator=ba,
                                                frameAggregators=agTmp,
                                                frameGraphs=fGraphs)
    }
    
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="margin events",
               type="margin events", frameIDs=frameIDs, summaryGraph=sgraph,
               frameProcesses=frameProcesses))
    
}    

        

## QA process indicating strange patterns over time
qaProcess.timeline <- function(set, channel, outdir, cutoff=0.1)
{
    ## create summary plot and its associated qaGraph object
    if(length(channel)!=1)
        stop("'channel' must be of length 1")
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    sfile <- file.path(tmp, "summary.pdf")
    pdf(file=sfile)
    summary <- timeLinePlot(set, channel)
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350)

    ## create graphs and aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in 1:length(set)){
        tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".pdf",
                                      sep=""))
        pdf(file=tfile)
        timeLinePlot(set[i], channel)
        dev.off()
        ba <- new("binaryAggregator", passed=summary[i]<cutoff)
        fg <- qaGraph(fileName=tfile, imageDir=idir, width=220)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <- qaProcessFrame(fid, ba, fg)
        cat(".")
    }

    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=paste("timeLine", channel),
               type="time line", frameIDs=frameIDs, summaryGraph=sgraph,
               frameProcesses=frameProcesses))
    
}    


allDens <- function(sets, channel){
    extDens <- function(x, channel)
        fsApply(x, function(x)
            {
                tmp <- density(x[,channel])
                return(rbind(tmp$x, tmp$y))
            }, use.exprs=TRUE)
    if(is.list(sets))
        dens <- sapply(sets, extDens, channel)
    else
        dens <- extDens(sets, channel)
    return(t(dens))
}


## QA process to compare KL distances
qaProcess.Similarity <- function(set, channel, outdir, groups=NULL, cutoff=0.1)
{
    ## create summary plot and its associated qaGraph object
    if(length(channel)!=1)
        stop("'channel' must be of length 1")
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    sfile <- file.path(tmp, "summary.pdf")
    dens <- t(allDens(set, channel))
    col <- if(is.null(grps)) "black" else as.integer(factor(groups))
    pdf(file=sfile)
    matplot(dens[seq(1, nrow(dens), by=2), ], dens[seq(2, nrow(dens), by=2), ],
            type="l", col=col)
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350)

    ## create graphs and aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in 1:length(set)){
        tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".pdf",
                                      sep=""))
        pdf(file=tfile)
        timeLinePlot(set[i], channel)
        dev.off()
        ba <- new("binaryAggregator", passed=summary[i]<cutoff)
        fg <- qaGraph(fileName=tfile, imageDir=idir, width=220)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <- qaProcessFrame(fid, ba, fg)
        cat(".")
    }

    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=paste("timeLine", channel),
               type="time line", frameIDs=frameIDs, summaryGraph=sgraph,
               frameProcesses=frameProcesses))
    
}    


