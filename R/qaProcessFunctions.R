## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address.
guid <- function()
    format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"])




## QA process indicating too many events on the margins in comparison to
## average number of margin events for a particular channel. The 'grouping'
## argument can be used to compare within groups. 'cFactor' is the cutoff
## value, essentially this is a factor in the confidence intervall defined
## by the standard deviation of the average number of margin events for a
## particular channel.
qaProcess.marginevents <- function(set, channels=NULL, grouping=NULL, outdir,
                                   cFactor=1, pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
       stop("'set' needs to be of class 'flowSet'")
    frameIDs <- sampleNames(set)
    ls <- length(set)
    if(is.null(channels)){
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
    }else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    if(!is.null(grouping))
        if(!is.character(grouping) || ! grouping %in% colnames(pData(set)))
            stop("'grouping' must be a character indicating one of the ",
                 "phenotypic variables in 'set'")
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cFactor) || length(cFactor)!=1)
        stop("'cFactor' must be numeric scalar")
    lp <- length(parms)
    
    ## count events on the margins using an expression filter
    ranges <- range(set[[1]], parms)
    perc <- matrix(ncol=ls, nrow=lp, dimnames=list(parms, frameIDs))
    cat("computing margin events...")
    for(p in parms){
        ef <- char2ExpressionFilter(paste("`", p, "`==", ranges[,p], sep="",
                     collapse="|"), filterId=p)
        ff <- filter(set, ef)
        perc[p,] <- sapply(ff, function(x) summary(x)$p)
        cat(".")
    }
    cat("\n")
    
    ## create summary plot
    require("lattice")
    gid <- guid()
    tmp <- tempdir()
    ## Necessary because the directory path in Windows returned by tmpdir is odd
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfile <- file.path(tmp, "summary.pdf")
    pdf(file=sfile)
    col.regions=colorRampPalette(c("white",  "darkblue"))(256)
    print(levelplot(t(perc)*100, scales = list(x = list(rot = 90)),
                    xlab="", ylab="", main="% margin events",
                    col.regions=col.regions))
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## deal with groups if there are any
    frameProcesses <- list()
    grps <- if(!is.null(grouping)) pData(set)[,grouping] else rep(1,ls)
    grps <- split(1:ls, grps, drop=TRUE)
    allData <- lapply(grps, function(x, set)
                      exprs(as(set[x,parms], "flowFrame")), set)
    
    ## create graphs and aggregators for each frame (and each channel)
    cat("creating frame plots...")
    for(i in 1:length(set)){
        fnames <- NULL
        agTmp <- aggregatorList()
        thisGrp <-
            if(is.null(grouping)) "1" else as.character(pData(set)[i,grouping])
        mv <- sv <- NULL
        for(j in 1:length(parms)){
            tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), "_",
                                          gsub("\\..*$", "", parms[j]), ".pdf",
                                          sep=""))
            pdf(file=tfile, height=3)
            par(mar=c(1,0,1,0))
            ll <- list(allData[[thisGrp]][,j], exprs(set[[i]][,parms[j]]))
            bw <- bw.nrd0(ll[[2]])
            dens <- lapply(ll, density, bw=bw)
            rl <- range(unlist(ll))
            xlim <- rl + c(-1,1) * (diff(rl)/7)
            ylim <- c(0, max(c(dens[[1]]$y, dens[[2]]$y)))
            plot(1,1,xlim=xlim, ylim=ylim, type="n", axes=FALSE, ann=FALSE)
            polygon(dens[[1]], col="gray", border="gray")
            lines(dens[[2]], col="darkred", lwd=3)
            dev.off()
            fnames <- c(fnames, tfile)
            
            m <- mean(perc[j,grps[[thisGrp]]])
            s <- sd(perc[j,grps[[thisGrp]]])
            if(is.na(s))
                s <- sd(perc[j,])
            passed <- perc[j,i] <= m+s*cFactor & perc[j,i] >= m-s*cFactor
            mv <- c(mv, m)
            sv <- c(sv,s)
            agTmp[[j]] <- new("rangeAggregator", passed=passed,
                              x=perc[j,i], min=0, max=1)
            cat(".")
        }

        ## summarize the results for separate channels and bundle things up
        names(agTmp) <- parms
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir, width=150, pdf=pdf)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                summaryAggregator=ba,
                                                frameAggregators=agTmp,
                                                frameGraphs=fGraphs,
                                                details=list(events=perc[,i],
                                                mevents=rowMeans(perc),
                                                m=mv, s=sv))
    }
    
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="margin events",
               type="margin events", summaryGraph=sgraph,
               frameProcesses=frameProcesses))
}    


## QA process indicating strange patterns in signal intensity over time.
## This is done for all channels at once now, with drilldown to the separate
## channel results. The cutoff is the variance cutoff used directly in function
## 'timeLinePlot'
qaProcess.timeline <- function(set, channels=NULL, outdir, cutoff=1,
                               name="time line",
                               sum.dimensions=NULL, det.dimensions=c(7,7),
                               pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    ls <- length(set)
    if(is.null(channels))
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
    else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    lp <- length(parms)
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cutoff) || length(cutoff)!=1)
        stop("'cutoff' must be numeric scalar")
    if((!is.null(sum.dimensions) && !all(is.numeric(sum.dimensions))) ||
       !all(is.numeric(det.dimensions)))
        stop("Plot dimensions must be provided as numerics")
    if(is.null(sum.dimensions))
        sum.dimensions <- c(4*lp, 7)
    det.dimensions <- rep(det.dimensions, 2)
    
    ## create summary plots for each channel
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    ## Necessary because the directory path in Windows returned by tmpdir is odd
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
    summary <- vector(lp, mode="list")
    for(j in seq_len(lp)){
        sfile <- file.path(tmp, paste("summary_", j, ".pdf", sep=""))
        pdf(file=sfile, width=sum.dimensions[1]/lp,
            height=sum.dimensions[2])
        binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
        summary[[j]] <- timeLinePlot(set, parms[j], binSize=binSize,
                                     varCut=cutoff)
        dev.off()
        sfiles <- c(sfiles, sfile)
        cat(".")
    }
    
    ##glue together the summary graphs and create a qaGraph object
    sfile <- paste(tmp, "summary.pdf", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=max(350,150*lp), pdf=pdf)

    ## create graphs and aggregators for each frame and channel
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in seq_len(ls)){
        fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp)){
            tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), "_",
                                          gsub("\\..*$", "", parms[j]), ".pdf",
                                          sep=""))
            pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
            timeLinePlot(set[i], parms[j],
                         main=paste("score=", signif(summary[[j]][i], 4), sep=""),
                         cex.main=2, binSize=binSize, varCut=cutoff)
            dev.off()
            fnames <- c(fnames, tfile)
            agTmp[[j]] <- new("numericAggregator", passed=summary[[j]][i]<=0,
                              x=summary[[j]][i])
            cat(".")
        }

        ## wrap graphs and aggregators in objects
        names(agTmp) <- parms
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
                               width=min(220, lp*120), pdf=pdf)
        fid <- frameIDs[i]
        dr <- lapply(seq_len(lp), function(x) attr(summary[[x]], "raw")[[i]])
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                summaryAggregator=ba,
                                                frameAggregators=agTmp,
                                                frameGraphs=fGraphs,
                                                details=list(raw=dr))
    }
            
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="time line",
                     type="time line", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}    

    


## Detect distrubances in the flow of cells over time
qaProcess.timeflow <- function(set, outdir, cutoff=2, name="time flow",
                               sum.dimensions=c(7,7), det.dimensions=c(7,7),
                               pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    ls <- length(set)
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cutoff) || length(cutoff)!=1)
        stop("'cutoff' must be numeric scalar")
    if(!all(is.numeric(sum.dimensions)) ||
       !all(is.numeric(det.dimensions)))
        stop("Plot dimensions must be provided as numerics")
    det.dimensions <- rep(det.dimensions, 2)
    sum.dimensions <- rep(sum.dimensions, 2)

    
    ## create summary plot and its associated qaGraph object
    sum.dimensions <- rep(sum.dimensions, 2)
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    ## Necessary because the directory path in Windows returned by tmpdir is odd
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
    sfile <- file.path(tmp, "summary.pdf")
    pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
    summary <- timeLinePlot(set, colnames(set)[[1]], binSize=binSize,
                            varCut=cutoff, type="frequency")
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## create graphs and aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    det.dimensions <- rep(det.dimensions, 2)
    for(i in seq_len(ls)){
        tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".pdf",
                                      sep=""))
        pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
        sum <- timeLinePlot(set[i], colnames(set)[[1]],
                     main=paste("score=", signif(summary[i], 4), sep=""),
                     cex.main=2, binSize=binSize, type="frequency",
                     varCut=cutoff)
        dev.off()
        ba <- new("binaryAggregator", passed=summary[i]<=0)
        fg <- qaGraph(fileName=tfile, imageDir=idir, width=220, pdf=pdf)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <-
            qaProcessFrame(fid, ba, fg, details=list(qaScore=sum))
        cat(".")
    }

    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name, type="time flow",
                     summaryGraph=sgraph, frameProcesses=frameProcesses))
}    


## Detect unusually low cell counts
qaProcess.cellnumber <- function(set, grouping=NULL, outdir, cFactor=0.5,
                                 name="cell number", sum.dimensions=c(7,7),
                                 pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    ls <- length(set)
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cFactor) || length(cFactor)!=1)
        stop("'cutoff' must be numeric scalar")
    if(!is.null(grouping))
        if(!is.character(grouping) || ! grouping %in% colnames(pData(set)))
            stop("'grouping' must be a character indicating one of the ",
                 "phenotypic variables in 'set'")
    if(!all(is.numeric(sum.dimensions)))
        stop("Plot dimensions must be provided as numerics")
    sum.dimensions <- rep(sum.dimensions, 2)

    ## deal with groups if there are any
    grps <- if(!is.null(grouping)) pData(set)[,grouping] else rep(1,ls)
    grpsi <- split(1:ls, grps, drop=TRUE)
    sset <- split(set, grps)
    
    ## create summary plot and its associated qaGraph object
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    ## Necessary because the directory path in Windows returned by tmpdir is odd
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    cellNumbers <- as.numeric(fsApply(set, nrow))
    sfile <- file.path(tmp, "summary.pdf")
    pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
    col <- "gray"
    par(mar=c(10.1, 4.1, 4.1, 2.1), las=2)
    barplot(cellNumbers, col=col, border=NA, names.arg=sampleNames(set),
            cex.names=0.8, cex.axis=0.8)
    abline(h=mean(cellNumbers), lty=3, lwd=2)
    dev.off()
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## create aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in seq_len(ls)){
        thisGrp <-
            if(is.null(grouping)) "1" else as.character(pData(set)[i,grouping])
        summary <- -(log(cellNumbers[i] / mean(cellNumbers[grpsi[[thisGrp]]])))
        ba <- new("numericAggregator", x=cellNumbers[i],
                  passed=summary<cFactor)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <-
            qaProcessFrame(fid, ba, details=list(qaScore=summary))
        cat(".")
    }

    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name, type="cell number",
                     summaryGraph=sgraph, frameProcesses=frameProcesses))
}    


multiSampleStats <- function(flowList,channels,func)
{ 
	alqLen<-length(flowList)
	if(is.null(channels))
		stop("Please specify the parameters to be compared")
	if(!all(channels %in% colnames(flowList[[1]])))
		stop("Invalid channel(s)")
	res<-data.frame()
	for( i in seq_len(alqLen))
	{ 
	    alqNo<-rep(i,length(flowList[[i]]))
	    temp1<-fsApply(flowList[[i]][,channels[1]],each_col,func)
	    temp2<-fsApply(flowList[[i]][,channels[2]],each_col,func)
	    newres<-data.frame(alqNo,temp1, temp2, patientID=rownames(temp1),check.names=FALSE)
	    res<-rbind(res,newres)
	}
	rownames(res)<-NULL
	patientID=sampleNames(flowList[[1]])
	formula<-paste("`",channels[1],"`"," ", "~"," ","`",channels[2],"`","|","patientID",sep="")
	
	graph<-xyplot(eval(parse(text=formula)),data=res, 
		    auto.key=list(space="right"),
		    par.settings = simpleTheme(col=1:alqLen),
		    groups=alqNo, type="p")
	
	return(list(graph=graph,stats=res))
}


qaProcess.summaryOutLiers <- function(flowList, channels=c("FSC-A","SSC-A"),
                             outdir="QAReport",cutoff=4,func=median,
                             det.dimensions=c(300,300),pdf=FALSE,...)
	{
	if(!is(flowList, "list"))
		stop("'flowList' needs to be of class 'List'")
	
	if(!is(flowList[[1]], "flowSet"))
		stop("'flowList' needs to be a List of flowSets")
	
	if(!all(channels %in% colnames(flowList[[1]])))
		stop("Invalid channel(s)")
	
	if(!is.character(outdir) || length(outdir)!=1)
		stop("'outdir' must be a valid file path")
	
	if(!is.numeric(cutoff) || length(cutoff)!=1)
		stop("'cutoff' must be numeric scalar")
	
	cat("creating summary plots...")
	summary<-multiSampleStats(flowList,channels=channels,func)
	cn <- levels(summary$stat[,4])
	colScheme<-list()
	thecol=""
	panelFlag=rep(FALSE,length(cn))
	outLierFlag<-rep(FALSE,nrow(summary$stat))
	nonOutLierCount=rep(0,length(cn))   
	summary$stat<-cbind(summary$stat,outLierFlag)
	
	for( i in seq_len(length(cn)))
	{
		subset <- summary$stats[summary$stats$patientID==cn[i], 2:3] 
		res<-sign2(subset)
		panelFlag[i]=all(res$x.dist<=cutoff)
		nonOutLierCount[i]<-sum(res$x.dist>cutoff)
		indx<-which(summary$stat$patientID==cn[i])
		val<-res$x.dist>cutoff
		summary$stat$outLierFlag[indx]<-val
		cat(".")
	}
	gid <- guid()
	tmp <- tempdir()
	tmp <- gsub("\\", "/", tmp, fixed=TRUE)
	
	formula<-paste("`",channels[1],"`"," ", "~"," ","`",channels[2],"`","|",colnames(summary$stat)[4],sep="")
	grph<-xyplot( eval(parse(text=formula)),data<-summary$stat,groups=summary$stat[,5],col=c("black","red"))
	
	sfile <- file.path(tmp, "summary.jpg")    
	jpeg(file=sfile)   
	print(grph)
	dev.off()

	idir <- file.path(outdir, "images", gid)
	sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350,pdf=pdf)
	
        frameProcesses <- list()
	cat("\nCreating frame plots...")
	ls <- length(cn)
	for(i in seq_len(ls))
	{  
		fid <- cn[i]
		tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".jpg", sep=""))
		jpeg(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
		print(update(grph, index.cond = list(i))) 
		dev.off()
		fGraph <- qaGraph(fileName=tfile, imageDir=idir, pdf=pdf)
		agTmp <- new("numericAggregator", passed=panelFlag[i],x=nonOutLierCount[i])
		frameProcesses[[fid]] <- qaProcessFrame(fid, agTmp, fGraph)
		cat(".")
	} 
	output<-qaProcess(id=gid, 
			name=paste(channels[1]," / ",channels[2]," ",substitute(func),sep=""),
			type=paste(channels[1]," / ",channels[2], " Comparision" ,sep=""),
			summaryGraph=sgraph, 
			frameProcesses=frameProcesses)
	return(output)
}	


qaProcess.densityOutLiers <- function(flowList,channel=c("FSC-A"),
		outdir="QAReport",cutoff=0.11,det.dimensions=c(300,300),
                func=euc,pdf=FALSE,...)
{   
	if(!is(flowList, "list"))
		stop("'flowList' needs to be of class 'List'")
	
	if(!is(flowList[[1]], "flowSet"))
		stop("'flowList' needs to be a List of flowSets")
	
	if(!all(channel %in% colnames(flowList[[1]])))
		stop("Invalid channel")
	
	if(length(channel)!=1)
	{
		stop("Only one channel is to be specified")
	}
	
	if(!is.character(outdir) || length(outdir)!=1)
		stop("'outdir' must be a valid file path")
	
	if(!is.numeric(cutoff) || length(cutoff)!=1)
		stop("'cutoff' must be numeric scalar")
	
	cat("creating summary plots...")
	gid <- guid()
	tmp <- tempdir()
	tmp <- gsub("\\", "/", tmp, fixed=TRUE)
	alqLen<-length(flowList)
	patientID=sampleNames(flowList[[1]])
	panelFlag=FALSE
	tempgrph<-list()
	tempDist<-list()
	formula<-paste("~","`",channel[1],"`","|","Patient",sep="")
	
	ymax<-0
	xmax<-0
	xmin<-100000
	
	for( i in patientID)
	{   res<-data.frame()
	    tempStat<-matrix(ncol=256,nrow=alqLen)
	    for(j in seq_len(alqLen))			
	    {   
		value=exprs(flowList[[j]]@frames[[i]][,channel])
		newres<-data.frame(Patient=rep(i,nrow(value)),
                                   Aliquot=rep(j,nrow(value)),
				   data=value,check.names=FALSE)
		res<-rbind(res,newres)
                valRange<-range(value)
                if(xmax<valRange[2])
		    xmax<-valRange[2]
		if(xmin> valRange[1])
		    xmin<-valRange[1]
                
                tempStat[j,]<- density(value,n=256,from=valRange[1],to=valRange[2])$y
		yrng<-range(tempStat[j,])[2]
		if(yrng>ymax)
		    ymax<-yrng
	    }
		tempDist[[i]]<-sum(euc(tempStat))
		tempgrph[[i]]<-densityplot(eval(parse(text=formula)),data=res,groups=Aliquot,plot.points=FALSE)
		cat(".")
	}
	
	xrange<-c(0.9*xmin,1.1*xmax)
	
        grph<-densityplot(~x|patientID, 
			  data = list(patientID = factor(names(tempgrph),
			  levels = names(tempgrph)),x = seq_along(tempgrph)), 
			  xlim =xrange,ylim=c(0,1.1*ymax),
			  xlab=as.character(channel[1]),
			  key=simpleKey(text=as.character(1:alqLen),space="right"),
			  panel = function(x, y, plot.list) {
			      do.call(panel.densityplot, trellis.panelArgs(tempgrph[[x]], 1))
			      }
			  )	
	sfile <- file.path(tmp, "summary.jpg")    
	jpeg(file=sfile)   
	print(grph)
	dev.off()
	
	idir <- file.path(outdir, "images", gid)
	sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350,pdf=pdf)
	frameProcesses <- list()
	cat("\nCreating frame plots...")
	
	ls <- length(patientID)
	for(i in seq_len(ls))
	{  
	    fid <- patientID[i]
	    tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".jpg", sep=""))
	    jpeg(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
	    print(update(grph, index.cond = list(i)))   
	    dev.off()
	    fGraph <- qaGraph(fileName=tfile, imageDir=idir, pdf=pdf)
	    panelFlag<-tempDist[fid]<=cutoff
	    agTmp <- new("binaryAggregator", passed=panelFlag)
	    frameProcesses[[fid]] <- qaProcessFrame(fid, agTmp,fGraph)
	    cat(".")
	} 
	output<-qaProcess(id=gid, 
			name=paste(channel[1],"  Densityplot",sep=""),
			type=paste(channel[1],"  Densityplot" ,sep=""),
			summaryGraph=sgraph, 
			frameProcesses=frameProcesses)	
	return(output)
}

qaProcess.ecdfOutLiers <- function(flowList,channel=c("FSC-A"),outdir="QAReport",
				   cutoff=0.5,det.dimensions=c(300,300),
                                   func=euc,pdf=FALSE, ...)
{
	if(!is(flowList, "list"))
		stop("'flowList' needs to be of class 'List'")

	if(!is(flowList[[1]], "flowSet"))
		stop("'flowList' needs to be a List of flowSets")

	if(!all(channel %in% colnames(flowList[[1]])))
		stop("Invalid channel")

	if(length(channel)!=1)
	{
		stop("Only one channel is to be specified")
	}
	
	if(!is.character(outdir) || length(outdir)!=1)
		stop("'outdir' must be a valid file path")

	if(!is.numeric(cutoff) || length(cutoff)!=1)
		stop("'cutoff' must be numeric scalar")

	cat("creating summary plots...")
	colScheme<-list()
	gid <- guid()
	tmp <- tempdir()
	tmp <- gsub("\\", "/", tmp, fixed=TRUE)

	alqLen<-length(flowList)
	patientID=sampleNames(flowList[[1]])
	panelFlag=FALSE
	nonOutLierCount=0
	tempgrph<-list()
	tempDist<-list()
	formula<-paste("~","`",channel[1],"`","|","Patient",sep="")
	pointCount<-120
	p<-ppoints(pointCount)
	q<-matrix(ncol=pointCount,nrow=alqLen)
	tempStat<-list()
	
        xmax<-0
        xmin<-100000
        for( i in patientID)
	{   res<-data.frame()
	    for(j in seq_len(alqLen))
	    {
		value=exprs(flowList[[j]]@frames[[i]][,channel])
		newres<-data.frame(Patient=rep(i,nrow(value)),
				  Aliquot=rep(j,nrow(value)),
				  data=value,check.names=FALSE)
		res<-rbind(res,newres)
                frameID<-paste(i,j,sep="_")
		q[j,]<-quantile(value,probs=p,names=FALSE) 

	    }
            tempRange<-range(res[channel])
            if(xmax<tempRange[2])
		xmax<-tempRange[2]
	    if(xmin> tempRange[1])
		xmin<-tempRange[1]                
        
	    tempgrph[[i]]<-ecdfplot(eval(parse(text=formula)),data=res,groups=Aliquot,plot.points=FALSE)
    	    tempStat[[i]]<-sum(euc(q))
	    cat(".")
	}
	xrange<-c(0.9*xmin,1.1*xmax)
	maxres<-range(tempStat)[2]
	nres<-lapply(tempStat,function(x){x/maxres})  
	flag<-rep(TRUE,length(patientID))
	flag[which(nres>cutoff)]<-FALSE

        grph<-ecdfplot(~x|patientID, data = list(patientID = factor(names(tempgrph), levels = names(tempgrph)), 
                   x = seq_along(tempgrph)),
                   xlim = xrange,ylim=c(0,1),
                   xlab=as.character(channel[1]),
                   key=simpleKey(text=as.character(1:alqLen),space="right"),
                   panel = function(x, y, plot.list) {
           do.call(panel.ecdfplot, trellis.panelArgs(tempgrph[[x]], 1))
        })

	sfile <- file.path(tmp, "summary.jpg")    
	jpeg(file=sfile)   
	print(grph)
	dev.off()	
	idir <- file.path(outdir, "images", gid)
	sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)


	frameProcesses <- list()
	cat("\nCreating frame plots...")
	ls <- length(patientID)
	for(i in seq_len(ls))
	{  
		fid <- patientID[i]
		tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".jpg", sep=""))
		jpeg(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
		print(update(grph, index.cond = list(i)))   #obtain individual panels from the multiplanel plots
		dev.off()
		fGraph <- qaGraph(fileName=tfile, imageDir=idir, pdf=pdf)
	        agTmp <- new("binaryAggregator", passed=flag[i])
		frameProcesses[[fid]] <- qaProcessFrame(fid, agTmp,fGraph)
		cat(".")
	} 
	output<-qaProcess(id=gid, 
			name=paste(channel[1],"  ECDF",sep=""),
			type=paste(channel[1]," ECDF" ,sep=""),
			summaryGraph=sgraph, 
			frameProcesses=frameProcesses)
	return(output)
}



## allDens <- function(sets, channel){
##     extDens <- function(x, channel)
##         fsApply(x, function(x)
##             {
##                 tmp <- density(x[,channel])
##                 return(rbind(tmp$x, tmp$y))
##             }, use.exprs=TRUE)
##     if(is.list(sets))
##         dens <- sapply(sets, extDens, channel)
##     else
##         dens <- extDens(sets, channel)
##     return(t(dens))
## }



## ## QA process to compare KL distances
## qaProcess.Similarity <- function(set, channel, outdir, groups=NULL, cutoff=0.1)
## {
##     ## create summary plot and its associated qaGraph object
##     if(length(channel)!=1)
##         stop("'channel' must be of length 1")
##     cat("creating summary plots...")
##     gid <- guid()
##     tmp <- tempdir()
##     sfile <- file.path(tmp, "summary.pdf")
##     dens <- t(allDens(set, channel))
##     col <- if(is.null(grps)) "black" else as.integer(factor(groups))
##     pdf(file=sfile)
##     matplot(dens[seq(1, nrow(dens), by=2), ], dens[seq(2, nrow(dens), by=2), ],
##             type="l", col=col)
##     dev.off()
##     idir <- file.path(outdir, "images", gid)
##     sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350)

##     ## create graphs and aggregators for each frame and wrap in object
##     frameIDs <- sampleNames(set)
##     frameProcesses <- list()
##     cat("\ncreating frame plots...")
##     for(i in 1:length(set)){
##         tfile <- file.path(tmp, paste("frame_", sprintf("%0.2d", i), ".pdf",
##                                       sep=""))
##         pdf(file=tfile)
##         timeLinePlot(set[i], channel)
##         dev.off()
##         ba <- new("binaryAggregator", passed=summary[i]<cutoff)
##         fg <- qaGraph(fileName=tfile, imageDir=idir, width=220)
##         fid <- frameIDs[i]
##         frameProcesses[[fid]] <- qaProcessFrame(fid, ba, fg)
##         cat(".")
##     }

##     ## create qaProcess object
##     cat("\n")
##     return(qaProcess(id=gid, name=name,
##                type="time line", frameIDs=frameIDs, summaryGraph=sgraph,
##                frameProcesses=frameProcesses))
    
## }    







