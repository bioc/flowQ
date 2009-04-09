## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address.
guid <- function()
    format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"])


locateParameter<-function(flowList,parm,flowSetIndx,flowFrameIndx){ 

    if(length(parm)!=1){
        stop("Only one parameter is to be specified")
    }
    temp <- pData(parameters(flowList[[flowSetIndx]][[flowFrameIndx]]))
    mIndx <- is.na(temp["desc"])
        temp["desc"][mIndx] <- temp["name"][mIndx]
        nameIndx <- temp["desc"]==parm
	if(length(temp["name"][nameIndx])>1){
	 stop("Multiple channels in a flowFrame stained for the same cell type")
        }
        return(temp["name"][nameIndx])   
}

locateDuplicatedParameters<-function(flowList){
	
    alqLen<- length(flowList)
    names <-data.frame()
    for( i in seq_len(alqLen)){
        temp <- pData(parameters(flowList[[i]][[1]]))	
        mIndx <- is.na(temp["desc"])
        temp["desc"][mIndx] <- temp["name"][mIndx]
        names<-rbind(names,temp[1:nrow(temp)-1,"desc",drop=FALSE])
    }
    dupes<- names[duplicated(names[,1]),1] 
    dIndx <- !duplicated(dupes)
    dupes<-dupes[dIndx]
    return(as.character(dupes))
	
}

normalizeSets <- function(flowList,dupes,peaks=NULL)
{   patientID<-sampleNames(flowList[[1]])
    alqLen<-length(flowList)
    for (cellType  in dupes){
        for( i in patientID){
            inList<-list()   
            alqTable <- rep(FALSE,alqLen)
            frameIndx<-1
            for(j in seq_len(alqLen)){
                        
                parm <- locateParameter(flowList,cellType,j,i)
                if(length(parm)!=0){
    
                    value<-flowList[[j]][[i]][,parm]
                    colnames(value) <- cellType
                    inList[[frameIndx]] <- value            
                    frameIndx <- frameIndx + 1
                    alqTable[j] <- TRUE
                }
            }
                  
            inSet <- flowSet(inList)
            inFrame<-as(inSet,"flowFrame")                
            cur1<-curv1Filter(cellType, filterId="myCurv1Filter")
	    if(is.null(peaks))	peakCount<- length(filter(inFrame,cur1)@filterDetails$myCurv1Filter$boundaries)
	    else  peakCount <- peaks

            norm1 <- normalization(normFun=function(x, parameters, ...) warpSet(x, parameters,peakNr=peakCount, ...),
                            parameters=as.character(cellType), arguments=list(grouping=NULL),
                            normalizationId="Norm")
            nData <- normalize(inSet,norm1)
            frameIndx<-1
            for(m in which(alqTable==T)){
                parm <- locateParameter(flowList,cellType,m,i)
                exprs(flowList[[m]]@frames[[i]])[,parm] <-
                                      exprs(nData[[frameIndx]])
                frameIndx <- frameIndx +1
            }
        
      }
  }
  return(flowList)
}

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



qaProcess.BoundaryPlot <- function(
flowList,
dyes=NULL,
outdir="QAReport",
cutoff=3,
det.dimensions=c(400,400),
pdf=FALSE,...)
{
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
    alqLen<- length(flowList)
    patientID=sampleNames(flowList[[1]])
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    parLbl<-vector(mode="character",length=alqLen)
    legend <-vector(mode="character",length=alqLen)
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    }else{
	dupes <- as.character(dyes)
    }
    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    panelFlag<-list()

    for (cellType in dupes ){
        res<-data.frame()
        panelFlag[[cellType]] <-list()
   	for( i in patientID){
            perc<-matrix(nrow=alqLen,ncol=1)
            panelFlag[[cellType]][[i]] <- TRUE
            for(j in seq_len(alqLen)){
                    par <- locateParameter(flowList,cellType,j,i)
                    if(length(par)!=0){
                        legend[j] <-"green"
                        parLbl[j] <- paste(j," ",par)
                        value=exprs(flowList[[j]]@frames[[i]][,par])
                        colnames(value) <- cellType
          		ranges <- t(range(flowList[[j]]@frames[[i]][,par]))
                        eps<- c(.Machine$double.eps, -.Machine$double.eps)
                        ranges<- ranges+ eps   
                        ef <- char2ExpressionFilter(
                              paste("`", par, "`<=", ranges[1]," | `",par,"`>=",ranges[2], sep="",
                          collapse=""), filterId=cellType)
                        ff <-filter(flowList[[j]]@frames[[i]][,par],ef)
                        perc[j,]<- summary(ff)$p*100                           
                    
                    }else{
                        legend[j] <-"white"
                        parLbl[j] <- paste(j," ")
                    }           
  
            } 
            colnames(perc)<-cellType
            #legend=rep("green",length(perc))
            legend[perc>cutoff]<-"red"
            panelFlag[[cellType]][[i]] <-all(perc[!is.na(perc)]<=cutoff)
            newres<-data.frame(Patient=rep(i,alqLen),Aliquot=seq_len(alqLen),
                              passed=factor(c(perc<=cutoff),levels=c(TRUE,FALSE)),
                              data=perc,check.names=F)
            res<-rbind(res,newres)      
            formula <- paste("`","Aliquot","`","~","`",cellType,"`",sep="")      
            tempgrph[[cellType]][[i]]<-
                         barchart( eval(parse(text=formula)), 
                                data=newres,origin = 0,
                                col=myCol[unique(newres[,"Aliquot"])],
                                key=list(space="right",points=list(pch=19,col=legend),
                                text=list(parLbl),col=myCol) 
                                    )                                                           
            cat(".")
	}
     
    	sfile <- file.path(tmp, paste("summary_", cellType, ".png", sep=""))
        png(file=sfile, width=det.dimensions[1]*2,height=det.dimensions[2]*2)
        formula <- paste("`","Aliquot","`","~","`",cellType,"`","|","Patient",sep="")   
        print(barchart( eval(parse(text=formula)), 
                data=res,col=(c("green","red")),origin = 0,drop.unused.levels=F,
		main="Percentage of margin events",
                groups=res[,"passed"],
 		key=simpleKey(text=as.character(c("failed","passed")),
			space="right",points=FALSE,col=c("red","green")))
        )
    	
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
	
    }
    sfile <- paste(tmp, "summary.png", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, 
		width=max(det.dimensions[1],det.dimensions[2]*lp),pdf=pdf)
	
    frameProcesses <- list()
    cat("\nCreating frame plots...")

    for(i in seq_len(ls)) ##over patient
    {  
	fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp)){  ##over nrow(dyes)
	      tfile <- file.path(tmp, paste("frame_", sprintf("%0.2s", i), "_",
                                          gsub("\\..*$", "", j), ".png",
                                          sep=""))
	      png(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
	      print(tempgrph[[j]][[patientID[i]]])
              dev.off()
	      fnames <- c(fnames, tfile)
	      agTmp[[j]] <- new("binaryAggregator", passed=panelFlag[[j]][[patientID[i]]])
              #   agTmp[[j]] <- new("discreteAggregator", passed=panelFlag[[j]][patientID[i]],x=1)
	      names(agTmp[[j]]) <-paste(dyes[j])
              cat(".")
	}
        names(agTmp) <-dupes
	nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail[!is.na(nfail)])==1) factor(2) else factor(0)
     	    if(sum(nfail[!is.na(nfail)])==0)
             val <- factor(1)
        ba <- new("discreteAggregator", x=val)
	fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
				  width=min(det.dimensions[1], lp*det.dimensions[2]), pdf=pdf)
	fid <- patientID[i]
	frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
						    summaryAggregator=ba,
						    frameAggregators=agTmp,
						    frameGraphs=fGraphs)
        cat(".")
    }
    cat("\n")
	output<-qaProcess(id=gid, 
		        name="Margin events ",
			#type=paste(dyes,"  ", " Margin events" ,sep=""),
                        type="BoundaryEvents",
			summaryGraph=sgraph, 
			frameProcesses=frameProcesses)
	return(output)
}

qaProcess.2DStatsPlot <- function(
		flowList,
		dyes=c("FSC-A","SSC-A"),
		outdir="QAReport",
		outBound=0.25, # numeric between 0/1 indicating outlier boundary defualt 0.25
		func=mean,	
		det.dimensions=c(400,400),
		pdf=FALSE,...
){
    if(is(dyes,"character")){
            if(length(dyes)!=2)
                    dyes<- matrix(dyes,byrow=TRUE,ncol=2)
            else
                    dyes <- matrix(dyes,ncol=2)
    }

    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
    alqLen<- length(flowList)
    patientID=sampleNames(flowList[[1]])
    ls <- length(patientID)
    lp <- nrow(dyes)
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    tempgrph<-list()
    parLbl<-vector(mode="character",length=alqLen)
    panelFlag<-list()
	
    for(cellType in seq_len(nrow(dyes))){
        tempgrph[[cellType]] <- list()
        panelFlag[[cellType]]<-list()
        outRes<-data.frame()
        for( i in patientID){
            res<-data.frame()
            panelFlag[[cellType]][i]<-TRUE            
            for(j in seq_len(alqLen)){
                par1 <- locateParameter(flowList,dyes[cellType,1],j,i)
                par2 <- locateParameter(flowList,dyes[cellType,2],j,i)
                
                if(length(par1)!=0 && length(par2)!=0){
                        parLbl[j] <- paste(j," ",par1,"/",par2)
                        value1 <- func(exprs(flowList[[j]]@frames[[i]][,par1]))
                        value2 <- func(exprs(flowList[[j]]@frames[[i]][,par2]))
                        newres<-data.frame(Aliquot=j,x=value1,y=value2, Patient=i,check.names=FALSE)
                        res<-rbind(res,newres)
                }else{
                        parLbl[j] <- paste(j," ") 
                }
            }
	    colnames(res) <-c("Aliquot",paste(dyes[cellType,1]),paste(dyes[cellType,2]),"Patient")
            if(nrow(res)==1){
		stop("\n Parameters",paste(colnames(res)[2:3],"should occur only in pairs in
                 in more than one Aliquot in the dataset \n"))
            }
	    tempIndx <- which(pcout(x=as.matrix(res[,2:3]),outbound=outBound)$wfinal01==0)
            if(length(tempIndx)==0){
                indx<-NA
	    }else{
                indx<-tempIndx
            }
            outLier <- rep(FALSE,nrow(res))
            legend <-rep("white",alqLen)
            legend[res[,"Aliquot"]] <-"green"
            if(!is.na(indx)[1]){
                outLier[indx] <- TRUE
                legend[res["Aliquot"][indx,]] <- "red" 
                panelFlag[[cellType]][[i]]<-FALSE
            }
            outRes <- rbind(outRes,cbind(res,outLier))
            formula<-paste("`",dyes[cellType,1],"`"," ", "~"," ","`",
			       dyes[cellType,2],"`","|","Patient",sep="")
            tempgrph[[cellType]][[i]] <- 
			xyplot(eval(parse(text=formula)),data=res,
                            groups=Aliquot,
                            ylim=extendrange(res[,dyes[cellType,1]],f=3),
                            xlim= extendrange(res[,dyes[cellType,2]],f=3),
                            col=myCol[unique(res[,"Aliquot"])],
                            key=list(space="right",points=list(pch=19,col=legend),
				text=list(parLbl),col=myCol),pch=19,cex=2
            		    )
            cat(".")
        }
        
        formula<-paste("`",dyes[cellType,1],"`"," ", "~"," ","`",dyes[cellType,2],"`","|","Patient",sep="")
        sfile <- file.path(tmp, paste("summary_", paste(dyes[cellType,1],"_",dyes[cellType,2],sep=""), ".png", sep=""))
        png(file=sfile, width=det.dimensions[1]*1.5,height=det.dimensions[2]*1.5)
        print(xyplot(eval(parse(text=formula)),data=outRes, 
                                        auto.key=list(space="right"),
                                        groups=outLier,
  					ylim=extendrange(outRes[,dyes[cellType,1]],f=1.8),
                            		xlim= extendrange(outRes[,dyes[cellType,2]],f=1.8),
                                        col=c("green","red"),
                                        key=simpleKey(text=as.character(c("OutLier","Non-outlier")),space="right",points=F,col=c("red","green")),
                                        plot.points=FALSE
                        )
        )
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
    }
    sfile <- paste(tmp, "summary.png", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                                    lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, 
		width=max(det.dimensions[1],det.dimensions[2]*lp),pdf=pdf)
    frameProcesses <- list()
    cat("\nCreating frame plots...")
    
    for(i in seq_len(ls)) ##over patient
    {  
        fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp)){  ##over nrow(dyes)
            tfile <- file.path(tmp, paste("frame_", sprintf("%0.2s", i), "_",
                                            gsub("\\..*$", "", j), ".png",
                                            sep=""))
            png(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
            print(tempgrph[[j]][[patientID[i]]])
            dev.off()
            fnames <- c(fnames, tfile)
            agTmp[[j]] <- new("binaryAggregator", passed=panelFlag[[j]][[patientID[i]]])
            names(agTmp[[j]]) <-paste(dyes[j,1],"/",dyes[j,2],sep="")
            cat(".")
        }
        dyeNames<-apply(dyes,1,function(x){
                                paste(" ",x[1]," / ",x[2]," ",sep="")
                        })
        
        names(agTmp) <-dyeNames
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
                val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
                        width=min(det.dimensions[1], lp*det.dimensions[2]), pdf=pdf)
        fid <- patientID[i]
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                        summaryAggregator=ba,
                        frameAggregators=agTmp,
                        frameGraphs=fGraphs)
        cat(".")
    } 
    cat("\n")	
    
    output<-qaProcess(id=gid, 
                    name=paste("Comparision of  ",substitute(func)," values",sep=""),
                    type="SummaryStatistic",
                    summaryGraph=sgraph, 
                    frameProcesses=frameProcesses)
    return(output)
}

qaProcess.DensityPlot <- function(
	flowList,
 	dyes=NULL,
	outdir="QAReport",
	alpha=0.05,
	det.dimensions=c(400,400),
       	pdf=FALSE,
	...
){
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
  
    alqLen<- length(flowList)
    patientID=sampleNames(flowList[[1]])
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    }else{
	dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL

    for (cellType in dupes ){
        formula<-paste("~","`",cellType,"`","|","Patient",sep="")
        tempgrph[[cellType]]<-list()
        tempDist[[cellType]]<-list()
        tempInput <-list()
        parLbl<-vector(mode="character",length=alqLen)
        ymax<-0
       
        for( i in patientID){
            res<-data.frame()
            tempStat<-matrix(ncol=256,nrow=alqLen)
            for(j in seq_len(alqLen)){
                par <- locateParameter(flowList,cellType,j,i)
                if(length(par)!=0){
                    parLbl[j] <- paste(j," ",par)
                    value=exprs(flowList[[j]]@frames[[i]][,par])
                    colnames(value) <- cellType
                    newres<-data.frame(Patient=rep(i,nrow(value)),
                                      Aliquot=rep(j,nrow(value)),
                                      data=value,check.names=FALSE)
                    res<-rbind(res,newres)
                    valRange<-range(flowList[[j]]@frames[[i]][,par])
                    tempStat[j,]<- density(value,n=256,
				    from=valRange[1,],
				    to=valRange[2,])$y
                    tempInput[[j]] <-value
                    yrng<-range(tempStat[j,])[2]
                    if(yrng>ymax)
                            ymax<-yrng
                }else{
                    parLbl[j] <- paste(j," ") 
                    tempInput[[j]] <-NA
                }
            }
            
	    dst <- KLdist.matrix(tempInput,symmetrize=TRUE)
            tempDist[[cellType]][[i]] <- sum(dst,na.rm=T)/
	                             length(which(!is.na(dst)==T))
            tempgrph[[cellType]][[i]]<-
                        densityplot(eval(parse(text=formula)),
                            data=res,
			    groups=Aliquot,plot.points=FALSE,
                            col=myCol[unique(res[,"Aliquot"])],
                            key=simpleKey(text=parLbl,space="right",
				    points=F,col=myCol),
                            lwd=2)
            cat(".")
            }
  
        xrange<-c(0.9*valRange[1,],1.1*valRange[2,])
        sfile <- file.path(tmp, paste("summary_", cellType, 
				".png", sep=""))
        png(file=sfile, width=det.dimensions[1]*1.5,
			height=det.dimensions[2]*1.5)
        print(densityplot(~x|patientID, 
              data = list(patientID = 
		      factor(names(tempgrph[[cellType]]),
              levels = names(tempgrph[[cellType]])),
		      x = seq_along(tempgrph[[cellType]])), 
              xlim =xrange,ylim=c(0,1.05*ymax),
              xlab=as.character(cellType),
              key=simpleKey(text=parLbl,space="right",
		      points=F,col=myCol),
              lwd=2,
              panel = function(x, y, plot.list) {
                  do.call(panel.densityplot, 
			  trellis.panelArgs(
				  tempgrph[[cellType]][[x]], 1))
                  }
              ))
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
   }
  
    sfile <- paste(tmp, "summary.png", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "),
			    " -geometry +0+0 -tile ",
                 lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, 
		    width=max(det.dimensions[1],
			    det.dimensions[2]*lp), pdf=pdf)

    frameProcesses <- list()
    cat("\ncreating frame plots...")

    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
         threshFlag[[i]]<-rep(TRUE,ls)
         tmpVal <- unlist(tempDist[[dupes[i]]])
         tmpIndx <- which(tmpVal < mean(tmpVal) )
         threshFlag[[i]][calout.detect(unlist(tempDist[[dupes[i]]]),
			 alpha=alpha,method="GESD")$ind]<-FALSE       
          threshFlag[[i]][tmpIndx] <- TRUE
    }

    for(i in seq_len(ls)){ #over patient
	fnames <- NULL
        agTmp <- aggregatorList()
	for(j in seq_len(lp)){   #over dupes
	    tfile <- file.path(tmp, paste("frame_", sprintf("%0.2s", i), "_",
                                          gsub("\\..*$", "", j), ".png",
                                          sep=""))
	    png(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
	    print(tempgrph[[dupes[j]]][[patientID[i]]])
	    dev.off()
	    fnames <- c(fnames, tfile)
	    agTmp[[j]] <- new("numericAggregator", passed=threshFlag[[j]][i],
                              x=tempDist[[dupes[j]]][[patientID[i]]]/
			      max(unlist(tempDist[[dupes[j]]])))
            cat(".")
	}
    
	names(agTmp) <- dupes
	nfail <- !sapply(agTmp, slot, "passed")
            val <- if(sum(nfail)==1) factor(2) else factor(0)
     	    if(sum(nfail)==0)
             val <- factor(1)

	ba <- new("discreteAggregator", x=val)
	fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
				  width=min(det.dimensions[1],
					  lp*det.dimensions[2]), pdf=pdf)
	fid <- patientID[i]
	frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
						    summaryAggregator=ba,
						    frameAggregators=agTmp,
						    frameGraphs=fGraphs)
    }
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="Density plot",
                     type="Density", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}

qaProcess.ECDFPlot <- function(flowList,
	dyes=NULL,
	outdir="QAReport",
	alpha = 0.05,
	det.dimensions=c(400,400),
	pdf=FALSE,...
){
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
    alqLen<- length(flowList)
    patientID=sampleNames(flowList[[1]])
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    }else{
	dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
# pointCount<-256
#   p<-ppoints(pointCount)
#   q<-matrix(ncol=pointCount,nrow=alqLen)

    for (cellType in dupes ){

	formula<-paste("~","`",cellType,"`","|","Patient",sep="")
	tempgrph[[cellType]]<-list()
	tempDist[[cellType]]<-list()
	xmax<-0
	xmin<-100000
        parLbl<-vector(mode="character",length=alqLen)
        tempInput <- list() 
	for( i in patientID){

	    res<-data.frame()
# tempStat<-matrix(ncol=pointCount,nrow=alqLen)
	    for(j in seq_len(alqLen)){
   
		par <- locateParameter(flowList,cellType,j,i)
		if(length(par)!=0){
    		    parLbl[j] <- paste(j," ",par)
		    value=exprs(flowList[[j]]@frames[[i]][,par])
                    colnames(value) <- cellType
		    newres<-data.frame(Patient=rep(i,nrow(value)),
				      Aliquot=rep(j,nrow(value)),
				      data=value,check.names=FALSE)
		    res<-rbind(res,newres)
		    valRange<-range(flowList[[j]]@frames[[i]][,par])
#	    tempStat[j,]<-quantile(value,probs=p,names=FALSE) 
		    tempInput[[j]] <-value
		}else{
#			tempStat[j,]=NA
			tempInput[[j]] <- NA
                        parLbl[j] <- paste(j," ") 
                }
	    }
  
	    tempgrph[[cellType]][[i]]<-
                         ecdfplot(eval(parse(text=formula)),
				data=res,groups=Aliquot,
                                plot.points=FALSE,
				col=myCol[unique(res[,"Aliquot"])],
				key=simpleKey(text=parLbl,
					space="right",
					points=F,
					col=myCol)
                                 )
            dst <- KLdist.matrix(tempInput,symmetrize=TRUE)
            if(length(which(!is.na(dst)==T))!=0)
	         tempDist[[cellType]][[i]]<-sum(dst,na.rm=T)/
		                          length(which(!is.na(dst)==T))
	    else     
		 tempDist[[cellType]][[i]]<- NA

	    cat(".")
	}

        xrange<-c(0.9*valRange[1,],1.1*valRange[2,])
	sfile <- file.path(tmp, paste("summary_", cellType, ".png", sep=""))
        png(file=sfile, width=det.dimensions[1],height=det.dimensions[2])
        print(ecdfplot(~x|patientID, 
		        data = list(patientID = factor(names(tempgrph[[cellType]]),
			            levels = names(tempgrph[[cellType]])),
				    x = seq_along(tempgrph[[cellType]])), 
			xlim =xrange,ylim=c(0,1.05),
			xlab=as.character(cellType),
			key=simpleKey(text=parLbl,space="right",
				points=F,col=myCol),
			panel = function(x, y, plot.list){
			      do.call(panel.ecdfplot, 
				      trellis.panelArgs(
					      tempgrph[[cellType]][[x]], 1))
			      }
			  ))
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
	
    }

    sfile <- paste(tmp, "summary.png", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, 
		      width=max(det.dimensions[1],
			      det.dimensions[2]*lp), pdf=pdf)

    frameProcesses <- list()
    cat("\ncreating frame plots...")
    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
	threshFlag[[i]]<-rep(TRUE,ls)
        tmpVal <- unlist(tempDist[[dupes[i]]])
        tmpIndx <- which(tmpVal < mean(tmpVal) )
        threshFlag[[i]][calout.detect(unlist(tempDist[[dupes[i]]]),
                         alpha=alpha,method="GESD")$ind]<-FALSE
        threshFlag[[i]][tmpIndx] <- TRUE
    }								      

    for(i in seq_len(ls)){ #over patient
            fnames <- NULL
            agTmp <- aggregatorList()
	    for(j in seq_len(lp)){   #over dupes
        
                tfile <- file.path(tmp, paste("frame_", sprintf("%0.2s", i), "_",
                                            gsub("\\..*$", "", j), ".png",
                                            sep=""))
                png(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
                print(tempgrph[[dupes[j]]][[patientID[i]]])
                dev.off()
                fnames <- c(fnames, tfile)
                agTmp[[j]] <- new("numericAggregator", 
                                   passed=threshFlag[[j]][i],
                                   x= tempDist[[dupes[j]]][[patientID[i]]]/
				              max(unlist(tempDist[[dupes[j]]])))
                cat(".")
            }
        
            names(agTmp) <- dupes
  	    nfail <- !sapply(agTmp, slot, "passed")
            val <- if(sum(nfail)==1) factor(2) else factor(0)
     	    if(sum(nfail)==0)
             val <- factor(1)
          
            ba <- new("discreteAggregator", x=val)
            fGraphs <- qaGraphList(imageFiles=fnames, 
			    imageDir=idir,
                            width=min(det.dimensions[1], 
			    lp*det.dimensions[2]), pdf=pdf)
            fid <- patientID[i]
            frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                      summaryAggregator=ba,
                                                      frameAggregators=agTmp,
                                                      frameGraphs=fGraphs)
    }
    cat("\n")
    return(qaProcess(id=gid, name="ECDF Plot",
                     type="ECDF", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}

qaProcess.KLDistPlot <- function(
	flowList,
	dyes=NULL,
	outdir="QAReport",
	alpha=0.05,
	det.dimensions=c(400,400),
 	pdf=FALSE, ...
){
    cat("creating summary plots...")
    gid <- guid()
    tmp <- tempdir()
    tmp <- gsub("\\", "/", tmp, fixed=TRUE)
    sfiles <- NULL
  
    alqLen<- length(flowList)
    patientID=sampleNames(flowList[[1]])
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    }else{
	dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    colorFun<-colorRampPalette(c("yellow","red"))
	
    for (cellType in dupes ){
     
	formula<-paste("~","`",cellType,"`","|","Patient",sep="")
	tempgrph[[cellType]]<-list()
	tempDist[[cellType]]<-list()
        parLbl<-vector(mode="character",length=alqLen)
        outRes<-data.frame()
	for( i in patientID){

                res<-data.frame()
                tempList<- list()
                for(j in seq_len(alqLen)){
                        par <- locateParameter(flowList,cellType,j,i)
                        if(length(par)!=0){
   			    parLbl[j] <- paste(j," ",par)
                            value=exprs(flowList[[j]]@frames[[i]][,par])
                            colnames(value) <- cellType
                            tempList[[j]]<-value
                        }
			else {
				tempList[[j]]<-NA
                                parLbl[j] <- paste(j," ")
			}
                }
	
	dst <- KLdist.matrix(tempList,symmetrize=TRUE) 
	tempDist[[cellType]][[i]]<-sum(dst,na.rm=T)/length(which(!is.na(dst)==T))
   
        pm<-as.matrix(dst)
	diag(pm)<-NA
        z<-as.matrix(as.vector(pm),ncol=1)
	 
        x <- as.character(sapply(parLbl,function(x){
           rep(x,length(parLbl))
        }))
        y <- rep(parLbl,length(parLbl))
        res <- data.frame(x=factor(x,levels=parLbl),y=factor(y,levels=parLbl),z=z)
	tempgrph[[cellType]][[i]]<-levelplot(z~x*y,data=res,xlab="Aliquot",ylab="Aliquot",
                                           main=cellType,scales = list(x = list(rot = 90)),
                                           col.regions=colorFun,colorkey=list(col=colorFun))    
    	Patient=rep(i,nrow(res))
        tm<-cbind(res,Patient)
        outRes<-rbind(tm,outRes)
        cat(".")
	}
        sfile <- file.path(tmp, paste("summary_", cellType, ".png", sep=""))
        png(file=sfile, width=det.dimensions[1]*2,height=det.dimensions[2]*2)
        print(grph<-levelplot(z~x*y|Patient,data=outRes,xlab="Aliquot",ylab="Aliquot",
                               main=cellType,scales = list(x = list(rot = 90)),
                               col.regions=colorFun,colorkey=list(col=colorFun)))
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
    }

    sfile <- paste(tmp, "summary.png", sep="/")
    system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 lp, "x1 ", sfile, sep=""))
    idir <- file.path(outdir, "images", gid)
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=max(det.dimensions[1],det.dimensions[2]*lp), pdf=pdf)

    frameProcesses <- list()
    cat("\ncreating frame plots...")

    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
       threshFlag[[i]]<-rep(TRUE,ls)
       tmpVal <- unlist(tempDist[[dupes[i]]])
       tmpIndx <- which(tmpVal < mean(tmpVal) )
       threshFlag[[i]][calout.detect(unlist(tempDist[[dupes[i]]]),alpha=alpha,method="GESD")$ind]<-FALSE       
       threshFlag[[i]][tmpIndx] <- TRUE
    }

    for(i in seq_len(ls)){ #over patient
	fnames <- NULL
        agTmp <- aggregatorList()
	for(j in seq_len(lp)){   #over dupes	
	    tfile <- file.path(tmp, paste("frame_", sprintf("%0.2s", i), "_",
                                          gsub("\\..*$", "", j), ".png",
                                          sep=""))
	    png(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
	    print(tempgrph[[dupes[j]]][[patientID[i]]])
	    dev.off()
	    fnames <- c(fnames, tfile)
	    agTmp[[j]] <- new("numericAggregator", passed=threshFlag[[j]][i],
                              x=tempDist[[dupes[j]]][[patientID[i]]]/max(unlist(tempDist[[dupes[j]]])))
            cat(".")
	}
    
	names(agTmp) <- dupes
	nfail <- !sapply(agTmp, slot, "passed")
            val <- if(sum(nfail)==1) factor(2) else factor(0)
     	    if(sum(nfail)==0)
             val <- factor(1)

	ba <- new("discreteAggregator", x=val)
	fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
				  width=min(det.dimensions[1], lp*det.dimensions[2]), pdf=pdf)
	fid <- patientID[i]
	frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
						    summaryAggregator=ba,
						    frameAggregators=agTmp,
						    frameGraphs=fGraphs)
    }
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="KL distance plot",
                     type="KLD", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))

}

