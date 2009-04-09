## Set up a HTML header
myOpenHtmlPage <- function(name, title = "", path="../") 
{
    name = gsub("\\.html$", "", name)
    con = file(paste(name, ".html", sep = ""), open = "wt")
    writeLines(paste("<html>\n<head>\n<title>", title, "</title>", sep=""),
               con)
    writeLines(paste("<link rel=\"stylesheet\" type=\"text/css\" ",
                     "href=\"", path, "qaReport.css\">", sep=""), con)
    writeLines(paste("<script src=\"", path, "qaReport.js\" ",
                     "type=\"text/javascript\"></script>", sep=""), con)
    writeLines(paste("</head>\n<body style=\"font-family: helvetica,arial",
                     ",sans-serif;\">", sep=""), con)
    return(con)
}

## include links to pdf versions of images if pdf=TRUE
pdfLink <- function(vimg, bimg, class, id, pdf=TRUE)
{
    if(pdf)
        paste("<a href='", vimg,  "' target='QAdetails'>\n<img class='", class,
              "' src='", bimg, "' id='img_", id, "'>\n</a>\n", sep="")
    else
        paste("<img class='", class, "Nolink' src='", bimg, "' id='img_", id, "'>\n",
              sep="")
}


## create HTML report for (lists of) QA processes 
writeQAReport  <- function(set, processes, outdir="./qaReport",
                           grouping=NULL, pagebreaks=TRUE,
                           pdf=TRUE)
{ 
    ## We only need panel tabs if 'set' is a list of 'flowSets'
    single <- FALSE
    if(!is.list(set)){
        set <- list(set)
        single <- TRUE
    }else if(! all(sapply(processes, is.list)) ||
       ! length(set) == length(processes))
        stop("Argument 'processes' must be a list of lists of ",
             "'qaProcess' objects for multiple panels")

    ## For the overview page, we need to match samples across panels
    sID <- all(sapply(set, function(x) "SampleID" %in% colnames(pData(x))))
    if(length(set)>1 && !sID)
        warning("Some of the panels in 'set' don't have a global ",
                "sample identifier.\nUnable to create overview.")
    ## A lot of sanity checking up front
    if(file.exists(file.path(outdir, "index.html")))
        warning("Target directory already exists. Content may be ",
                    "overwritten")
    if(!is(processes, "list") ||
       !all(sapply(processes, function(x) is(x, "qaProcess") ||
                   sapply(x, is, "qaProcess"))))
        stop("'processes' must be a list or list of lists of objects of ",
             "class 'qaProcess'")
    if(!all(sapply(set, is, "flowSet")))
        stop("'set' must be a flow set")

    ## copy infrastructure
    sdir <- system.file("htmlTemplates", package = "flowQ")
    file.copy(dir(sdir, full.names=TRUE), file.path(outdir, "images"),
              overwrite=TRUE)
    
    ## iterate over panels
    for(s in seq_along(set)){
        
        ## rearange set according to grouping if necessary
        grps <- NULL
        if(!is.null(grouping)){
            if(!is.character(grouping) ||
               !grouping %in% varLabels(phenoData(set[[s]])))
                stop("'grouping' must be a factor variable in the flowSet ",
                     "phenoData")
            sord <- order(pData(set[[s]])[, grouping])
            if(!all(sord==seq_len(length(set[[s]]))))
                set[[s]] <- set[[s]][sord]
            grps <- pData(set[[s]])[, grouping]
        }

        ## open a file connection  
        ifile <- if(!sID && s==1) "index" else if(!sID)
            paste("index", s-1, sep="") else  paste("index", s, sep="")
        con <- myOpenHtmlPage(file.path(outdir, ifile), "qatest", "images/")
        
        ## setup of table and table header row
        process <- if(length(set)>1) processes[[s]] else processes
        writeLines("<table class=\"QA\">", con)
        pIDs <- sapply(process, slot, "id")
        pNames <- sapply(process, slot, "name")
        pTypes <- sapply(process, slot, "type")
        relBase <- file.path("images", pIDs)
        sumLinks <- file.path(relBase, sapply(process, function(x)
                                              names(x@summaryGraph)))
        sumVecLinks <- gsub("\\..*$", ".pdf", sumLinks)
        nrAggr <- sapply(process, function(x)
                         length(x@frameProcesses[[1]]@frameAggregators))+1
        th <- paste("<th class=\"QAHeader\" colspan=\"", nrAggr, "\" ",
                    "id=\"", pIDs, "_sumHeader\">\n",
                    "<div class=\"QASumButton\" id=\"", pIDs, "_button",
                    "\" onClick=\"toggleImage('", pIDs, "')\">\n",
                    pNames, "\n</div>\n</th>", sep="")
        esel <- sapply(process, function(x)
                       length(x@summaryGraph@fileNames))==0
        th[esel] <- paste("<th class=\"QAHeader\" colspan=\"",
                          nrAggr[esel], "\" ",
                          "id=\"", pIDs[esel], "_sumHeader\">\n",
                          "<div class=\"QASumButton\" id=\"", pIDs[esel],
                          "_button",
                          "\">\n", pNames[esel], "\n</div>\n</th>", sep="",
                          collapse="\n")
        th <- paste(th, collapse="\n")
        pd <- pData(set[[s]][[1]]@parameters)[,c("name", "desc", "minRange",
                                                 "maxRange")]
        writeLines(paste("<tr class=\"QAHeader\">\n<th class=\"QAHeader\">\n",
                         "<div class=\"QASumButton\" id=\"parameters_button",
                         "\" onClick=\"toggleImage('parms')\">\n",
                         "flow set details\n</div>\n</th>",  
                         "</th>\n", th, "\n</tr>", sep=""), con)
        td <- paste("<td class=\"QASummary\" colspan=\"", nrAggr, "\"",
                    " id=\"", pIDs, "_sumBack\">\n",
                    pdfLink(sumVecLinks, sumLinks, "QASummary", pIDs, pdf=pdf),
                    "</td>", sep="", collapse="\n")
        writeLines(paste("<tr class=\"QASummary\">\n<th class=\"QASummary\">",
                         "<span id=\"img_parms\" style=\"display:none;\">",
                         sep=""), con)
        writeLines(pd, con)
        writeLines(paste("<span>\n</th>\n", td, "\n</tr>", sep=""), con)
        
        frameIDs <- sampleNames(set[[s]])
        classes <- paste("QAFrameHeader",
                         c("Even", "Odd")[(seq_along(frameIDs)%%2)+1], sep="")
        names(classes) <- frameIDs
        fpp <- if(pagebreaks) 14 else 10e20
        lf <- length(frameIDs)
        nrPages <- (lf %/% fpp) +1
        counter <- 1
        page <- 1
        lastGrp <- grps[1]
        for(f in frameIDs){
            showRow <- ifelse(counter>fpp, "none", "table-row")
            pd <- pData(set[[s]])[f,]
            phi <- paste("<span id=\"img_pd_", counter,
                         "\" style=\"display:none;\"",
                         ">", sep="", collapse="\n")
            ## new table row and column header for aggregators
            if(is.null(grouping)){## no grouping is specified
                writeLines(paste("<tr class=\"", classes[f], "\" id=\"frow1_",
                                 counter, "\" style=\"display:",
                                 showRow, ";\">\n<th ",
                                 "class=\"QAFrameHeader\" id=\"",
                                 f, "_sumHeader\">\n<div class=\"QARowButton\"",
                                 " id=\"", f, "_button\" onClick=\"toggle",
                                 "Image('pd_", counter, "')\">\n<div class=\"",
                                 "QAFrameHeaderNr\">", counter, "</div><span ",
                                 "class=\"QAFrameHeaderID\">", f, "</span>",
                                 "\n", phi, sep=""), con)
                writeLines(pd, con)
            }else{##grouping
                caption <-  paste("<div class=\"QAFrameHeaderNr\">",
                                  counter, "</div><span class=\"QAFrameHeader",
                                  "ID\">", f, "</span>",
                                  "<span class=\"QAFrameHeaderGrp\">",
                                  grps[counter], "</span>", sep="")
                if(counter==1 || grps[counter]!=lastGrp){
                    writeLines(paste("<tr class=\"", classes[f],
                                     "\" id=\"frow1_",
                                     counter, "\" style=\"display:", showRow,
                                     ";\">\n<th class=\"QAFrameHeader\" id=\"",
                                     f, "_sumHeader\">\n<div class=\"QAGrp",
                                     "Header\">", grps[counter],
                                     "</div>\n<div ",
                                     "class=\"QARowButton\" id=\"", f,
                                     "_button\" onClick=\"toggleImage(",
                                     "'pd_", counter, "')\">", caption,
                                     "\n</div>\n", phi, sep=""), con)
                    writeLines(pd, con)                      
                }else{  
                    writeLines(paste("<tr class=\"", classes[f],
                                     "\" id=\"frow1_",
                                     counter, "\" style=\"display:",
                                     showRow, ";\">\n<th ",
                                     "class=\"QAFrameHeader\" id=\"",
                                     f, "_sumHeader\">\n<div class=\"QARow",
                                     "Button\" id=\"", f, "_button\" ",
                                     "onClick=\"toggle",
                                     "Image('pd_", counter, "')\">", caption,
                                     "\n</div>\n", phi, sep=""), con)
                    writeLines(pd, con)
                }
            }
            writeLines(paste("</span></th>"), con)
            
            ## aggregators first, each process is one column
            for(p in seq_along(process)){
                
                thisProcess <- process[[p]]@frameProcesses[[f]]
                pid <- thisProcess@id
                sid <- thisProcess@summaryGraph@id
                writeLines(paste("<td class=\"QASumAggr\" id=\"", pid,
                                 "_sumHeader\" align=\"center\">", sep=""), con)
                nrDetAgr <- length(thisProcess@frameAggregators)
                offset <- ""
                ## add trigger to access details
                if((counter==1 || (counter-1) %% fpp == 0) && nrDetAgr>0){
                    writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                     pIDs[p], "_detTriggerIn_", page, "\" ",
                                     "onClick=\"",
                                     "toggleDetails(", nrDetAgr, ", ",
                                     length(frameIDs), ", ", p, 
                                     ", '", pIDs[p], "', ", nrPages, ")\">",
                                     "\n+\n</div>", sep=""), con)
                    writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                     pIDs[p], "_detTriggerOut_", page, "\" ",
                                     "onClick=\"",
                                     "toggleDetails(", nrDetAgr, ", ",
                                     length(frameIDs), ", ", p, 
                                     ", '",  pIDs[p], "', ", nrPages, ")\"",
                                     " style=\"display:none;\">\n&#150;\n</div>",
                                     sep=""), con)
                    offset <- paste(" style=\"position:relative; left:-7px;",
                                    "margin-left:15px; margin-right:0px;",
                                    "z-index:0;\"")
                    page <- page+1
                }
                ## add summary aggregator and link to image if necessary
                if(length(sid)>0)
                    writeLines(paste("<div class=\"QAFrameButton\" ",
                                     "id=\"", pid, "_button\" onClick=\"",
                                     "toggleImage('", sid, "')\"", offset,
                                     ">", sep=""), con)
                else
                    writeLines(paste("<div class=\"QAFrameButtonNoSel\"",
                                     offset,
                                     ">", sep=""), con)
                writeLines(thisProcess@summaryAggregator, con)
                writeLines("</div>\n</td>", con)
                ## add detailed aggregators and links to images if necessary
                for(d in seq_along(thisProcess@frameAggregators)){
                    fname <- names(thisProcess@frameAggregators)[d]
                    fname  <- ifelse(is.null(fname), "", paste(fname,
                                                               "\n<br>\n"))
                    did <- thisProcess@frameGraphs[[d]]@id
                    id <- paste(p, d, counter, sep="_")
                    writeLines(paste("<td class=\"QADetAggr\" id=\"row_", id,
                                     "_1\" align=\"center\">", sep=""), con)
                    if(length(did)>0)
                        writeLines(paste("<div class=\"",
                                         "QADetButton\" ", "id=\"button_",id,
                                         "\" onClick=\"toggleImage('", id,
                                         "')\">", fname, sep=""), con)
                    else
                        writeLines("<div>", con)
                    writeLines(thisProcess@frameAggregators[[d]], con)
                    writeLines("</div>\n</td>", con)
                }## end for d
            }## end for p
            writeLines("</tr>", con)
            
            ## new table row and column header for images
            writeLines(paste("<tr class=\"", classes[f], "\" id=\"frow2_",
                             counter, "\" style=\"display:", showRow,
                             ";\">\n<th ",
                             "class=\"QAFrameHeader\">\n</th>", sep=""), con)
            ## now the images, each process is one column
            for(p in seq_along(process)){
                thisProcess <- process[[p]]@frameProcesses[[f]]
                pid <- thisProcess@id
                sid <- thisProcess@summaryGraph@id
                writeLines(paste("<td class=\"QASumGraph\" id=\"",
                                 pid, "_sumBack\" align=\"center\">",
                                 sep=""), con)
                if(length(sid)>0){
                    sGraph <- file.path(relBase[p],
                                        names(thisProcess@summaryGraph))
                    sVecGraph <- gsub("\\..*$", ".pdf", sGraph)
                    writeLines(pdfLink(sVecGraph, sGraph, "QASumGraph", sid,
                                       pdf=pdf), con)
                }
                writeLines("</td>", con)
                for(d in seq_along(thisProcess@frameAggregators)){
                    did <- thisProcess@frameGraphs[[d]]@id
                    id <- paste(p, d, counter, sep="_")
                    writeLines(paste("<td class=\"QADetGraph\" id=\"row_", id,
                                     "_2\" align=\"center\">", sep=""), con)
                    if(length(did)>0){
                        fGraph <- file.path(relBase[p],
                                            names(thisProcess@frameGraphs[[d]]))
                        fVecGraph <- gsub("\\..*$", ".pdf", fGraph)
                        writeLines(pdfLink(fVecGraph, fGraph, "QADetGraph", id,
                                           pdf=pdf), con)
                    }
                    writeLines("</td>", con)
                }## end for d
            }## end for p
            lastGrp <- grps[counter]
            counter <- counter+1
            writeLines("</tr>", con)
        }## end for f
        writeLines("</table>", con)
        
        ## the page navigation
        writeLines(paste("<div class=\"QAPagesTile\"><table width=\"100%\"",
                         " style=\"padding-right:35px;\"><tr><td align=",
                         "\"left\">", sep=""), con)
        if(lf>fpp){
            from <- c(seq(1, lf, fpp))
            nt <- length(from)
            to <- c(seq(fpp, lf, fpp), lf)[1:nt]
            writeLines(paste("<span class=\"QAPages\" id=\"pages_", 1:nt,
                             "\" onClick=\"togglePages(", from, ", ", to, ", ",
                             nt, ", ", lf, ")\">Frames ", from, "-", to,
                             "</span>", sep=""), con)
        }
        np <- length(set)
        iFiles <- if(!sID && np>1) c("", 1:(np-1)) else as.character(1:np)
        if(np>1){
            writeLines("</td><td align=\"right\">", con)
            panels <- paste("<span class=\"QAPanels\" id=\"panels_", 1:np,
                            "\"n><a class=\"QAPanels\" href=\"index",
                            iFiles, ".html\">Panel ",
                            1:np, "</a></span>", sep="")
            if(!is.null(names(set)))
                panels[s] <- gsub(paste("Panel", s),
                                  paste("Panel ", s, " <i><small>(",
                                        names(set)[s],
                                        ")</i></small>", sep=""),
                                  gsub("QAPanels","QAPanelsAct", panels[s]))
            if(sID)
                panels <- c(paste("<span class=\"QAPanels\" id=\"panels_0",
                                  "\"n><a class=\"QAPanels\" href=\"index",
                                  ".html\">Summary</a></span>", sep=""),
                            panels)
            writeLines(panels, con) 
        }
        writeLines("</td></tr></table></div>", con)  
        closeHtmlPage(con)
    }##end s

    ## We create an overview page if we have multiple panels
    if(sID){
        if(single)
          processes <- list(processes)
        ifile <- "index"
        con <- myOpenHtmlPage(file.path(outdir, ifile), "qatest", "images/")
        on.exit(closeHtmlPage(con))
        summary <- failedProcesses(processes, set)
        writeLines(summary, con)
    }
    return(file.path(outdir, "index.html"))
}



## Create QA output based on a flowSet and a list of QA functions.
## This is a very basic convenience function, for more complex experiments
## including panels use writeQAReport directly
qaReport <- function(set, qaFunctions, outdir="./qaReport", argLists,
                     grouping=NULL, ...)
{
    processes <- list()
    for(i in seq_along(qaFunctions)){
        cat(paste("quality process ", i, ":\n", sep=""))
        if(missing(argLists))
            processes[[i]] <- do.call(qaFunctions[i],
                                      list(set=set, outdir=outdir,
                                           grouping=grouping))
        else{
            argLists[[i]]$set <- set
            argLists[[i]]$outdir <- outdir
             argLists[[i]]$grouping <- grouping
            processes[[i]] <- do.call(qaFunctions[i], argLists[[i]])
        }
        save(processes, file=file.path(outdir, "processes.rda"))
    }
    writeQAReport(set, processes, outdir=outdir, grouping=grouping, ...)
}





## count numbers of failed qaProcesses in a list of lists of such objects.
## Each item in the outer list is a list of qaProcess objects for a single
## panel. The output is an object for which HTML output can be generated
## via a writeLines method.
## The phenoData of the list of flowSets that gets passed as the second
## argument has to contain a column sampleIDs which provides the mapping
## of samples over panels.
failedProcesses <- function(processes, set)
{
    ## some sanity checking first
    sampleIDs <- lapply(set, function(x) pData(x)$SampleID)
    if(!all(listLen(lapply(sampleIDs, unique))==listLen(sampleIDs)))
        stop("'SampleIDs' must be unique in each panel")
    sids <- unlist(sampleIDs)
    comSids <- unique(sids)
    fids <- unlist(lapply(set, function(x) rownames(pData(x))))
    #if(any(duplicated(fids)))
    #    stop("The 'FrameIDs' in the whole experiment are not unique")
    allChannels <- c(unique(unlist(lapply(set, colnames))), "global")
    res <- ranges <- mapping <- vector(length(set), mode="list")
    ## iterate over panels
    for(i in seq_along(processes)){
        fmat <- matrix(0, ncol=length(allChannels), nrow=length(comSids),
                           dimnames=list(comSids, allChannels))
        clist <- slist <- NULL
        nrSum <- 0
        ## iterate over qaProcess objects for one panel
        for(j in seq_along(processes[[i]])){
            nrSamp <- 0
            mlist <- NULL
            ## iterate over samples in the qaProcess object
            for(pro in processes[[i]][[j]]@frameProcesses){
                ## match frameID to sampleID
                samp <- sids[match(pro@frameID, fids)]
                mlist <- rbind(mlist, c(samp, pro@frameID, nrSamp+1))
                ## check for the multiple channels and iterate over those 
                if(length(pro@frameAggregators) && !processes[[i]][[j]]@type %in% c("Density","ECDF","KLD","SummaryStatistic","BoundaryEvents")){
                    channels <- names(pro@frameAggregators)
                    clist <- c(clist, channels)
                    
                    for(chan in seq_along(pro@frameAggregators)){
                        fmat[samp, channels[chan]] <-
                            fmat[samp, channels[chan]] + 
                                as.numeric(!pro@frameAggregators[[chan]]@passed)
                    }
                ## count the "global" qaProcesses    
                }else{
                    fmat[samp, "global"] <- fmat[samp, "global"] +
                        as.numeric(!pro@summaryAggregator@passed)
                    nrSum <- nrSum+1
                }
                slist <- c(slist, samp)
                nrSamp <- nrSamp+1
            }
        }
        colnames(mlist) <- c("sample", "frame", "number")
        res[[i]] <- fmat
        mapping[[i]] <- mlist
        ## We want to sum up over qaProcesses for one panel
        rtemp <- numeric(length(allChannels))
        names(rtemp) <- allChannels
        if(!is.null(clist)){
            tc <- table(clist)/nrSamp
            rtemp[names(tc)] <- tc
        }
        rtemp["global"] <- nrSum/nrSamp
        ranges[[i]] <- rtemp
    }
    ## we also want an overall summary
    sum <- res[[1]]
    if(length(res)>1)
        for(i in 2:length(res))
            sum <- sum+res[[i]]
    names(res) <- names(ranges) <- names(mapping) <- names(set)
    return(new("qaProcessSummary", panels=res, summary=sum, ranges=ranges,
               mapping=mapping))
}
