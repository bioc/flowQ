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


## create HTML report for (lists of) QA processes 
writeQAReport  <- function(set, processes, outdir="./qaReport",
                           grouping=NULL, pagebreaks=TRUE)
{
    if(!is.list(set))
        set <- list(set)
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
        file.copy(dir(sdir, full.names=TRUE), file.path(outdir, "images"))
    
    ## iterate over panels
    for(s in seq_along(set)){
       
        ## rearange set according to grouping
        grps <- NULL
        if(!is.null(grouping)){
            if(!is.character(grouping) ||
               !grouping %in% varLabels(phenoData(set[[s]])))
                stop("'grouping' must be a factor variable in the flowSet ",
                     "phenoData")
            set[[s]] <- set[[s]][order(pData(set[[s]])[, grouping])]
            grps <- pData(set[[s]])[, grouping]
        }
        
        
        ## open a file connection
        ifile <- ifelse(s==1, "index", paste("index", s, sep=""))
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
        esel <- sapply(process, function(x)  length(x@summaryGraph@fileNames))==0
        th[esel] <- paste("<th class=\"QAHeader\" colspan=\"", nrAggr[esel], "\" ",
                        "id=\"", pIDs[esel], "_sumHeader\">\n",
                        "<div class=\"QASumButton\" id=\"", pIDs[esel], "_button",
                        "\">\n", pNames[esel], "\n</div>\n</th>", sep="",
                          collapse="\n")
        th <- paste(th, collapse="\n")
        pd <- pData(set[[s]][[1]]@parameters)[,c("name", "desc", "minRange",
                                            "maxRange")]
        fh <- paste("<td class=\"QAParameter\">", pd[,1], "</td>\n",
                    "<td class=\"QAParameter\">", pd[,2], "</td>\n",
                    "<td class=\"QAParameter\">", pd[,3], " - ",
                    signif(pd[,4],2),
                    "</td>\n</tr>", sep="", collapse="\n")
        writeLines(paste("<tr class=\"QAHeader\">\n<th class=\"QAHeader\">\n",
                         "<div class=\"QASumButton\" id=\"parameters_button",
                         "\" onClick=\"toggleImage('parms')\">\n",
                         "flow set details\n</div>\n</th>",  
                         "</th>\n", th, "\n</tr>", sep=""), con)
        td <- paste("<td class=\"QASummary\" colspan=\"", nrAggr, "\"",
                    " id=\"", pIDs, "_sumBack\">\n<a href=\"",
                    sumVecLinks, "\" target=\"QAdetails\">\n",
                    "<img class=\"QASummary\" src=\"", sumLinks, "\"",
                    " id=\"img_", pIDs, "\">\n</a>\n</td>", sep="",
                    collapse="\n")
        writeLines(paste("<tr class=\"QASummary\">\n<th class=\"QASummary\">",
                         "<span id=\"img_parms\" style=\"display:none;\">",
                         "<table class=\"QAParameter\" align=\"center\">",
                         "<tr>\n<th class=\"QAParameter\">channel</th>\n",
                         "<th class=\"QAParameter\">fluorochrome</th>\n",
                         "<th class=\"QAParameter\">range</th>\n</tr>\n",
                         fh, "</table><span>\n</th>\n", td, "\n</tr>",
                         sep=""), con)
        
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
            ph <- paste("<span id=\"img_pd_", counter,
                        "\" style=\"display:none;\"",
                        "><table class=\"",
                        classes[f], "Pheno\">",
                        "<tr>", paste("<th class=\"QAPheno\">", names(pd),
                                      "</th>", sep="", collapse="\n"), "</tr>",
                        "<tr>", paste("<td class=\"QAPheno\">", pd, "</td>",
                                      sep="", collapse="\n"), "</tr>\n</table>",
                        "</span>", sep="", collapse="\n")
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
                                 "\n", ph, "</th>", sep=""), con)
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
                                     "\n</div>\n", ph, "\n</th>", sep=""), con)
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
                                     "\n</div>\n", ph, "</th>", sep=""), con)
                }
            }
            
            
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
                    writeLines(paste("<a href=\"", sVecGraph, "\" target=\"",
                                     "QAdetails\">\n<img class=\"QASumGraph\" ",
                                     "src=\"", sGraph, "\" id=\"img_", sid,
                                     "\">\n</a>", sep=""), con)
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
                        writeLines(paste("<a href=\"",
                                         fVecGraph, "\" target=\"",
                                   "QAdetails\">\n<img class=\"QADetGraph\" ",
                                         "src=\"", fGraph, "\" id=\"img_", id,
                                         "\">\n</a>", sep=""), con)
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
                         " style=\"padding-right:35px;\"><tr><td align=\"left\">", sep=""), con)
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
        if(np>1){
            writeLines("</td><td align=\"right\">", con)
            writeLines(paste("<span class=\"QAPanels\" id=\"panels_1,",
                             "\"n><a href=\"index.html\">Panel 1</a>",
                             "</span>", sep=""), con) 
            writeLines(paste("<span class=\"QAPanels\" id=\"panels_", 2:np,
                             "\"n><a href=\"index", 2:np, ".html\">Panel ",
                             2:np, "</a></span>", sep=""), con)
        }
        writeLines("</td></tr></table></div>", con)  
        closeHtmlPage(con)
    }##end s
}



## Create QA output based on a flowSet and a list of QA functions.
## This is a very basic convenience function, for more complex experiments
## including panels use writeQAReport directly
qaReport <- function(set, qaFunctions, outdir="./qaReport", argLists,
                     grouping=NULL)
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
    }
    writeQAReport(set, processes, outdir)
}



