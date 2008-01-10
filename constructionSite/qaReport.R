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





writeQAReport  <- function(set, processes, outdir="./qaReport")
{
    if(!is(set, "flowSet"))
        stop("'set' must be a flow set")
    if(!is(processes, "list") || !all(sapply(processes, is, "qaProcess")))
        stop("'processes' must be a list of objects of class 'qaProcess'")
    if(file.exists(file.path(outdir, "index.html")))
        warning("Target directory already exists. Content may be overwritten")

    ## copy infrastructure (FIXME: generalize)
    system(paste("cp ~/Rpacks/flowQ/constructionSite/images/*.png",
                 file.path(outdir, "images")))
    system(paste("cp ~/Rpacks/flowQ/constructionSite/images/*.css",
                 file.path(outdir, "images")))
    system(paste("cp ~/Rpacks/flowQ/constructionSite/images/*.js",
                  file.path(outdir, "images")))

    ## open a file connection
    con <- myOpenHtmlPage(file.path(outdir, "index"), "qatest", "images/")

    ## setup of table and table header row
    writeLines("<table class=\"QA\">", con)
    pIDs <- sapply(processes, slot, "id")
    pNames <- sapply(processes, slot, "name")
    pTypes <- sapply(processes, slot, "type")
    relBase <- file.path("images", pIDs)
    sumLinks <- file.path(relBase, sapply(processes, function(x)
                                          names(x@summaryGraph)))
    sumVecLinks <- gsub("\\..*$", ".pdf", sumLinks)
    nrAggr <- sapply(processes, function(x)
                     length(x@frameProcesses[[1]]@frameAggregators))+1
    th <- paste("<th class=\"QAHeader\" colspan=\"", nrAggr, "\" ",
                "id=\"", pIDs, "_sumHeader\">\n",
                "<div class=\"QASumButton\" id=\"", pIDs, "_button",
                "\" onClick=\"toggleImage('", pIDs, "')\">\n",
                pNames, "\n</div>\n</th>", sep="", collapse="\n")
    writeLines(paste("<tr class=\"QAHeader\">\n<th class=\"QAHeader\">\n",
                     "</th>\n", th, "\n</tr>", sep=""), con)
    td <- paste("<td class=\"QASummary\" colspan=\"", nrAggr, "\"",
                " id=\"", pIDs, "_sumBack\">\n<a href=\"",
                sumVecLinks, "\" target=\"QAdetails\">\n",
                "<img class=\"QASummary\" src=\"", sumLinks, "\"",
                " id=\"img_", pIDs, "\">\n</a>\n</td>", sep="",
                collapse="\n")
    writeLines(paste("<tr class=\"QASummary\">\n<th class=\"QASummary\">",
                     "\n</th>\n", td, "\n</tr>", sep=""), con)
   
    frameIDs <- sampleNames(set)
    classes <- paste("QAFrameHeader",
                     c("Even", "Odd")[(seq_along(frameIDs)%%2)+1], sep="")
    names(classes) <- frameIDs
    counter <- 1
    for(f in frameIDs){
        ## new table row and column header for aggregators
        writeLines(paste("<tr class=\"", classes[f], "\">\n<th ",
                         "class=\"QAFrameHeader\" id=\"",
                         f, "_sumHeader\">\n<div class=\"QARowButton\"",
                         " id=\"", f, "_button\">frame ", f,
                         "\n</div>\n</th>", sep=""), con)
        ## aggregators first, each process is one column
        for(p in seq_along(processes)){
          
            thisProcess <- processes[[p]]@frameProcesses[[f]]
            pid <- thisProcess@id
            sid <- thisProcess@summaryGraph@id
            writeLines(paste("<td class=\"QASumAggr\" id=\"", pid,
                             "_sumHeader\" align=\"center\">", sep=""), con)
            nrDetAgr <- length(thisProcess@frameAggregators)
            offset <- ""
            ## add trigger to access details
            if(counter==1 && nrDetAgr>0){
                writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                 pid, "_detTriggerIn\" onClick=\"",
                                 "toggleDetails(", nrDetAgr, ", ",
                                 length(frameIDs), ", ", p, 
                                 ", this, '", pid, "_detTriggerOut')\">",
                                 "\n+\n</div>", sep=""), con)
                writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                 pid, "_detTriggerOut\" onClick=\"",
                                 "toggleDetails(", nrDetAgr, ", ",
                                 length(frameIDs), ", ", p, 
                                 ", this, '", pid, "_detTriggerIn')\"",
                                 " style=\"display:none;\">\n&#150;\n</div>",
                                 sep=""), con)
                offset <- paste(" style=\"position:relative; left:-7px;",
                                "margin-left:15px; margin-right:0px;",
                                "z-index:0\"")
            }
            ## add summary aggregator and link to image if necessary
            if(length(sid)>0)
                writeLines(paste("<div class=\"QAFrameButton\" ",
                                 "id=\"", pid, "_button\" onClick=\"",
                                 "toggleImage('", sid, "')\"", offset,
                                 ">", sep=""), con)
            else
                writeLines(paste("<div class=\"QAFrameButtonNoSel\"", offset,
                                 ">", sep=""), con)
            writeLines(thisProcess@summaryAggregator, con)
            writeLines("</div>\n</td>", con)
            ## add detailed aggregators and links to images if necessary
            for(d in seq_along(thisProcess@frameAggregators)){
                did <- thisProcess@frameGraphs[[d]]@id
                id <- paste(p, d, counter, sep="_")
                writeLines(paste("<td class=\"QADetAggr\" id=\"row_", id,
                                 "_1\" align=\"center\">", sep=""), con)
                if(length(did)>0)
                    writeLines(paste("<div class=\"",
                                     "QADetButton\" ", "id=\"button_",id,
                                     "\" onClick=\"toggleImage('", id,
                                     "')\">", sep=""), con)
                else
                    writeLines("<span>", con)
                writeLines(thisProcess@frameAggregators[[d]], con)
                writeLines("</div>\n</td>", con)
            }## end for d
        }## end for p
        writeLines("</tr>", con)
        
        ## new table row and column header for images
        writeLines(paste("<tr class=\"", classes[f], "\">\n<th ",
                         "class=\"QAFrameHeader\">\n</th>", sep=""), con)
        ## now the images, each process is one column
        for(p in seq_along(processes)){
            thisProcess <- processes[[p]]@frameProcesses[[f]]
            pid <- thisProcess@id
            sid <- thisProcess@summaryGraph@id
            writeLines(paste("<td class=\"QASumGraph\" id=\"",
                             pid, "_sumBack\" align=\"center\">",
                             sep=""), con)
            if(length(sid)>0){
                sGraph <- file.path(relBase[p], names(thisProcess@summaryGraph))
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
                    writeLines(paste("<a href=\"", fVecGraph, "\" target=\"",
                                     "QAdetails\">\n<img class=\"QADetGraph\" ",
                                     "src=\"", fGraph, "\" id=\"img_", id,
                                     "\">\n</a>", sep=""), con)
                }
                writeLines("</td>", con)
            }## end for d
        }## end for p
        counter <- counter+1
        writeLines("</tr>", con)
    }## end for f
    ## table footer
    writeLines("</table>", con)
    closeHtmlPage(con)
}



## create QA output based on a flowSet and a list of QA functions 
qaReport <- function(set, qaFunctions, outdir="./qaReport", argLists)
{
    processes <- list()
    for(i in seq_along(qaFunctions)){
        cat(paste("quality process ", i, ":\n", sep=""))
        if(missing(argLists))
            processes[[i]] <- do.call(qaFunctions[i],
                                      list(set=set, outdir=outdir))
        else{
            argLists[[i]]$set <- set
            argLists[[i]]$outdir <- outdir
            processes[[i]] <- do.call(qaFunctions[i], argLists[[i]])
        }
    }
    writeQAReport(set, processes, outdir)
}
