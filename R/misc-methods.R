## ===========================================================================
## writeLines for data frames
## ---------------------------------------------------------------------------
## write method to create HTML output
setMethod("writeLines", signature("data.frame", "file", "missing"),
          function(text, con){
              modes <- sapply(text, mode)
              for(i in which(modes == "numeric"))
                  text[,i] <-  round(as.numeric(as.character(text[,i])),3)
              th <- paste("<th class=\"QAParameter\">",
                          colnames(text), "</th>", sep="")
              th[ncol(text)] <- gsub("QAParameter", "QAParameterRight", th[ncol(text)])
              th <- paste("<tr class=\"QAParameter\">\n",
                          paste(th, sep="", collapse="\n"), "\n</tr>", sep="",
                          collapse="\n")
              td <- apply(text, 1, function(x) 
                          paste("<td class=\"QAParameter\">", x, "</td>", sep=""))
              td[ncol(text),] <-
                  gsub("QAParameter", "QAParameterRight", td[ncol(text),])
              td <- paste(apply(td,2,function(x)
                                paste("<tr class=\"QAParameter\">\n",
                                      paste(x, collapse="\n", sep=""), "</tr>", sep="",
                                      collapse="/n")), collapse="\n", sep="")
              paste("<tr class=\"QAParameter\">\n", td, "</tr>", sep="",
                          collapse="/n")
              writeLines(paste("<table class=\"QAParameter\" align=\"center\">\n",
                               th, td, "\n</table>\n", sep=""), con)
          })




## display details about qaProcessSummary and write HTML object
setMethod("show", signature("qaProcessSummary"),
          function(object)
      {
          cat("qaProcess summary for", ncol(object@summary)-1,
              "parameters and", nrow(object@summary),
              "common samples\n") 
      })


## write method to create HTML output
setMethod("writeLines", signature("qaProcessSummary", "file", "missing"),
          function(text, con)
      {
          samples <- rownames(text@summary)
          channels <- colnames(text@summary)
          panels <- length(text@panels)
          ## we want a summary of all failed qaChecks across panels
          sumRanges <- text@ranges[[1]]
          if(panels>1)
              for(r in 2:panels)
                  sumRanges <- sumRanges + text@ranges[[r]]
          ## The bounding table: columns are panels, rows are samples
          writeLines("<table class=\"qaPanelBound\"><tr><td>", con)
          writeLines(paste("\n\n<table class=\"qaPanelSummary\" ",
                           "align=\"center\">"), con)
          writeLines(paste("<tr class=\"qaPanelSummaryEven\">\n<th ",
                           "class=\"qaPanelSummaryRight\">\n</th>"),con)
          if(panels>1)
              writeLines(paste("<th class=\"qaPanelSummaryRight\">\n<span ",
                               "class=\"qaPanelSummary",
                               "\">summary<span>\n</th>\n", sep=""), con)
          writeLines(paste(paste("<th class=\"qaPanelSummaryTop\">",
                                 "<a href=\"index", 1:panels, ".html\">",
                                 "panel ", 1:panels, "\n<br><span class=\"",
                                 "qaPanelSummarySub\">",
                                 shortNames(names(text@panels), n=12),
                                 "</span></a>\n</th>",
                                 sep="", collapse="\n"), "</tr>", sep=""), con)
          ## iterate over rows (samples)
          for(s in seq_along(samples)){
              thisSamp <- samples[s]
              class <- ifelse((s %% 2)==0, "Even", "Odd") 
              ## first the summary over all pannels
             
              writeLines(paste("\n<tr class=\"qaPanelSummary", class,
                               "\">\n<th ", "class=\"qaPanelSummaryRight\">",    
                               thisSamp, "</th>\n", sep=""), con)
              if(panels>1){
                writeLines(paste("<td class=\"qaPanelSummarySum\">",
                                 sep=""), con)
                writeLines(htmlBarplot(text@summary[s,], sumRanges, class=class), con)
                writeLines("</td>", con)
              }
              ## now iterate over each pannel
              for(p in seq_len(panels)){
                  writeLines("<td class=\"qaPanelSummary\">", con)
                  writeLines(htmlBarplot(text@panels[[p]][s,], text@ranges[[p]],
                                         p, match(thisSamp,
                                                  text@mapping[[p]][,"sample"]),
                                         FALSE, class=class), con)
                  writeLines("</td>", con)
              }
              writeLines("</tr>", con)
          }
          writeLines("</table></table></td></tr>", con)
      })


## create a HTML code "barplot" for a single sample
htmlBarplot <- function(data, ranges, panel=NULL, frame=NULL,
                        rownames=TRUE, ignoreChannel="time", class)
{
    sel <- match(tolower(ignoreChannel), tolower(names(data)))
    if(!is.na(sel)){
        data <- data[-sel]
        ranges <- ranges[-sel]
    }
    ld <- length(data)
    maxRange <- max(ranges)
    width <- ifelse(rownames, maxRange*20 + max(nchar(names(data)))*10, 0)
    out <- paste("<table class=\"", ifelse(rownames, "qaBar", "qaSumBar"),
                 "\" align=\"center\" cellspacing=\"0\"",
                 " cellpadding=\"0\" ",
                 ifelse(is.null(panel), "", paste("onclick=\"link2Panel(", panel,
                 ", ", frame, ")\"", sep="")), "width=\"", width, "\">", sep="")
    rows <- character(ld)
    for(i in seq_len(ld))
    {
        failed <- ifelse(is.na(data[i]), 0, data[i])
        passed <- ranges[i] - failed
        empty <- maxRange - ranges[i]
       
        rows[i] <- paste(ifelse(rownames,
                                paste("<tr class=\"qaBar", class, 
                                      "\">\n<th class=\"qaBar\"><b>",
                                      names(data)[i], "</b></th>\n", collapse="", sep=""),
                                ""),
                         paste(rep(paste("<td class=\"qaBarPassed", class,
                                   "\"></td>\n", sep=""), passed), sep="",
                               collapse=""),
                         paste(rep(paste("<td class=\"qaBarFailed", class,
                                   "\"></td>\n", sep=""), failed), sep="",
                               collapse=""),
                         paste(rep(paste("<td class=\"qaBarEmpty", class,
                                   "\"></td>\n", sep=""), empty), sep="",
                               collapse=""),
                         "</tr>", sep="", collapse="")
    }
    out <- c(out, rows, "</table>")
    return(out)
}

## truncate character vectors to 'n' chars for pretty names plotting
shortNames <- function(x, n=13)
{
    if(any(nchar(x)>n))
        x <- paste(sapply(x, substring, 1,n), "...", sep="")
    return(x)
}

