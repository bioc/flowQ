## ===========================================================================
## writeLines for data frames
## ---------------------------------------------------------------------------
## write method to create HTML output
setMethod("writeLines", signature("data.frame", "file", "missing"),
          function(text, con){
              modes <- sapply(text, mode)
              for(i in which(modes == "numeric"))
                  text[,i] <-  round(as.numeric(as.character(text[,i])),3)
              th <- paste("<th>",
                          colnames(text), "</th>", sep="")
              th[ncol(text)] <- gsub("<th", "<th class=\"right\"", th[ncol(text)])
              th <- paste("<tr>\n",
                          paste(th, sep="", collapse="\n"), "\n</tr>\n", sep="",
                          collapse="\n")
              td <- apply(text, 1, function(x) 
                          paste("<td>", x, "</td>", sep=""))
              td[ncol(text),] <-
                  gsub("<td", "<td class=\"right\"", td[ncol(text),])
              td <- paste(apply(td,2,function(x)
                                paste("<tr>\n",
                                      paste(x, collapse="\n", sep=""), "\n</tr>", sep="",
                                      collapse="/n")), collapse="\n", sep="")
              paste("<tr>\n", td, "</tr>\n", sep="",
                    collapse="/n")
              writeLines(paste("<table class=\"QAParameter\" align=\"center\">\n",
                               th, td, "\n</table>", sep=""), con)
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
          writeLines("<table class=\"qaPanelBound\">\n<tr>\n<td>\n", con)
          writeLines(paste("<table class=\"qaPanelSummary\">\n<tr",
                           "class=\"even\">\n<th class=\"left\">\n</th>"), con)
          if(panels>1)
              writeLines("<th class=\"left\"><span>Summary</span></th>", con)
          writeLines(paste(paste("<th class=\"right\"><a href=\"index", 1:panels, ".html\">",
                                 "panel ", 1:panels, "\n<br><span>",
                                 shortNames(names(text@panels), n=12), "</span></a>",
                                 "\n</th>", sep="", collapse="\n"), "\n</tr>", sep=""), con)
          ## iterate over rows (samples)
          for(s in seq_along(samples))
          {
              thisSamp <- samples[s]
              class <- ifelse((s %% 2)==0, "even", "odd") 
              ## first the summary over all pannels
              writeLines(paste("\n<tr class=\"", class, "\">\n<th class=\"left\">",
                               thisSamp, "</th>", sep=""), con)
              if(panels>1)
              {
                  writeLines("\n<td class=\"sum\">", con)
                  writeLines(htmlBarplot(text@summary[s,], sumRanges), con)
                  writeLines("</td>", con)
              }
              ## now iterate over each pannel
              for(p in seq_len(panels)){
                  writeLines("\n<td class=\"bars\">", con)
                  writeLines(htmlBarplot(text@panels[[p]][s,], text@ranges[[p]],
                                         p, match(thisSamp,
                                                  text@mapping[[p]][,"sample"]),
                                         FALSE), con)
                  writeLines("</td>", con)
              }
              writeLines("</tr>", con)
          }
          writeLines("</table>\n</td>\n</tr>\n</table>", con)
      })


## create a HTML code "barplot" for a single sample
htmlBarplot <- function(data, ranges, panel=NULL, frame=NULL,
                        rownames=TRUE, ignoreChannel="time")
{
    sel <- match(tolower(ignoreChannel), tolower(names(data)))
    if(!is.na(sel)){
        data <- data[-sel]
        ranges <- ranges[-sel]
    }
    ld <- length(data)
    maxRange <- max(ranges)
    width <- ifelse(rownames, maxRange*20 + max(nchar(names(data)))*10, 0)
    out <- paste("\n<table class=\"qaBar\" align=\"center\"",
                 ifelse(is.null(panel), "", paste("onclick=\"link2Panel(", panel,
                                                  ", ", frame, ")\"", sep="")), ">", sep="")
    rows <- character(ld)
    for(i in seq_len(ld))
    {
        failed <- ifelse(is.na(data[i]), 0, data[i])
        passed <- ranges[i] - failed
        empty <- maxRange - ranges[i]
        lclass <- if(rownames) "" else "link "
        bclass <- if(ranges[i]==0) rep("", maxRange) else
        if(maxRange==1 && ranges[i]==1) "single " else
        c("left ", rep("", maxRange-2), "right ")
        td <-  c(rep(sprintf("<td class=\"%%s%spassed\"></td>\n", lclass), passed),
                 rep(sprintf("<td class=\"%%s%sfailed\"></td>\n", lclass), failed),
                 rep("<td class=\"%sempty\"></td>\n", empty))
        rn <- ifelse(rownames, paste("<th>", names(data)[i], "</th>", collapse="",
                                    sep=""), "")
        rows[i] <- sprintf("<tr>\n%s%s</tr>", rn,
                           paste(mapply(sprintf, td, bclass), sep="", collapse=""))
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

