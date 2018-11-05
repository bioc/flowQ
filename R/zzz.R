## check for a valid ImageMagick installation
.onAttach <- function(...)
{
    vers <- suppressWarnings(system("convert -list configure", intern=TRUE, ignore.stderr=TRUE))
    path <- gsub("Path: *", "", grep("Path:", vers, value=TRUE))[1]
    vers <- gsub("VERSION *", "", grep("^VERSION", vers, value=TRUE))[1]
    if(!length(vers) || !length(path))
        stop("\n\nUnable to find installation of ImageMagick on this system.\n",
             "Please install the ImageMagick library (http://www.imagemagick.org)\n",
             call.=FALSE)
    path <- dirname(dirname(path))
    packageStartupMessage(sprintf("Using ImageMagic library at %s\n(version %s)\n",
                  path, vers))
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.10")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
