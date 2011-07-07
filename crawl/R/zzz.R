".First.lib" <- function(library, pkgname)
{
    ## Return a list, each element of which is a vector
    ## the first element of the vector is the stuff before the colon in info[[1]]
    ## the second element is the stuff after the colon (can get > 2 elements some
    ## times but ignore)
    info <- strsplit(library(help=pkgname, character.only=TRUE)$info[[1]], "\\:[ ]+")
    ## Go through the list, pulling out the Package, Version and Built strings
    l <- length(info)
    package <- version <- built <- ""
    for (i in 1:l) {
        if(info[[i]][1] == "Package") package <- info[[i]][2]
        if(info[[i]][1] == "Version") version <- info[[i]][2]
        if(info[[i]][1] == "Built") built <- info[[i]][2]
    }
    ## Print these out
    cat(paste("This is", package, version, "\nBuilt:", built, "\n"))
    ## uncomment for fortran/c code
    ## library.dynam("filenameForDll", pkgname)
    library.dynam("crawl", pkgname)
}

".Last.lib" <- function(libpath)
{
library.dynam.unload("crawl", libpath)
cat("\nBye-Bye from crawl\n\n")
return(invisible())
}

