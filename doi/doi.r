library('stringr')
library('readr')
library('rcrossref')

doiref <- function(path, style = 'apa'){
        mystring <- readr::read_file(path)
        doi <- unlist(stringr::str_extract_all(mystring, "\\b10\\.(\\d+\\.*)+[\\/](([^\\s\\.])+\\.*)+\\b"))
        doi <- unique(doi)
        ref <- vector()
        for (i in 1:length(doi)){
                temp <- try(rcrossref::cr_cn(dois = doi[i], format = "text", style = style), T)
                ref <- c(ref,temp)
        }
        readr::write_lines(ref, path = 'bibliography.txt')
        return(ref)
}