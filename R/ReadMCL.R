#' Read an MCL output file into R.
#' 
#' @param FileName A character vector specifying the location of an MCL file to read.
#' @keywords MCL, Orthology
#' @export
#' @examples
#' ReadMCL()

ReadMCL <- function(FileName,
                    Verbose = FALSE) {
  ######
  # This function is currently only being tested on a single MCL output format
  # it is currently not extensible and does not have error checking
  ######
  
  if (Verbose) {
    pBar <- txtProgressBar(style = 1L)
    TimeStart <- Sys.time()
  }
  
  RawInput <- readLines(FileName)
  Clusters <- sapply(RawInput,
                     function(x) unlist(strsplit(x, "\t")),
                     USE.NAMES = FALSE,
                     simplify = FALSE)
  ExpandedPairs <- vector("list",
                          length = length(Clusters))
  for (i in seq_along(ExpandedPairs)) {
    if (length(Clusters[[i]]) > 2L) {
      Count <- 1L
      Total <- (length(Clusters[[i]])^2L - length(Clusters[[i]])) / 2L
      o <- order(as.integer(str_extract(Clusters[[i]],
                                        "^(\\d)+")))
      Clusters[[i]] <- Clusters[[i]][o]
      ExpandedPairs[[i]] <- vector("list",
                                   length = Total)
      for (m1 in seq_len(length(Clusters[[i]]) - 1L)) {
        for (m2 in (m1 + 1L):length(Clusters[[i]])) {
          ExpandedPairs[[i]][[Count]] <- c(Clusters[[i]][m1],
                                           Clusters[[i]][m2])
          Count <- Count + 1L
        }
      }
    } else if (length(Clusters[[i]] <= 2L)) {
      ExpandedPairs[[i]] <- list(Clusters[[i]])
    }
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = i / length(ExpandedPairs))
    }
  }
  CollapsedPairs <- unlist(ExpandedPairs,
                           recursive = FALSE)
  AllPairs <- vector("character",
                     length = length(CollapsedPairs))
  if (Verbose) {
    cat("\n")
  }
  
  for (i in seq_along(CollapsedPairs)) {
    AllPairs[i] <- paste(CollapsedPairs[[i]], collapse = " ")
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = i / length(CollapsedPairs))
    }
  }
  if (Verbose) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  
  return(AllPairs)
}