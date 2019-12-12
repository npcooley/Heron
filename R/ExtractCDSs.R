#' Extract coding regions of a gene to return a DNAStringSet for that gene.

#' @param Match A character string of the format ..X..Y..X.. where X indicates a gap between starts and stops and Y indicates a gap between Starts and Starts.
#' @param Start Gene Start if breaks are given for the coding sequence.
#' @param Stop Gene Stop if no breaks are given for the coding sequence.
#' @param Strand Strand of the coding sequence.
#' @param SourceString DNAString from which to pull the gene.
#' @param Mode Argument for future use to allow multi-exon genes to be extracted in other formats.

ExtractCDSs <- function(Match,
                        Start,
                        Stop,
                        Strand = 0L,
                        SourceString,
                        Mode = "Character") {
  if (grepl(pattern = "Y",
            x = Match)) {
    PosMat <- do.call(rbind,
                      strsplit(x = strsplit(x = Match,
                                            split = "Y",
                                            fixed = TRUE)[[1L]],
                               split = "X",
                               fixed = TRUE))
    Gene <- DNAStringSet(unlist(extractAt(x = SourceString,
                                          at = IRanges(start = as.integer(PosMat[, 1L]),
                                                       end = as.integer(PosMat[, 2L])))))
    if (Strand) {
      Gene <- reverseComplement(Gene)
    }
  } else {
    Gene <- extractAt(x = SourceString,
                      at = IRanges(start = Start,
                                   end = Stop))
    if (Strand) {
      Gene <- reverseComplement(Gene)
    }
  }
  return(Gene)
}
