#' A function for aligning core genomes.
#' 
#' Core Aligner is a wrapper for `DECIPHER`'s `AlignTranslation` function. Given a subset of the output of `Catalog` or `SubCatalog`, it will align genes from the complete linkage clusters and then concatenate those alignments.
#' @param ListOfSets A list of matrices
#' @param PATH A character string identifying a sqlite database
#' @param GeneCalls A list of dataframes, one for each genome, with named columns "Start", "Stop", and "Strand"
#' @param IDs A vector of characters indicating the identifiers for which genomes to select
#' @keywords Core Genome
#' @export
#' @examples
#' CoreAligner()

CoreAligner <- function(ListOfSets,
                        PATH,
                        GeneCalls,
                        IDs,
                        Verbose = FALSE) {
  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
  }
  TotalGenomes <- ncol(ListOfSets[[1]])
  stopifnot(TotalGenomes == length(GeneCalls))
  Genomes <- vector("list",
                    length = TotalGenomes)
  for (i in seq_along(GeneCalls)) {
    Genomes[[i]] <- SearchDB(dbFile = PATH,
                             tblName = "Seqs",
                             identifier = IDs[i],
                             type = "XStringSet",
                             nameBy = "identifier",
                             limit = -1L,
                             verbose = FALSE)
  }
  ######
  # Use a DNAStringSetList
  ######
  Genomes <- DNAStringSetList(Genomes)
  GeneSets <- vector("list",
                     length = length(ListOfSets))
  if (Verbose == TRUE) {
    pBar <- txtProgressBar(style = 1L)
  }
  for (i in seq_along(ListOfSets)) {
    ######
    # Extract Gene IDs from matrices
    ######
    GeneIDs <- apply(ListOfSets[[i]],
                     2L,
                     function(x) unique(x[!is.na(x)]))
    GeneSets[[i]] <- vector("list",
                            length = TotalGenomes)
    for (j in seq_along(GeneIDs)) {
      ######
      # use gene IDs to collect genes, ReverseComplement when necessary
      ######
      GeneSets[[i]][[j]] <- unlist(extractAt(x = Genomes[[j]][GeneCalls[[j]][GeneIDs[j], "Index"]],
                                             at = IRanges(start = GeneCalls[[j]][GeneIDs[j], "Start"],
                                                          end = GeneCalls[[j]][GeneIDs[j], "Stop"])))
      # GeneSets[[i]][[j]] <- DNAStringSet(Genomes[[j]][GeneCalls[[j]][GeneIDs[j], "Start"]:GeneCalls[[j]][GeneIDs[j], "Stop"]])
      if (GeneCalls[[j]][GeneIDs[j], "Strand"] == 1L) {
        GeneSets[[i]][[j]] <- reverseComplement(GeneSets[[i]][[j]])
      }
    }
    GeneSets[[i]] <- do.call(c,
                             GeneSets[[i]])
    GeneSets[[i]] <- AlignTranslation(myXStringSet = GeneSets[[i]],
                                      readingFrame = 1L,
                                      verbose = FALSE)
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = i/length(GeneSets))
      }
  }
  ######
  # Concatenate into a single DNAStringSet with multiple strings, they will be in the order that
  # the genomes exist in the database
  ######
  CoreGenome <- do.call(xscat, GeneSets)
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(CoreGenome)
}

