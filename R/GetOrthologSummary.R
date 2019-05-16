#' Summarize statistics about ortholog predictions.
#' 
#' @param OrthologsObject An object of class "Orthologs"
#' @param GeneCalls A list of dataframes, one for each genome in SyntenyObject, with named columns "Index", "Start", "Stop", "Strand", and an optional "Annotation" column.
#' @param DBPath A character vector specifying a sqlite database built by DECIPHER
#' @param Verbose Run with progress bar, return total upon completion.
#' @param SimilarityScores A logical indicating whether to get the global alignment scores for all predicted orthologs
#' @param Type provide alignment similarity scores in nucleotide or amino acid space. AAStringSet or DNAStringSet
#' @keywords Orthology, Synteny, Gene Length, Alignment
#' @export
#' @examples
#' GetOrthologSummary()

GetOrthologSummary <- function(OrthologsObject,
                               GeneCalls,
                               DBPath,
                               SimilarityScores = FALSE,
                               Type = "DNAStringSet",
                               Verbose = FALSE) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ######
  # make sure that gene calls are in
  # ascending order
  ######
  
  if (any(lengths(OrthologsObject[lower.tri(OrthologsObject)]) == 0L)) {
    stop ("Orthologs object must not be in 'Sparse format'")
  }
  if (length(GeneCalls) != ncol(OrthologsObject)) {
    stop ("Orthologs object and gene predictions are not compatible")
  }
  if (!("DECIPHER" %in% .packages())) {
    stop ("Required package DECIPHER is not loaded")
  }
  if (!is(OrthologsObject, "Orthologs")) {
    stop ("Object is not an Orthologs object.")
  }
  
  if (SimilarityScores) {
    Genomes <- vector("list",
                        length = length(GeneCalls))
    for (i in seq_along(Genomes)) {
      Genomes[[i]] <- SearchDB(dbFile = DBPath,
                               identifier = names(GeneCalls[i]),
                               nameBy = "identifier",
                               verbose = FALSE)
    }
    Genomes <- DNAStringSetList(Genomes)
  }
  
  Size <- dim(OrthologsObject)[1]
  Total <- (Size^2 - Size) / 2
  
  TotalCoverage <- vector("list",
                          length = Total)
  MinCoverage <- vector("list",
                        length = Total)
  MaxCoverage <- vector("list",
                        length = Total)
  PairMatrix <- vector("list",
                       length = Total)
  IndexMatrix <- vector("list",
                        length = Total)
  QueryGeneLength <- vector("list",
                            length = Total)
  SubjectGeneLength <- vector("list",
                              length = Total)
  CombinedGeneLength <- vector("list",
                               length = Total)
  NormGeneDiff <- vector("list",
                           length = Total)
  AbsStartDelta <- vector("list",
                          length = Total)
  AbsStopDelta <- vector("list",
                         length = Total)
  NormDeltaStart <- vector("list",
                           length = Total)
  NormDeltaStop <- vector("list",
                          length = Total)
  Scores <- vector("list",
                   length = Total)
  QueryCharacter <- vector("list",
                           length = Total)
  SubjectCharacter <- vector("list",
                             length = Total)
  LabelNames <- vector("list",
                       length = Total)
  Count <- 1L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      
      o <- order(GeneCalls[[m1]][, "Index"],
                 GeneCalls[[m1]][, "Start"],
                 decreasing = FALSE)
      GeneCalls[[m1]] <- GeneCalls[[m1]][o, ]
      o <- order(GeneCalls[[m2]][, "Index"],
                 GeneCalls[[m2]][, "Start"],
                 decreasing = FALSE)
      GeneCalls[[m2]] <- GeneCalls[[m2]][o, ]
      ######
      # Find the distance from the closest hit to the start of the gene
      # in nucleotide space
      # for both the subject and the query
      # take the delta of these differences, and convert to a percentage of the combined gene length
      # do the same for stops
      ######
      AbsStartDelta[[Count]] <- abs(OrthologsObject[m2, m1][[1]][, "QueryStartDisplacement"] - OrthologsObject[m2, m1][[1]][, "SubjectStartDisplacement"])
      AbsStopDelta[[Count]] <- abs(OrthologsObject[m2, m1][[1]][, "QueryStopDisplacement"] - OrthologsObject[m2, m1][[1]][, "SubjectStopDisplacement"])
      PairMatrix[[Count]] <- cbind(OrthologsObject[m1, m2][[1]][, "QueryGene"],
                                   OrthologsObject[m1, m2][[1]][, "SubjectGene"])
      IndexMatrix[[Count]] <- cbind(OrthologsObject[m1, m2][[1]][, "QueryIndex"],
                                    OrthologsObject[m1, m2][[1]][, "SubjectIndex"])
      QueryGeneLength[[Count]] <- GeneCalls[[m1]][PairMatrix[[Count]][, 1L], "Stop"] - GeneCalls[[m1]][PairMatrix[[Count]][, 1L], "Start"] + 1L
      SubjectGeneLength[[Count]] <- GeneCalls[[m2]][PairMatrix[[Count]][, 2L], "Stop"] - GeneCalls[[m2]][PairMatrix[[Count]][, 2L], "Start"] + 1L
      CombinedGeneLength[[Count]] <- QueryGeneLength[[Count]] + SubjectGeneLength[[Count]]
      NormGeneDiff[[Count]] <- abs(QueryGeneLength[[Count]] - SubjectGeneLength[[Count]]) / CombinedGeneLength[[Count]]
      NormDeltaStart[[Count]] <- AbsStartDelta[[Count]] / CombinedGeneLength[[Count]]
      NormDeltaStop[[Count]] <- AbsStopDelta[[Count]] / CombinedGeneLength[[Count]]
      TotalCoverage[[Count]] <- (OrthologsObject[m1, m2][[1]][, "ExactOverlap"] * 2L) / CombinedGeneLength[[Count]]
      MinCoverage[[Count]] <- OrthologsObject[m1, m2][[1]][, "ExactOverlap"] / apply(cbind(QueryGeneLength[[Count]],
                                                                                           SubjectGeneLength[[Count]]),
                                                                                     MARGIN = 1L,
                                                                                     FUN = function(x) min(x))
      MaxCoverage[[Count]] <- OrthologsObject[m1, m2][[1]][, "ExactOverlap"] / apply(cbind(QueryGeneLength[[Count]],
                                                                                           SubjectGeneLength[[Count]]),
                                                                                     MARGIN = 1L,
                                                                                     FUN = function(x) max(x))
      QueryCharacter[[Count]] <- vector("character",
                                        length = nrow(OrthologsObject[m1, m2][[1]]))
      SubjectCharacter[[Count]] <- vector("character",
                                          length = nrow(OrthologsObject[m1, m2][[1]]))
      LabelNames[[Count]] <- vector("character",
                                    length = nrow(OrthologsObject[m1, m2][[1]]))
      if (SimilarityScores) {
        Scores[[Count]] <- vector("list",
                                  length = nrow(OrthologsObject[m1, m2][[1]]))
      }
      
      for (i in seq_len(nrow(OrthologsObject[m1, m2][[1]]))) {
        QueryCharacter[[Count]][i] <- paste(m1,
                                            IndexMatrix[[Count]][i, 1L],
                                            PairMatrix[[Count]][i, 1L],
                                            sep = "_")
        SubjectCharacter[[Count]][i] <- paste(m2,
                                              IndexMatrix[[Count]][i, 2L],
                                              PairMatrix[[Count]][i, 2L],
                                              sep = "_")
        LabelNames[[Count]][i] <- paste(QueryCharacter[[Count]][i],
                                        SubjectCharacter[[Count]][i],
                                        sep = " ")
        if (SimilarityScores) {
          Scores[[Count]][[i]] <- unlist(extractAt(x = c(Genomes[[m1]][OrthologsObject[m1, m2][[1]][i, "QueryIndex"]],
                                                         Genomes[[m2]][OrthologsObject[m1, m2][[1]][i, "SubjectIndex"]]),
                                                   at = IRangesList(IRanges(start = GeneCalls[[m1]][OrthologsObject[m1, m2][[1]][i, "QueryGene"], "Start"],
                                                                            end = GeneCalls[[m1]][OrthologsObject[m1, m2][[1]][i, "QueryGene"], "Stop"]),
                                                                    IRanges(start = GeneCalls[[m2]][OrthologsObject[m1, m2][[1]][i, "SubjectGene"], "Start"],
                                                                            end = GeneCalls[[m2]][OrthologsObject[m1, m2][[1]][i, "SubjectGene"], "Stop"]))))
          if (GeneCalls[[m1]][OrthologsObject[m1, m2][[1]][i, "QueryGene"], "Strand"] == 1L) {
            Scores[[Count]][[i]][1] <- reverseComplement(Scores[[Count]][[i]][1])
          }
          if (GeneCalls[[m2]][OrthologsObject[m1, m2][[1]][i, "SubjectGene"], "Strand"] == 1L) {
            Scores[[Count]][[i]][2] <- reverseComplement(Scores[[Count]][[i]][2])
          }
          Scores[[Count]][[i]] <- 1 - DistanceMatrix(myXStringSet = AlignTranslation(myXStringSet = Scores[[Count]][[i]],
                                                                                     sense = "+",
                                                                                     direction = "5' to 3'",
                                                                                     readingFrame = 1L,
                                                                                     verbose = FALSE,
                                                                                     type = Type),
                                                     includeTerminalGaps = TRUE,
                                                     verbose = FALSE)[1, 2]
        } # end similarity scores conditional
      } # end similarity scores loop
      
      ######
      # Go to next matrix position,
      # Assign next list position
      ######
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = Count / Total)
      }
      Count <- Count + 1L
    } # end of columns loop
  } # end of rows loop
  if (SimilarityScores) {
    DF <- data.frame("TotalCoverage" = unlist(TotalCoverage),
                     "MaxCoverage" = unlist(MaxCoverage),
                     "MinCoverage" = unlist(MinCoverage),
                     "NormDeltaStart" = unlist(NormDeltaStart),
                     "NormDeltaStop" = unlist(NormDeltaStop),
                     "NormGeneDiff" = unlist(NormGeneDiff),
                     "Similarity" = unlist(Scores),
                     stringsAsFactors = FALSE)
    rownames(DF) <- unlist(LabelNames)
  } else {
    DF <- data.frame("TotalCoverage" = unlist(TotalCoverage),
                     "MaxCoverage" = unlist(MaxCoverage),
                     "MinCoverage" = unlist(MinCoverage),
                     "NormDeltaStart" = unlist(NormDeltaStart),
                     "NormDeltaStop" = unlist(NormDeltaStop),
                     "NormGeneDiff" = unlist(NormGeneDiff),
                     stringsAsFactors = FALSE)
    rownames(DF) <- unlist(LabelNames)
  }
  
  if (Verbose) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(DF)
}

