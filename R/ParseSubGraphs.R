######
# Parse Groups extracted from a Graph
# as subgraphs?
# If it's a valid graph, build the stringset, alignment, and dendrogram
# Origin is an igraph graph
# DBPATH is a character vector pointing to a DECIPHER DB
# GeneCalls are a list of genecalls
# Groups are a list of character vectors pointing to nodes in the graph
######


ParseSubGraphs <- function(Groups,
                           DBPATH,
                           GeneCalls,
                           IgnoreDefault = FALSE,
                           Verbose = TRUE) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    cat("\n")
  }
  
  ReturnList <- vector(mode = "list",
                       length = length(Groups))
  SourceStringSets <- vector(mode = "list",
                             length = length(GeneCalls))
  for (m1 in seq_along(SourceStringSets)) {
    SourceStringSets[[m1]] <- SearchDB(dbFile = DBPATH,
                                       identifier = names(GeneCalls)[m1], # identifier should match genecall names
                                       nameBy = "description", # description should just be the index
                                       verbose = FALSE)
  }
  
  for (m1 in seq_along(ReturnList)) {
    ReturnList[[m1]] <- vector(mode = "list",
                               length = 3L)
    CurrentIDs <- Groups[[m1]]
    IDMat <- do.call(rbind,
                     strsplit(CurrentIDs,
                              split = "( |_)"))
    IDMat <- matrix(as.integer(IDMat),
                    nrow = nrow(IDMat))
    IDMat <- IDMat[order(IDMat[, 1L],
                         IDMat[, 2L],
                         IDMat[, 3L]), ]
    CurrentCalls <- vector(mode = "list",
                           length = nrow(IDMat))
    for (m2 in seq_along(CurrentCalls)) {
      CurrentCalls[[m2]] <- GeneCalls[[IDMat[m2, 1L]]][IDMat[m2, 3L], ]
    }
    CurrentCalls <- do.call(rbind,
                            CurrentCalls)
    CurrentCalls <- cbind(CurrentCalls,
                          "GenomeName" = IDMat[, 1L])
    CurrentGenes <- vector(mode = "list",
                           length = nrow(CurrentCalls))
    for (m2 in seq_along(CurrentGenes)) {
      CurrentGenes[[m2]] <- ExtractCDSs(Match = CurrentCalls$Match[m2],
                                        Start = CurrentCalls$Start[m2],
                                        Stop = CurrentCalls$Stop[m2],
                                        Strand = CurrentCalls$Strand[m2],
                                        SourceString = SourceStringSets[[CurrentCalls$GenomeName[m2]]][[CurrentCalls$Index[m2]]])
    }
    # CurrentGenes <- mapply(function(z, y, x, w, v) ExtractCDSs(Match = v,
    #                                                            Start = w,
    #                                                            Stop = x,
    #                                                            Strand = y,
    #                                                            SourceString = z[[1]]),
    #                        v = CurrentCalls$Match,
    #                        w = CurrentCalls$Start,
    #                        x = CurrentCalls$Stop,
    #                        y = CurrentCalls$Strand,
    #                        z = mapply(function(u, s) SourceStringSets[[u]][[s]],
    #                                   u = IDMat[, 1L],
    #                                   s = IDMat[, 2L]))
    CurrentGenes <- do.call(c,
                            CurrentGenes)
    names(CurrentGenes) <- paste(IDMat[, 1],
                                 IDMat[, 2],
                                 IDMat[, 3],
                                 sep = "_")
    print(CurrentGenes)
    print(CurrentCalls)
    ReturnList[[m1]][[1L]] <- CurrentGenes
    if (all(CurrentCalls$Coding) & !IgnoreDefault) {
      AlignedGenes <- AlignTranslation(myXStringSet = CurrentGenes,
                                       sense = "+",
                                       direction = "5' to 3'",
                                       type = "AAStringSet",
                                       readingFrame = 1L,
                                       verbose = FALSE)
      DistMat <- DistanceMatrix(myXStringSet = AlignedGenes,
                                type = "matrix",
                                includeTerminalGaps = TRUE,
                                verbose = FALSE)
      SetDend <- IdClusters(myDistMatrix = DistMat,
                            myXStringSet = AlignedGenes,
                            method = "NJ",
                            type = "dendrogram",
                            showPlot = FALSE,
                            verbose = FALSE)
    } else {
      AlignedGenes <- AlignSeqs(myXStringSet = CurrentGenes,
                                verbose = FALSE)
      DistMat <- DistanceMatrix(myXStringSet = AlignedGenes,
                                type = "matrix",
                                includeTerminalGaps = TRUE,
                                correction = "Jukes-Cantor",
                                verbose = FALSE)
      SetDend <- IdClusters(myDistMatrix = DistMat,
                            myXStringSet = AlignedGenes,
                            method = "ML",
                            type = "dendrogram",
                            showPlot = FALSE,
                            verbose = FALSE)
    }
    ReturnList[[m1]][[2L]] <- AlignedGenes
    ReturnList[[m1]][[3L]] <- SetDend
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / length(Groups))
    }
  }
  
  if (Verbose) {
    cat("\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
    cat("\n")
  }
  
  return(ReturnList)
}









