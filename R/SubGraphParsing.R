#' Collect alignments, distance matrices, dendrograms, and phylograms for subgraphs of genes
#'
#' @param SubGraphs A vector either subgraphs from a large adjacency graph, or characters identifying positions on an adjacency graph
#' @param DBPATH A character vector specifying a sqlite database built by DECIPHER
#' @param GeneCalls A named list of dataframes providing columns named "Start" "Stop" "Index" and "Strand"
#' @param CopyType "Single" to select on clusters with no competing nodes from the same genome "Multi" allows for all possible clusters
#' @param InputType "Groups" indicates "Subgraphs" are a list of character vectors "Graphs" indicates a list of subgraphs
#' @param OriginGraph The graph from which the child subraphs were generated
#' @param Verbose Print progress bars to screen
#' @keywords Graph, SubGraph, Synteny, Orthology, Alignment
#' @export
#' @examples
#' ParseSubGraphs()

ParseSubGraphs <- function(SubGraphs,
                           DBPATH,
                           GeneCalls,
                           CopyType = "Single",
                           InputType = "Groups",
                           OriginGraph,
                           Verbose = TRUE) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  Graphs <- vector("list",
                   length = length(SubGraphs))
  
  if (InputType == "Groups") {
    if (Verbose) {
      print("Generating SubGraphs ... ")
    }
    for (m1 in seq_along(SubGraphs)) {
      Graphs[[m1]] <- induced_subgraph(graph = OriginGraph,
                                       vids = SubGraphs[[m1]])
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / length(SubGraphs))
      }
    }
    if (Verbose) {
      cat("\n")
    }
    
  } else if (InputType == "Graphs") {
    Graphs <- SubGraphs
  }
  
  MultiCopy <- sapply(SubGraphs,
                      function(x) any(duplicated(str_extract(string = names(x[1]),
                                                             pattern = "(^[0-9]+)"))))
  if (CopyType == "Single") {
    Graphs <- Graphs[!MultiCopy]
  } else if (CopyType == "Multi") {
    Graphs <- Graphs[MultiCopy]
  } else if (CopyType == "All") {
    Graphs <- Graphs
  }
  
  
  
  KeepSet <- which(sapply(Graphs,
                          function(x) length(names(x[1]))) > 3L)
  ######
  # phylo objects are meant to be used for RF distance, and unrooted
  # phylo objects of size 3 only have 1 node
  ######
  
  Graphs <- Graphs[KeepSet]
  
  SourceStringSets <- vector("list",
                             length = length(GeneCalls))
  
  ReturnList <- vector("list",
                       length = length(Graphs))
  
  for (m1 in seq_along(SourceStringSets)) {
    SourceStringSets[[m1]] <- SearchDB(dbFile = DBPATH,
                                       identifier = names(GeneCalls)[m1],
                                       nameBy = "description",
                                       verbose = FALSE)
  }
  
  print("Generating SubGraph Data ... ")
  
  for (m1 in seq_along(Graphs)) {
    
    ReturnList[[m1]] <- vector("list",
                               length = 9L)
    
    CurrentIDs <- names(Graphs[[m1]][1])
    
    ReturnList[[m1]][[1]] <- Graphs[[m1]]
    
    CurrentDataMat <- do.call(rbind,
                              strsplit(CurrentIDs,
                                       split = "( |_)"))
    CurrentDataMat <- matrix(as.integer(CurrentDataMat),
                             nrow = nrow(CurrentDataMat))
    CurrentIDs <- CurrentIDs[order(CurrentDataMat[, 1L],
                                   CurrentDataMat[, 2L])]
    CurrentDataMat <- CurrentDataMat[order(CurrentDataMat[, 1L],
                                           CurrentDataMat[, 2L]), ]
    Positions <- t(mapply(function(x, y) sapply(GeneCalls[x],
                                                function(z) unlist(z[y, c("Strand", "Start", "Stop")]),
                                                USE.NAMES = FALSE,
                                                simplify = TRUE),
                          x = CurrentDataMat[, 1L],
                          y = CurrentDataMat[, 3L]))
    CurrentDataMat <- cbind(CurrentDataMat,
                            Positions)
    rownames(CurrentDataMat) <- NULL
    UnAlignedSeqs <- unlist(extractAt(x = DNAStringSet(unlist(mapply(function(y, z) SourceStringSets[[y]][[z]],
                                                                     y = CurrentDataMat[, 1L],
                                                                     z = CurrentDataMat[, 2L]))),
                                      at = IRangesList(mapply(function(y, z) IRanges(start = y,
                                                                                     end = z),
                                                              y = CurrentDataMat[, 5L],
                                                              z = CurrentDataMat[, 6L]))))
    UnAlignedSeqs <- DNAStringSet(ifelse(test = CurrentDataMat[, 4L],
                                         yes = reverseComplement(UnAlignedSeqs),
                                         no = UnAlignedSeqs))
    
    names(UnAlignedSeqs) <- CurrentIDs
    
    ReturnList[[m1]][[2]] <- AlignTranslation(myXStringSet = UnAlignedSeqs,
                                              readingFrame = 1L,
                                              verbose = FALSE,
                                              type = "AAStringSet")
    ReturnList[[m1]][[3]] <- AlignTranslation(myXStringSet = UnAlignedSeqs,
                                              readingFrame = 1L,
                                              verbose = FALSE,
                                              type = "DNAStringSet")
    
    ReturnList[[m1]][[4]] <- DistanceMatrix(myXStringSet = ReturnList[[m1]][[2]],
                                            includeTerminalGaps = TRUE,
                                            verbose = FALSE)
    ReturnList[[m1]][[5]] <- DistanceMatrix(myXStringSet = ReturnList[[m1]][[3]],
                                            includeTerminalGaps = TRUE,
                                            correction = "Jukes-Cantor",
                                            verbose = FALSE)
    
    ReturnList[[m1]][[6]] <- IdClusters(myDistMatrix = ReturnList[[m1]][[4]],
                                        myXStringSet = ReturnList[[m1]][[2]],
                                        method = "NJ",
                                        verbose = FALSE,
                                        showPlot = FALSE,
                                        type = "dendrogram")
    TempTree <- "~/TempTree"
    WriteDendrogram(x = ReturnList[[m1]][[6]],
                    file = TempTree,
                    quoteLabels = FALSE)
    ReturnList[[m1]][[7]] <- read.tree(file = TempTree)
    ReturnList[[m1]][[7]]$tip.label <- str_extract(string = ReturnList[[m1]][[7]]$tip.label,
                                                   pattern = "(^[0-9]+)")
    ReturnList[[m1]][[7]] <- unroot(ReturnList[[m1]][[7]])
    system("rm ~/TempTree")
    
    ReturnList[[m1]][[8]] <- IdClusters(myDistMatrix = ReturnList[[m1]][[5]],
                                        myXStringSet = ReturnList[[m1]][[3]],
                                        method = "ML",
                                        verbose = FALSE,
                                        showPlot = FALSE,
                                        type = "dendrogram")
    TempTree <- "~/TempTree"
    WriteDendrogram(x = ReturnList[[m1]][[8]],
                    file = TempTree,
                    quoteLabels = FALSE)
    ReturnList[[m1]][[9]] <- read.tree(file = TempTree)
    ReturnList[[m1]][[9]]$tip.label <- str_extract(string = ReturnList[[m1]][[9]]$tip.label,
                                                   pattern = "(^[0-9]+)")
    ReturnList[[m1]][[9]] <- unroot(ReturnList[[m1]][[9]])
    system("rm ~/TempTree")
    
    names(ReturnList[[m1]]) <- c("Cluster",
                                 "AAAlignment",
                                 "DNAAlignment",
                                 "AADist",
                                 "DNADist",
                                 "AADend",
                                 "AAPhylo",
                                 "DNADend",
                                 "DNAPhylo")
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / length(Graphs))
    }
  }
  if (Verbose) {
    cat("\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  if (length(ReturnList) < 1) {
    ReturnList <- "No clusters fit requirements"
  }
  return(ReturnList)
}

