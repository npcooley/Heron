#' Resolve conflicting ortholog predictions
#' 
#' @param SummaryObject A dataframe with statistics of a predicted ortholog pair, the labels of that specific pair represented as the rowname
#' @export
#' @examples
#' ResolveConflicts()

ResolveTesting <- function(SummaryObject,
                           ResolveBy = "Coverage",
                           Verbose = FALSE) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ######
  # Pull labels from the rownames of the summary object
  # Convert them to an integer matrix for easy comparisons
  ######
  
  LabelCharacters <- rownames(SummaryObject)
  
  LabelList <- strsplit(LabelCharacters,
                        "( |_)",
                        fixed = FALSE)
  
  LabelMatrix <- matrix(as.integer(unlist(LabelList)),
                        nrow = nrow(SummaryObject),
                        byrow = TRUE)
  
  SummarizedPairs <- LabelMatrix[, c(1L, 4L)]
  SummarizedPairs <- unique(SummarizedPairs)
  
  ######
  # Create Lists for query side duplicate genes
  # and subject side duplicate genes
  # within each comparison subset
  ######
  
  DuplicatedQueries <- vector("list",
                              length = nrow(SummarizedPairs))
  DuplicatedSubjects <- vector("list",
                               length = nrow(SummarizedPairs))
  
  ######
  # loop through all comparison subsets and find the duplicated indices
  ######
  
  for (i in seq_len(nrow(SummarizedPairs))) {
    SubLabels <- LabelMatrix[which(LabelMatrix[, 1L] == SummarizedPairs[i, 1L] &
                                     LabelMatrix[, 4L] == SummarizedPairs[i, 2L]), ]
    QDuplicates <- SubLabels[duplicated(SubLabels[, 1:3]), 1:3, drop = FALSE]
    SDuplicates <- SubLabels[duplicated(SubLabels[, 4:6]), 4:6, drop = FALSE]
    QDuplicates <- unique(QDuplicates)
    SDuplicates <- unique(SDuplicates)
    
    DuplicatedQueries[[i]] <- vector("list",
                                     length = nrow(QDuplicates))
    DuplicatedSubjects[[i]] <- vector("list",
                                      length = nrow(SDuplicates))
    
    ######
    # for each duplicated index, collect the labels that conflict
    # perform for both query and subject labels that are duplicates
    ######
    
    for (j in seq_len(nrow(QDuplicates))) {
      CurrentPair <- SubLabels[which(apply(SubLabels,
                                           MARGIN = 1L,
                                           function(x) identical(x[1:3], QDuplicates[j, ]))), ]
      
      if (diff(range(CurrentPair[, 5L])) < .Machine$double.eps ^ 0.5) {
        
        ######
        # if the current conflicting labels all occur on the same index
        # this list is further nested for expansion purposes later
        ######
        
        DuplicatedQueries[[i]][[j]] <- list(apply(CurrentPair,
                                                  MARGIN = 1L,
                                                  function(x) paste(x,
                                                                    c("_", "_", " ", "_", "_", ""),
                                                                    sep = "",
                                                                    collapse = "")))
      } else if (!(diff(range(CurrentPair[, 5L])) < .Machine$double.eps ^ 0.5) &
                 length(CurrentPair[, 5L]) > 2L) {
        
        ######
        # if there are more than 2 conflicting labels
        # and they occur on more than 1 index
        ######
        
        ResolveMultiIndex <- vector("list",
                                    length = length(CurrentPair[, 5L]))
        for (k in seq_along(ResolveMultiIndex)) {
          ResolveMultiIndex[[k]] <- which(CurrentPair[, 5L] == CurrentPair[k, 5L])
        }
        
        ######
        # keep only unique sets that are longer than 1L
        ######
        
        ResolveMultiIndex <- unique(ResolveMultiIndex)
        ResolveMultiIndex <- ResolveMultiIndex[which(lengths(ResolveMultiIndex) > 1L)]
        
        ######
        # If labels remain
        ######
        
        if (length(ResolveMultiIndex) > 0L) {
          
          ######
          # This list is further nested to keep all possible outcomes of ResolveMultiIndex
          ######
          
          DuplicatedQueries[[i]][[j]] <- vector("list",
                                                length = length(ResolveMultiIndex))
          for (k in seq_along(ResolveMultiIndex)) {
            DuplicatedQueries[[i]][[j]][[k]] <- apply(CurrentPair[ResolveMultiIndex[[k]], ],
                                                      MARGIN = 1L,
                                                      function(x) paste(x,
                                                                        c("_", "_", " ", "_", "_", ""),
                                                                        sep = "",
                                                                        collapse = ""))
          }
        } # end of resolve multi index conditional
      } # end of CurrentPair conditional
    } # end of query duplicates loop
    
    ######
    # repeat for subject
    ######
    
    for (j in seq_len(nrow(SDuplicates))) {
      CurrentPair <- SubLabels[which(apply(SubLabels,
                                           MARGIN = 1L,
                                           function(x) identical(x[4:6], SDuplicates[j, ]))), ]
      if (diff(range(CurrentPair[, 2L])) < .Machine$double.eps ^ 0.5) {
        DuplicatedSubjects[[i]][[j]] <- list(apply(CurrentPair,
                                                   MARGIN = 1L,
                                                   function(x) paste(x,
                                                                     c("_", "_", " ", "_", "_", ""),
                                                                     sep = "",
                                                                     collapse = "")))
      } else if (!(diff(range(CurrentPair[, 2L])) < .Machine$double.eps ^ 0.5) &
                 length(CurrentPair[, 2L]) > 2L) {
        
        ######
        # if there are more than 2 conflicting labels
        # and they occur on more than 1 index
        ######
        
        ResolveMultiIndex <- vector("list",
                                    length = length(CurrentPair[, 2L]))
        for (k in seq_along(ResolveMultiIndex)) {
          ResolveMultiIndex[[k]] <- which(CurrentPair[, 2L] == CurrentPair[k, 2L])
        }
        
        ######
        # keep only unique sets that are longer than 1L
        ######
        
        ResolveMultiIndex <- unique(ResolveMultiIndex)
        ResolveMultiIndex <- ResolveMultiIndex[which(lengths(ResolveMultiIndex) > 1L)]
        
        
        ######
        # If labels remain
        ######
        
        if (length(ResolveMultiIndex) > 0L) {
          
          ######
          # This list is further nested to keep all possible outcomes of ResolveMultiIndex
          ######
          
          DuplicatedSubjects[[i]][[j]] <- vector("list",
                                                 length = length(ResolveMultiIndex))
          
          for (k in seq_along(ResolveMultiIndex)) {
            DuplicatedSubjects[[i]][[j]][[k]] <- apply(CurrentPair[ResolveMultiIndex[[k]], ],
                                                       MARGIN = 1L,
                                                       function(x) paste(x,
                                                                         c("_", "_", " ", "_", "_", ""),
                                                                         sep = "",
                                                                         collapse = ""))
          }
        } # end of resolve multi index conditional
      } # end of CurrentPair conditional
    } # end of subject duplicates loop
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = i / nrow(SummarizedPairs))
    }
  } # end i loop
  
  ######
  # Combine labels, remove empty label slots, keep only unique label sets
  ######
  
  AllRemovals <- c(DuplicatedQueries,
                   DuplicatedSubjects)
  
  AllRemovals <- unlist(AllRemovals,
                        recursive = FALSE)
  
  AllRemovals <- unlist(AllRemovals,
                        recursive = FALSE)
  
  if (length(which(lengths(AllRemovals) > 0)) > 0) {
    
    ######
    # if there are labels to remove
    # loop through the AllRemovals list
    # which is a list of conflicting labels
    # these labels are character vectors
    ######
    
    AllRemovals <- AllRemovals[which(lengths(AllRemovals) > 0L)]
    AllRemovals <- unique(AllRemovals)
    
    if (Verbose) {
      TimeStop <- Sys.time()
      cat("\n")
      print("Removal Labels Collected")
      print(TimeStop - TimeStart)
      cat("\n")
    }
    
    RemovePositions <- vector("list",
                              length = length(AllRemovals))
    
    if (ResolveBy == "Coverage" |
        ResolveBy == "Similarity") {
      
      ######
      # for these positions, keep the maximum position
      ######
      
      for (i in seq_along(AllRemovals)) {
        
        CurrentSubSet <- SummaryObject[rownames(SummaryObject) %in% AllRemovals[[i]], ]
        SubSetPositions <- which(rownames(SummaryObject) %in% AllRemovals[[i]])
        
        TargetPositions <- CurrentSubSet[, ResolveBy]
        KeepPosition <- which.max(TargetPositions)
        RemovePositions[[i]] <- SubSetPositions[-KeepPosition]
        
        if (Verbose) {
          setTxtProgressBar(pb = pBar,
                            value = i / length(AllRemovals))
        }
      } # end of removals loop
    } else if (ResolveBy == "NormDeltaStart" |
               ResolveBy == "NormDeltaStop" |
               ResolveBy == "NormGeneDiff") {
      
      ######
      # for these metrics, keep the minimum position
      ######
      
      for (i in seq_along(AllRemovals)) {
        
        CurrentSubSet <- SummaryObject[rownames(SummaryObject) %in% AllRemovals[[i]], ]
        SubSetPositions <- which(rownames(SummaryObject) %in% AllRemovals[[i]])
        
        TargetPositions <- CurrentSubSet[, ResolveBy]
        KeepPosition <- which.min(TargetPositions)
        RemovePositions[[i]] <- SubSetPositions[-KeepPosition]
        
        if (Verbose) {
          setTxtProgressBar(pb = pBar,
                            value = i / length(AllRemovals))
        }
      } # end of removals loop
    } # end of ResolveBy conditional
    
    SubSetSummary <- SummaryObject[-(unique(unlist(RemovePositions))), ]
    
    if (Verbose) {
      TimeStop <- Sys.time()
      cat("\n")
      print("Labels Removed")
      print(TimeStop - TimeStart)
      cat("\n")
    }
    
    return(SubSetSummary)
    
  } else {
    
    ######
    # if no indices to remove
    ######
    
    if (Verbose) {
      TimeStop <- Sys.time()
      cat("\n")
      print(TimeStop - TimeStart)
      print("No Conflicts Present")
    }
    
    return(SummaryObject)
  }
  
}






