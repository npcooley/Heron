#' A function for doing simple network analysis on sets of homologs.
#' 
#' Catalog takes in an object of type `Orthology` and builds a list of all complete linkage clusters. This list is in the form of matrices that denote how transitive the predicted clusters are.
#' @param OrthologyObject a matrix of lists, the upper triange of which is filled with matrices of putative homologs
#' @keywords Homology
#' @export
#' @examples
#' Catalog()

Catalog <- function(OrthologyObject,
                    Verbose = FALSE) {
  stopifnot(nrow(OrthologyObject) >= 3L)
  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  ######
  # Function to Cycle correctly through rows and columns
  # Depending on which part of the vector is the actual origin
  Cycler <- function(Origin,
                     Vector) {
    V.Length <- length(Vector)
    O.Position <- which(Vector == Origin)
    if (O.Position == 1L) {
      Cycle <- 1L:length(Vector)
    } else if (length(Vector) > 2L) {
      Cycle <- c(O.Position:V.Length, 1:(O.Position - 1L))
    } else if (length(Vector) == 2L &
               O.Position == 2L) {
      Cycle <- c(2, 1)
    }
    return(Cycle)
  }
  ######
  
  ######
  # UpFront Math
  Size <- nrow(OrthologyObject)
  CoreSize <- ((Size - 1L) * Size) / 2L
  MaxUsefulIterations <- Size - 1L # This would really be -2, but the iterator is created beforehand
  FilledPositions <- OrthologyObject[upper.tri(OrthologyObject)]
  TotalRows <- sapply(FilledPositions,
                      function(x) nrow(x),
                      simplify = TRUE)
  PairsMatrix <- matrix(data = NA_integer_,
                        nrow = sum(TotalRows),
                        ncol = Size)
  ######
  
  ######
  # Create a matrix matrix of edges
  Count <- 1L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      CurrentRows <- nrow(OrthologyObject[m1, m2][[1]])
      PairsMatrix[(Count:(Count + CurrentRows - 1L)), m1] <- as.integer(OrthologyObject[m1, m2][[1]][, 1])
      PairsMatrix[(Count:(Count + CurrentRows - 1L)), m2] <- as.integer(OrthologyObject[m1, m2][[1]][, 2])
      Count <- Count + CurrentRows
    }
  }
  ######
  
  ######
  # Begin work
  GeneSets <- vector("list",
                     length = nrow(PairsMatrix))
  for (z1 in seq_along(GeneSets)) {
    # Set While Loop Conditions
    Remain <- TRUE
    # Grab starting matrix
    CurrentMatrix <- PairsMatrix[z1,
                                 ,
                                 drop = FALSE]
    # Start Infinite Loop Counter
    Counter <- 1L
    # Set starting column for loops
    CurrentOriginColumn <- which(!is.na(PairsMatrix[z1,
                                                    ,
                                                    drop = TRUE]))[1L]
    while (Remain) {
      # Set Starting MatchVector
      CurrentMatchVector <- apply(CurrentMatrix,
                                  MARGIN = 2L,
                                  function(x) if (all(is.na(x))) {
                                    NA
                                  } else {
                                    unique(x)[!is.na(unique(x))]
                                  })
      CurrentColumns <- which(!is.na(CurrentMatchVector))
      CurrentIDs <- CurrentMatchVector[CurrentColumns]
      
      PreviousMatrix <- CurrentMatrix
      RowsToAdd <- vector("list",
                          length = length(CurrentColumns))
      for (z3 in seq_along(RowsToAdd)) {
        RowsToAdd[[z3]] <- which(PairsMatrix[, CurrentColumns[z3]] == CurrentIDs[z3])
      }
      RowsToAdd <- sort(unique(unlist(RowsToAdd)))
      RowCycle <- Cycler(Origin = z1,
                         Vector = RowsToAdd)
      CurrentMatrix <- PairsMatrix[RowsToAdd[RowCycle],
                                   ,
                                   drop = FALSE]
      rownames(CurrentMatrix) <- RowsToAdd[RowCycle]
      ColCycle <- Cycler(Origin = CurrentOriginColumn,
                         Vector = seq_len(Size))
      for (z2 in seq_along(ColCycle)) {
        Column <- CurrentMatrix[, ColCycle[z2]]
        RemoveRow <- integer(length = 0L)
        if (all(is.na(Column))) {
          # do nothing
        } else if (length(!is.na(Column)) > 1L) {
          if (is.na(CurrentMatchVector[ColCycle[z2]])) {
            # If the current column does not have an ID to match to assign one
            CorrectID <- Column[which(!is.na(Column))[1L]]
          } else {
            # If it does, use that ID
            CorrectID <- CurrentMatchVector[ColCycle[z2]]
          }
          RemoveRow <- which(Column != CorrectID &
                               !is.na(Column))
          if (length(RemoveRow) != 0L) {
            CurrentMatrix <- CurrentMatrix[-RemoveRow,
                                           ,
                                           drop = FALSE]
          } # end of removal conditional
        } # end of column check conditional
      } # end of loop through columns
      Counter <- Counter + 1L
      if (nrow(CurrentMatrix) == CoreSize) {
        Remain <- FALSE
      } else if (all(dim(CurrentMatrix) == dim(PreviousMatrix))) {
        CurrentCheck <- CurrentMatrix[!is.na(CurrentMatrix)]
        PreviousCheck <- PreviousMatrix[!is.na(PreviousMatrix)]
        if (all(CurrentCheck == PreviousCheck)) {
          Remain <- FALSE
        }
      } else if (Counter >= MaxUsefulIterations) {
        Remain <- FALSE
      }
    } # End While Loop
    o <- order(as.integer(rownames(CurrentMatrix)))
    GeneSets[[z1]] <- CurrentMatrix[o,
                                    ,
                                    drop = FALSE]
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = z1/length(GeneSets))
    }
  }
  GeneSets <- unique(GeneSets)
  ######
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(GeneSets)
}











