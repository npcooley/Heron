#' A function for doing simple network analysis on sets of homologs
#' 
#' The function `Catalog` is not parallelizable, this is a parallelizable version `Catalog`, but requires an EdgeMatrix instead of an Orthology Object. Lists of these subcatalogs can then be stitched together.
#' @param EdgeMatrix A matrix of predicted orthologs, with each row being an "edge" consisting of two nodes, each being genes in a predicted orthologous pair. The nodes are located in a column that specifies a genome.
#' @param SubDivision Selected rows to build an edge network for
#' @keywords Homology
#' @export
#' @examples
#' SubCatalog()

SubCatalog <- function(EdgeMatrix,
                       SubDivision,
                       Verbose = FALSE) {
  ######
  # Determine the linkages for paired homologs,
  # but not for the full set of homologs,
  # this function is specifically designed to send to
  # the grid to complete a homologs list in parallel, instead of in serial
  ######
  
  ######
  # Set Progress Bar
  ######
  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  ######
  # Cycler Function
  ######
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
  # Upfront work:
  # What are the max number of useful iterations
  # What is max size of the a "core" genome
  # These are cut-offs for the while loop
  ######
  Size <- ncol(EdgeMatrix)
  CoreSize <- ((Size - 1L) * Size) / 2L
  MaxUsefulIterations <- Size - 1L # This would really be -2, but the iterator is created beforehand
  
  GeneSets <- vector("list",
                     length = length(SubDivision))
  for (z1 in seq_along(GeneSets)) {
    # Set While Loop Conditions
    Remain <- TRUE
    # Grab starting matrix
    CurrentMatrix <- EdgeMatrix[SubDivision[z1],
                                 ,
                                 drop = FALSE]
    # Start Infinite Loop Counter
    Counter <- 1L
    # Set starting column for loops
    CurrentOriginColumn <- which(!is.na(EdgeMatrix[SubDivision[z1],
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
        RowsToAdd[[z3]] <- which(EdgeMatrix[, CurrentColumns[z3]] == CurrentIDs[z3])
      }
      RowsToAdd <- sort(unique(unlist(RowsToAdd)))
      RowCycle <- Cycler(Origin = SubDivision[z1],
                         Vector = RowsToAdd)
      CurrentMatrix <- EdgeMatrix[RowsToAdd[RowCycle],
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
    } # end of while loop
    o <- order(as.integer(rownames(CurrentMatrix)))
    GeneSets[[z1]] <- CurrentMatrix[o,
                                    ,
                                    drop = FALSE]
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = z1/length(GeneSets))
    }
  } # end of for loop
  GeneSets <- unique(GeneSets)
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(GeneSets)
}