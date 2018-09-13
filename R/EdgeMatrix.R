#' A Function for creating a matrix of `edges` or predicted orthologs
#' 
#' EdgeMatrix exists internally in `Catalog`, but `Catalog` can be an intensely time consuming function for very large datasets. `EdgeMatrix` combined with `SubCatalog` allows for the construction of very large complete linkage networks to be passed to the grid, or otherwise parallelized.
#' @param OrthologsObject An object of predicted Orthologs with various integers relating to that prediction
#' @keywords Homology
#' @export
#' @examples
#' Catalog()


EdgeMatrix <- function(OrthologsObject) {
  
  
  if (!is(OrthologsObject, "Orthologs")) {
    stop ("Object is not an Orthologs object.")
  }
  
  ######
  # Create a matrix matrix of edges
  ######
  Size <- nrow(OrthologsObject)
  # MaxUsefulIterations <- Size - 1L # This would really be -2, but the iterator is created beforehand
  FilledPositions <- OrthologsObject[upper.tri(OrthologsObject)]
  TotalRows <- sapply(FilledPositions,
                      function(x) nrow(x),
                      simplify = TRUE)
  EdgeMatrix <- matrix(data = NA_integer_,
                        nrow = sum(TotalRows),
                        ncol = Size)
  
  Count <- 1L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      CurrentRows <- nrow(OrthologsObject[m1, m2][[1]])
      EdgeMatrix[(Count:(Count + CurrentRows - 1L)), m1] <- as.integer(OrthologsObject[m1, m2][[1]][, "QueryGene"])
      EdgeMatrix[(Count:(Count + CurrentRows - 1L)), m2] <- as.integer(OrthologsObject[m1, m2][[1]][, "SubjectGene"])
      Count <- Count + CurrentRows
    }
  }
  return(EdgeMatrix)
}