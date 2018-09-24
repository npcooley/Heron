#' A function for parsing genomic General Format Files directly from the NCBI ftp server
#' 
#' This parser is specifically built for gff files that have been produced by NCBI's annotation pipeline
#' @param GFFAddress a vector of characters for ftp addresses of gff.gz files
#' @keywords GFF
#' @export
#' @examples
#' GFFParser()


NCBIGFFParser <- function(GFFAddress,
                          Verbose = FALSE) {
  
  ######
  # error checking
  ######
  
  if (!("stringr" %in% .packages())) {
    stop ("Required package stringr is not loaded.")
  }
  if (is.null(names(GFFAddress))) {
    stop ("Addresses provided must be named.")
  }
  if (any(is.na(names(GFFAddress)))) {
    stop ("One or more addresses is unnamed.")
  }
  if (length(names(GFFAddress)) != length(GFFAddress)) {
    stop ("Length of names does not match length of object.")
  }
  
  GeneCalls <- vector("list",
                      length = length(GFFAddress))
  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  for (i in seq_along(GFFAddress)) {
    z1 <- gzcon(url(GFFAddress[i]))
    z2 <- textConnection(readLines(z1))
    z3 <- readLines(z2)
    z4 <- strsplit(z3,
                   split = "\t")
    TotalIndices <- sapply(z4,
                           function(x) ifelse(test = length(x) == 9L,
                                              yes = x[1],
                                              no = NA))
    TotalIndices <- unique(TotalIndices)[!is.na(unique(TotalIndices))]
    Index <- sapply(z4,
                    function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                       yes = x[1],
                                       no = NA))
    Start <- sapply(z4,
                    function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                       yes = x[4],
                                       no = NA),
                    USE.NAMES = FALSE)
    Stop <- sapply(z4,
                   function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                      yes = x[5],
                                      no = NA),
                   USE.NAMES = FALSE)
    Strand <- sapply(z4,
                     function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                        yes = x[7],
                                        no = NA),
                     USE.NAMES = FALSE)
    Strand <- ifelse(test = Strand == "+",
                     yes = 0L,
                     no = 1L)
    Product <- sapply(z4,
                      function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                         yes = x[9],
                                         no = NA),
                      USE.NAMES = FALSE)
    Product <- str_extract(Product,
                           "(?<=product=)(.*)(?=;protein_id)")
    Index <- Index[which(!is.na(Start) &
                           !is.na(Stop) &
                           !is.na(Strand) &
                           !is.na(Product))]
    Indices <- sapply(Index,
                      function(x) which(TotalIndices == x),
                      USE.NAMES = FALSE,
                      simplify = TRUE)
    
    z5 <- as.integer(Start[which(!is.na(Start) &
                                   !is.na(Stop) &
                                   !is.na(Strand) &
                                   !is.na(Product))])
    z6 <- as.integer(Stop[which(!is.na(Start) &
                                  !is.na(Stop) &
                                  !is.na(Strand) &
                                  !is.na(Product))])
    z7 <- Strand[which(!is.na(Start) &
                         !is.na(Stop) &
                         !is.na(Strand) &
                         !is.na(Product))]
    z8 <- Product[which(!is.na(Start) &
                          !is.na(Stop) &
                          !is.na(Strand) &
                          !is.na(Product))]
    GeneCalls[[i]] <- data.frame("Index" = Indices,
                                 "Start" = z5,
                                 "Stop" = z6,
                                 "Strand" = z7,
                                 "Annotation" = z8,
                                 stringsAsFactors = FALSE)
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = i/length(GFFAddress))
    }
  }
  names(GeneCalls) <- names(GFFAddress)
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    TotalTime <- (TimeStop - TimeStart)
    cat("\n")
    print(TotalTime)
  }
  return(GeneCalls)
}

