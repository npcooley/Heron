######
# TODO
# 1: if genes contain multiple CDSs with the same stop or start
# keep only the longer CDS
# 2: if a CDS is VERY short
# remove it
# 3: add product line from Note=, if no other product line exists, and a note exists
# 4: how to deal with mobile_genetic_elemets and CRISPR arrays?
######


gffToDataFrame <- function(gff,
                           additionalAttrs = NULL,
                           additionalTypes = NULL,
                           rawTableOnly = FALSE) {
  
  if (!("stringr" %in% .packages())) {
    stop ("Required package stringr is not loaded.")
  }
  
  if (grepl(pattern = "www.|http:|https:|ftp:|ftps:",
            x = gff)) {
    CONN <- gzcon(url(gff))
    Z01 <- readLines(CONN)
    close.connection(con = CONN)
  } else {
    Z01 <- readLines(gff)
  }
  
  T01 <- strsplit(Z01,
                  split = "\t",
                  fixed = TRUE)
  T01 <- T01[which(lengths(T01) == 9L)]
  T01 <- do.call(rbind,
                 T01)
  T02 <- strsplit(T01[, 9L],
                  split = ";",
                  fixed = TRUE)
  T01 <- T01[, -9L]
  
  # create a DF with the correct raw formats
  
  T03 <- as.data.frame(T01,
                       stringsAsFactors = FALSE)
  
  T03[, 4L] <- as.integer(T03[, 4L])
  T03[, 5L] <- as.integer(T03[, 5L])
  T03[, 7L] <- ifelse(test = T03[, 7L] == "+",
                      yes = 0L,
                      no = 1L)
  
  CollectAttr <- c("ID",
                   "Parent",
                   "Name",
                   "gbkey",
                   "gene",
                   "product",
                   "protein_id",
                   "gene_biotype")
  
  if (!is.null(additionalAttrs)) {
    CollectAttr <- c(CollectAttr,
                     additionalAttrs)
  }
  
  # rectangularize attribute fields
  
  T04 <- sapply(T02,
                function(y) sapply(CollectAttr,
                                   function(z) if (length(grep(pattern = paste(z,
                                                                             "=",
                                                                             sep = ""),
                                                             x = y)) == 1L) {
                                     str_extract(string = y[grepl(pattern = paste(z,
                                                                                  "=",
                                                                                  sep = ""),
                                                                  x = y)],
                                                 pattern = paste("(?<=",
                                                                 z,
                                                                 "=)",
                                                                 "(.*)",
                                                                 sep = ""))
                                   } else if (length(grep(pattern = paste(z,
                                                                          "=",
                                                                          sep = ""),
                                                          x = y)) == 0L) {
                                     NA
                                   } else if (length(grep(pattern = paste(z,
                                                                          "=",
                                                                          sep = ""),
                                                          x = y)) > 1L) {
                                     str_extract(string = y[grepl(pattern = paste(z,
                                                                                  "=",
                                                                                  sep = ""),
                                                                  x = y)],
                                                 pattern = paste("(?<=",
                                                                 z,
                                                                 "=)",
                                                                 "(.*)",
                                                                 sep = ""))[1L]
                                     warning("An attribute field is repeated, only the first is collected.")
                                   }))
  
  T04 <- t(T04)
  T04 <- as.data.frame(T04,
                       stringsAsFactors = FALSE)
  
  colnames(T03) <- c("Contig",
                     "Source",
                     "Type",
                     "Start",
                     "Stop",
                     "Score",
                     "Strand",
                     "Phase")
  
  CompleteTable <- cbind(T03,
                         T04,
                         stringsAsFactors = FALSE)
  
  Contigs <- CompleteTable[!duplicated(CompleteTable$Contig), ]
  ContigMaxes <- Contigs$Stop[order(Contigs$Stop,
                                    decreasing = TRUE)]
  Contigs <- Contigs[order(Contigs$Stop,
                           decreasing = TRUE), ]
  Contigs <- Contigs$Contig
  
  Index <- sapply(CompleteTable$Contig,
                  function(x) which(Contigs == x),
                  USE.NAMES = FALSE,
                  simplify = TRUE)
  
  CompleteTable <- cbind(CompleteTable,
                         "Index" = Index)
  
  CompleteTable <- CompleteTable[order(CompleteTable$Index,
                                       CompleteTable$Start), ]
  
  if (rawTableOnly) {
    return(CompleteTable)
  }
  
  MatchTypes <- c("gene",
                  "pseudogene")
  
  if (!is.null(additionalTypes)) {
    MatchTypes <- c(MatchTypes,
                    additionalTypes)
  }
  
  MatchTable <- CompleteTable[CompleteTable$Type %in% MatchTypes, ]
  MatchTable <- MatchTable[order(MatchTable$Index,
                                 MatchTable$Start), ]
  
  # Correctly match Coding sequences to genes ... 
  CodingSelect <- vector(mode = "logical",
                         length = nrow(MatchTable))
  MatchLine <- vector(mode = "character",
                      length = length(CodingSelect))
  ProductLine <- vector(mode = "character",
                        length = length(MatchLine))
  for (m1 in seq_len(nrow(MatchTable))) {
    if (!is.na(MatchTable$gene_biotype[m1])) {
      # biotype is not NA ... evaluate further
      CurrentLeft <- MatchTable$Start[m1]
      CurrentRight <- MatchTable$Stop[m1]
      CurrentIndex <- MatchTable$Index[m1]
      CurrentStrand <- MatchTable$Strand[m1]
      CurrentID <- MatchTable$ID[m1]
      SelectLines <- which(CompleteTable$Start >= CurrentLeft &
                             CompleteTable$Stop <= CurrentRight &
                             CompleteTable$Index == CurrentIndex &
                             CompleteTable$Strand == CurrentStrand &
                             CompleteTable$Type == "CDS")
      if (any(!is.na(CompleteTable[SelectLines, "Parent"]))) {
        # attempt to make sure that overlapping genes are not including CDSs from
        # each other
        ParentCheck <- !is.na(CompleteTable[SelectLines, "Parent"]) & CompleteTable[SelectLines, "Parent"] != CurrentID
        if (any(ParentCheck)) {
          SelectLines <- SelectLines[!ParentCheck]
        }
      }
      if (length(SelectLines) > 1L) {
        # if more than one CDSs are selected
        CodingSelect[m1] <- TRUE
        CurrentLines <- CompleteTable[SelectLines, c("Start", "Stop")]
        MatchLine[m1] <- paste(apply(X = CurrentLines,
                                     MARGIN = 1,
                                     FUN = function(x) paste(x[1],
                                                             x[2],
                                                             sep = "X")),
                               collapse = "Y")
        PossibleProducts <- CompleteTable[SelectLines, "product"]
        if (any(!is.na(PossibleProducts))) {
          PossibleProducts <- PossibleProducts[!is.na(PossibleProducts)]
          if (length(unique(PossibleProducts)) == 1L) {
            ProductLine[m1] <- unique(PossibleProducts)
          } else {
            ProductLine[m1] <- PossibleProducts[1L]
            warning("A gene contains CDSs with differing product annotations.")
          }
        }
      } else if (length(SelectLines) == 1L) {
        # if only one CDS is captured
        CodingSelect[m1] <- TRUE
        MatchLine[m1] <- paste(CompleteTable[SelectLines, "Start"],
                               CompleteTable[SelectLines, "Stop"],
                               sep = "X")
        if (!is.na(CompleteTable[SelectLines, "product"])) {
          ProductLine[m1] <- CompleteTable[SelectLines, "product"]
        }
      } else if (length(SelectLines) < 1L) {
        # if no CDSs are captured
        CodingSelect[m1] <- FALSE
        MatchLine[m1] <- paste(CurrentLeft,
                               CurrentRight,
                               sep = "X")
        # if there are any associated lines that have a product field, include
        CheckLines <- which(CompleteTable$Start >= CurrentLeft &
                              CompleteTable$Stop <= CurrentRight &
                              CompleteTable$Index == CurrentIndex &
                              CompleteTable$Strand == CurrentStrand)
        LineSubSet <- which(!is.na(CompleteTable[CheckLines ,"product"]))
        if (length(LineSubSet) >= 1L) {
          PotentialProducts <- unique(CompleteTable[CheckLines[LineSubSet], "product"])
          if (length(PotentialProducts) == 1L) {
            ProductLine[m1] <- PotentialProducts
          } else {
            ProductLine[m1] <- PotentialProducts[1L]
            warning("A gene without CDSs contains differing annotations.")
          }
        }
      }
    } else {
      # biotype is NA just go with line start and stops
      CodingSelect[m1] <- FALSE
      MatchLine[m1] <- paste(MatchTable[m1, "Start"],
                             MatchTable[m1, "Stop"],
                             sep = "X")
      if (!is.na(MatchTable[m1, "product"])) {
        ProductLine[m1] <- MatchTable[m1, "product"]
      }
    }
  }
  MatchTable <- cbind(MatchTable,
                      "Match" = MatchLine,
                      "Coding" = CodingSelect,
                      "Product" = ProductLine,
                      stringsAsFactors = FALSE)
  MatchTable <- MatchTable[, c("Index",
                               "Strand",
                               "Start",
                               "Stop",
                               "Type",
                               "gene_biotype",
                               "Match",
                               "Product",
                               "Coding",
                               "Contig")]
  rownames(MatchTable) <- NULL
  
  # rewrite any matches where genes extend over the end of an index
  
  for (m1 in seq_len(nrow(MatchTable))) {
    if (MatchTable$Stop[m1] > ContigMaxes[MatchTable$Index[m1]]) {
      CurrentBounds <- sapply(strsplit(strsplit(MatchTable$Match[m1],
                                                split = "Y",
                                                fixed = TRUE)[[1]],
                                       split = "X",
                                       fixed = TRUE),
                              function(x) as.integer(x))
      BreakPoint <- which(CurrentBounds[1, ] <= ContigMaxes[MatchTable$Index[m1]] &
                            CurrentBounds[2, ] > ContigMaxes[MatchTable$Index[m1]])
      Section <- CurrentBounds[, BreakPoint]
      NewSection <- matrix(data = c(Section[1],
                                    ContigMaxes[MatchTable$Index[m1]],
                                    1L,
                                    Section[2] - ContigMaxes[MatchTable$Index[m1]]),
                           ncol = 2L)
      if (BreakPoint < ncol(CurrentBounds) &
          BreakPoint != 1L) {
        # if breakpoint is not either the first or last bound set
        AdjustedBounds <- CurrentBounds[,
                                        (BreakPoint + 1L):ncol(CurrentBounds),
                                        drop = FALSE] - ContigMaxes[MatchTable$Index[m1]]
        NewBoundSet <- cbind(CurrentBounds[, 1L:(BreakPoint - 1L)],
                             NewSection,
                             AdjustedBounds)
      } else if (BreakPoint == 1L &
                 ncol(CurrentBounds) > 1L) {
        # if breakboint is the first bound set
        CurrentBounds <- CurrentBounds - ContigMaxes[MatchTable$Index[m1]]
        NewBoundSet <- cbind(NewSection,
                             CurrentBounds[, (BreakPoint + 1L):ncol(CurrentBounds)])
      } else if (BreakPoint == ncol(CurrentBounds) &
                 ncol(CurrentBounds) > 1L) {
        # if breakpoint is the last bound set
        NewBoundSet <- cbind(CurrentBounds[, 1:(BreakPoint - 1L)],
                             NewSection)
      } else {
        # breakpoint is the only bound set
        NewBoundSet <- NewSection
      }
      NewMatchSet <- paste(apply(X = NewBoundSet,
                                 MARGIN = 2L,
                                 FUN = function(x) paste(x,
                                                         collapse = "X")),
                           collapse = "Y")
      MatchTable$Match[m1] <- NewMatchSet
    }
  }
  
  return(MatchTable)
}






