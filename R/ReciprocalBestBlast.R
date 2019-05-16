######
# Perform reciprocal best blast as needed given:
# a series of genomes
# a DBPath
# and a series of predictions for those genes
######

ReciprocalBestBlast <- function(PATH,
                                ListOfGeneCalls,
                                GenomeIDs,
                                Align.Fast = FALSE) {
  TotalTimeStart <- Sys.time()
  stopifnot("DECIPHER" %in% .packages())
  stopifnot(length(GenomeIDs) == length(ListOfGeneCalls))
  source("20180501_HSSP_v01.R",
         echo = FALSE)
  ######
  # Initialize Results
  ######
  L <- length(GenomeIDs)
  ResultsMatrix <- matrix(data = vector("list",
                                        length = 1L),
                          nrow = L,
                          ncol = L)
  pBar <- txtProgressBar(style = 1L)
  ######
  # Collect current genomes
  ######
  for (m1 in 1L:(L - 1L)) {
    for (m2 in (m1 + 1L):L) {
      Q.Genome <- SearchDB(dbFile = PATH,
                           tblName = "Seqs",
                           identifier = as.character(GenomeIDs[m1]),
                           type = "XStringSet",
                           nameBy = "identifier",
                           verbose = FALSE)
      S.Genome <- SearchDB(dbFile = PATH,
                           tblName = "Seqs",
                           identifier = as.character(GenomeIDs[m2]),
                           type = "XStringSet",
                           nameBy = "identifier",
                           verbose = FALSE)
      ######
      # Collect all genes
      ######
      cat("\n")
      cat("Getting Query Genes... ")
      cat("\n")
      Q.Genes <- unlist(extractAt(x = DNAStringSet(sapply(ListOfGeneCalls[[m1]][, "Index"],
                                                          function(x) Q.Genome[[x]])),
                                  at = IRangesList(mapply(function(z, y) IRanges(start = z,
                                                                                 end = y),
                                                          z = ListOfGeneCalls[[m1]][, "Start"],
                                                          y = ListOfGeneCalls[[m1]][, "Stop"]))))
      # Q.Genes <- extractAt(x = Q.Genome[[1]],
      #                      at = IRanges(start = ListOfGeneCalls[[m1]][, "Start"],
      #                                   end = ListOfGeneCalls[[m1]][, "Stop"]))
      RC <- which(ListOfGeneCalls[[m1]][, "Strand"] == "1")
      Q.Genes[RC] <- reverseComplement(Q.Genes[RC])
      names(Q.Genes) <- as.character(1:length(Q.Genes))
      cat("...Getting Subject Genes")
      cat("\n")
      S.Genes <- unlist(extractAt(x = DNAStringSet(sapply(ListOfGeneCalls[[m2]][, "Index"],
                                                          function(x) S.Genome[[x]])),
                                  at = IRangesList(mapply(function(z, y) IRanges(start = z,
                                                                                 end = y),
                                                          z = ListOfGeneCalls[[m2]][, "Start"],
                                                          y = ListOfGeneCalls[[m2]][, "Stop"]))))
      # S.Genes <- extractAt(x = S.Genome[[1]],
      #                      at = IRanges(start = ListOfGeneCalls[[m2]][, "Start"],
      #                                   end = ListOfGeneCalls[[m2]][, "Stop"]))
      RC <- which(ListOfGeneCalls[[m2]][, "Strand"] == "1")
      S.Genes[RC] <- reverseComplement(S.Genes[RC])
      names(S.Genes) <- as.character(1:length(S.Genes))
      Q.Prot <- translate(Q.Genes,
                          if.fuzzy.codon = "solve")
      S.Prot <- translate(S.Genes,
                          if.fuzzy.codon = "solve")
      ######
      # Make BLAST DB
      ######
      cat("\n")
      cat("Performing Forward BLAST... ")
      dir.create(paste0("results"))
      writeXStringSet(S.Prot,
                      filepath = "temp1.fa")
      writeXStringSet(Q.Prot,
                      filepath = "temp2.fa")
      system("makeblastdb -in temp1.fa -dbtype prot -out temp1")
      system(paste0("blastp -query temp2.fa -outfmt 6 -evalue 1e-6 -use_sw_tback -db temp1 -out ./results",
                    "/",
                    "Forward",
                    ".txt"))
      ######
      # Read in results as DF, discard .txt output from BLAST
      ######
      ForwardResult <- read.table(paste0("./results",
                                         "/",
                                         "Forward",
                                         ".txt"),
                                  sep="\t",
                                  header=F)
      unlink(paste0("./results",
                    "/",
                    "Forward",
                    ".txt"))
      ######
      # Reciprocal BLAST with previous query as DB
      ######
      cat("\n")
      cat("...Performing reciprocal BLAST")
      system("makeblastdb -in temp2.fa -dbtype prot -out temp2")
      system(paste0("blastp -query temp1.fa -outfmt 6 -evalue 1e-6 -use_sw_tback -db temp2 -out ./results",
                    "/",
                    "Reciprocal",
                    ".txt"))
      ######
      # Read in result as DF and discard blast .txt output
      ######
      ReciprocalResult <- read.table(paste0("./results",
                                            "/",
                                            "Reciprocal",
                                            ".txt"),
                                     sep="\t",
                                     header=F)
      unlink(paste0("./results",
                    "/",
                    "Reciprocal",
                    ".txt"))
      ######
      # Remove extraneous files and things
      ######
      system("rm temp1* temp2*")
      system("rm -rf ./results")
      cat("\n")
      cat("Reciprocal BLASTS Completed... ")
      ######
      # Find Best Match from filtered result
      ######
      begin <- 1L
      u <- unique(ForwardResult$V1)
      BestForwardHit <- vector("list",
                               length = length(u))
      for (k in seq_along(u)) {
        w <- .Call("multiMatch",
                   ForwardResult$V1,
                   u[k],
                   begin)
        BestForwardHit[[k]] <- w[which.max(ForwardResult$V12[w])]
        begin <- w[length(w)] + 1L
      }
      ######
      # Filter ForwardResult further to solely best hits for each sequence
      ######
      CondensedForward <- ForwardResult[unlist(BestForwardHit), ]
      ######
      # Find Best Match from filtered result
      ######
      begin <- 1L
      u <- unique(ReciprocalResult$V1)
      BestReciprocalHit <- vector("list",
                                  length = length(u))
      for (k in seq_along(u)) {
        w <- .Call("multiMatch",
                   ReciprocalResult$V1,
                   u[k],
                   begin)
        BestReciprocalHit[[k]] <- w[which.max(ReciprocalResult$V12[w])]
        begin <- w[length(w)] + 1L
      }
      ######
      # Filter ReciprocalResult further to solely best hits for each sequence
      ######
      CondensedReciprocal <- ReciprocalResult[unlist(BestReciprocalHit), ]
      ######
      # Which indices match to the reciprocal indices in the DF for the opposing direction
      ######
      cat("Filtered... ")
      y <- vector("list",
                  length = nrow(CondensedReciprocal))
      
      z <- vector("list",
                  length = nrow(CondensedForward))
      for(i in seq_along(y)) {
        y[[i]] <- c(CondensedReciprocal$V1[i], CondensedReciprocal$V2[i])
      }
      for(i in seq_along(z)) {
        z[[i]] <- c(CondensedForward$V2[i], CondensedForward$V1[i])
      }
      y1 <- y %in% z
      z1 <- z %in% y
      ######
      # Initialize vectors for HSSP & EVal, for both directions
      # Fill in the HSSP and EVals that change for hits that are
      # reciprocal pairs
      ######
      ForHSSP <- rep(NA_real_,
                     length(Q.Genes))
      RecHSSP <- rep(NA_real_,
                    length(S.Genes))
      ForEVal <- rep(NA_real_,
                     length(Q.Genes))
      RecEVal <- rep(NA_real_,
                     length(S.Genes))
      for (i in seq_along(CondensedForward[z1, "V1"])) {
        ForHSSP[CondensedForward[z1, "V1"]][i] <- HSSPScore(PID = CondensedForward[z1, "V3"][i],
                                                            MatchLength = CondensedForward[z1, "V4"][i])
        ForEVal[CondensedForward[z1, "V1"]][i] <- CondensedForward[z1, "V11"][i]
      }
      for (i in seq_along(CondensedReciprocal[y1, "V1"])) {
        RecHSSP[CondensedReciprocal[y1, "V1"]][i] <- HSSPScore(PID = CondensedReciprocal[y1, "V3"][i],
                                                               MatchLength = CondensedReciprocal[y1, "V4"][i])
        RecEVal[CondensedReciprocal[y1, "V1"]][i] <- CondensedReciprocal[y1, "V11"][i]
      }
      ForHSSP <- ForHSSP[!is.na(ForHSSP)]
      ForEVal <- ForEVal[!is.na(ForEVal)]
      RecHSSP <- RecHSSP[!is.na(RecHSSP)]
      RecEVal <- RecEVal[!is.na(RecEVal)]
      HSSP <- mapply(max,
                     ForHSSP,
                     RecHSSP)
      EVal <- mapply(min,
                     ForEVal,
                     RecEVal)
      cat("HSSP and EVals Calculated... ")
      ######
      # Grab sequences that are paired
      # Get DECIPHER alignments for those pairs
      ######
      ReciprocallyPairedSeqs <- vector("list",
                                       length = nrow(CondensedForward[z1, ]))
      Q.IDs <- CondensedForward[z1, "V1"]
      S.IDs <- CondensedForward[z1, "V2"]
      cat("Collecting Best Blast Pairs")
      cat("\n")
      if (Align.Fast) {
        for (k in seq_along(ReciprocallyPairedSeqs)) {
          ReciprocallyPairedSeqs[[k]] <- c(DECIPHER:::.subset(x = Q.Prot,
                                                              i = Q.IDs[k]),
                                           DECIPHER:::.subset(x = S.Prot,
                                                              i = S.IDs[k]))
          setTxtProgressBar(pb = pBar,
                            value = k/length(ReciprocallyPairedSeqs))
        }
        cat("\n")
        cat("Aligning Best Blast Pairs")
        cat("\n")
        Alignments <- vector("list",
                             length = length(ReciprocallyPairedSeqs))
        for (k in seq_along(Alignments)) {
          Alignments[[k]] <- AlignProfiles(pattern = DECIPHER:::.subset(x = ReciprocallyPairedSeqs[[k]],
                                                                        i = 1L),
                                           subject = DECIPHER:::.subset(x = ReciprocallyPairedSeqs[[k]],
                                                                        i = 2L))
        }
        
        DistMatrices <- sapply(Alignments,
                               function(x) DistanceMatrix(x,
                                                          includeTerminalGaps = TRUE,
                                                          verbose = FALSE),
                               simplify = FALSE)
        PIDs <- sapply(DistMatrices,
                       function(x) (1 - x[3]),
                       simplify = TRUE)
      } else {
        for (k in seq_along(ReciprocallyPairedSeqs)) {
          ReciprocallyPairedSeqs[[k]] <- c(DECIPHER:::.subset(x = Q.Genes,
                                                              i = Q.IDs[k]),
                                           DECIPHER:::.subset(x = S.Genes,
                                                              i = S.IDs[k]))
          setTxtProgressBar(pb = pBar,
                            value = k/length(ReciprocallyPairedSeqs))
        }
        cat("\n")
        cat("Aligning Best Blast Pairs")
        cat("\n")
        Alignments <- sapply(ReciprocallyPairedSeqs,
                             function(x) AlignTranslation(x,
                                                          type = "AAStringSet",
                                                          readingFrame = 1L,
                                                          verbose = FALSE),
                             simplify = FALSE)
        
        DistMatrices <- sapply(Alignments,
                               function(x) DistanceMatrix(x,
                                                          includeTerminalGaps = TRUE,
                                                          verbose = FALSE),
                               simplify = FALSE)
        PIDs <- sapply(DistMatrices,
                       function(x) (1 - x[3]),
                       simplify = TRUE)
      }
      ###### 
      # Return both sets of scores that matter,
      # Return both condensed DFs,
      # if columns 1 and 2 are not identical, there is a problem
      # other columns can and will differ
      ######
      colnames(CondensedForward) <- c("qseqid",
                                      "sseqid",
                                      "pident",
                                      "length",
                                      "mismatch",
                                      "gapopen",
                                      "qstart",
                                      "qend",
                                      "sstart",
                                      "send",
                                      "evalue",
                                      "bitscore")
      colnames(CondensedReciprocal) <- c("qseqid",
                                         "sseqid",
                                         "pident",
                                         "length",
                                         "mismatch",
                                         "gapopen",
                                         "qstart",
                                         "qend",
                                         "sstart",
                                         "send",
                                         "evalue",
                                         "bitscore")
      ResultsMatrix[m1, m2] <- list(list(Q.Genes,
                                         S.Genes,
                                         Q.Prot,
                                         S.Prot,
                                         HSSP,
                                         EVal,
                                         ReciprocallyPairedSeqs,
                                         Alignments,
                                         DistMatrices,
                                         PIDs,
                                         CondensedForward[z1, ],
                                         CondensedReciprocal[y1, ]))
      cat("Results Matrix Updated")
    }
  }
  TotalTimeEnd <- Sys.time()
  print(TotalTimeEnd - TotalTimeStart)
  return(ResultsMatrix)
}



