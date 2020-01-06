\name{PairSummaries}
\alias{PairSummaries}
\title{Summarize syntenically genes and genomic features.}
\usage{
PairSummaries(SyntenyLinks, GeneCalls, BPATH, PIDs = TRUE, IgnoreDefaultStringSet = FALSE, Verbose = TRUE, GapPenalty = TRUE, TerminalPenalty = TRUE, Model = "Global", Correction = "none")
}
\arguments{

\item{SyntenyLinks}{An object specifying where syntenic hits connect pairs of genes.}

\item{GeneCalls}{A list of data.frames.}

\item{DBPATH}{A PATH to the associated DECIPHER sqlite database.}

\item{PIDS}{Logical specifying whether to align pairs. TRUE by default}

\item{IgnoreDefaultStringSet}{Logical, if TRUE all alignments will be nucleotide alignments.}

\item{Verbose}{Print progress bar?}

\item{GapPenalty}{Argument to be passed to DECIPHER's AlignSeqs.}

\item{TerminalPenalty}{Argument to be passed to DECIPHER's AlignSeqs.}

\item{Model}{A character specifying whether to use a basic model for rejecting links that appear to be spurious. If NULL, all links will be kept. If a user specified model has been created, it can be specified here as an object.}

\item{Correction}{Argument to be passed to DECIPHER's DistanceMatrix.}

}

\description{
Summarize syntenic links.
}

\examples{
PairSummaries()
}

\keyword{Synteny}
\keyword{Orthologs}