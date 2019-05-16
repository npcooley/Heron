\name{ReciprocalBestBlast}
\alias{ReciprocalBestBlast}
\title{A generic reciprocal best blast search.}
\usage{
ReciprocalBestBlast(PATH, ListOfGeneCalls, GenomeIDs, Align.Fast = FALSE)
}
\arguments{
\item{PATH}{A character vector specifying the location of a sqlite database built by DECIPHER.}

\item{ListOfGeneCalls}{A named list of data.frames, one for each genome in SyntenyObject, with named columns "Index", "Start", "Stop", "Strand", and an optional "Annotation" column.}

\item{GenomeIDs}{A vector of characters specifying the identifiers to select from the database specified by PATH, and the list locations to select from in ListOfGeneCalls.}

\item{Align.Fast}{Logical, run semi-redundant DECIPHER alignment of paired sequences.}

}
\description{
Perform reciprocal best BLAST searches on set of genomes and their associated genomes specified by GenomeIDs. Uses an e-value cutoff of 1e-6, and the smith-waterman traceback options. BLAST must be installed locally and accessible by R.
}
\examples{
ReciprocalBestBlast()
}
\keyword{BLAST}
\keyword{Alignment}
