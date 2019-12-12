\name{ParseSubGraphs}
\alias{ParseSubGraphs}
\title{A function for parsing subgraphs of a large adjacency graph.}
\usage{
ParseSubGraphs(Groups, DBPATH, GeneCalls, IgnoreDefault = FALSE, Verbose = TRUE)
}
\arguments{

\item{Groups}{A list of either subgraphs from a large adjacency graph, or character vectors identifying positions on an adjacency graph.}

\item{DBPATH}{A character vector specifying a sqlite database built by DECIPHER.}

\item{GeneCalls}{A named list of dataframes providing columns named Start, Stop, Index, and Strand.}

\item{IgnoreDefault}{Logical, if TRUE all alignments will be DNA Alignments.}

\item{Verbose}{Print progress bars to screen}

}

\description{
Collect alignments, distance matrices, dendrograms, and phylograms for subgraphs of genes
}
\examples{
ParseSubGraphs()
}
\keyword{Graph}
\keyword{SubGraph}
\keyword{Alignment}
\keyword{Orthology}
