\name{}
\alias{}
\usage{
ParseSubGraphs(SubGraphs, DBPATH, GeneCalls, CopyType = "Single",
  InputType = "Groups", OriginGraph, Verbose = TRUE)
}
\arguments{

\item{SubGraphs}{A list of either subgraphs from a large adjacency graph, or character vectors identifying positions on an adjacency graph.}

\item{DBPATH}{A character vector specifying a sqlite database built by DECIPHER.}

\item{GeneCalls}{A named list of dataframes providing columns named Start, Stop, Index, and Strand.}

\item{CopyType}{Single to select on clusters with no competing nodes from the same genome Multi allows for all possible clusters.}

\item{InputType}{Groups indicates Subgraphs are a list of character vectors Graphs indicates a list of subgraphs.}

\item{OriginGraph}{The graph from which the child subraphs were generated.}

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
