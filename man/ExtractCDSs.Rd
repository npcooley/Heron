\name{ExtractCDSs}
\alias{ExtractCDSs}
\title{Extract coding regions of a gene.}
\usage{
ExtractCDSs(Match, Start, Stop, Strand = 0L, SourceString, Mode = "Character")
}
\arguments{

\item{Match}{A character string of the format ..X..Y..X.. where X indicates a gap between starts and stops and Y indicates a gap between Starts and Starts.}

\item{Start}{Gene Start if breaks are given for the coding sequence.}

\item{Stop}{Gene Stop if no breaks are given for the coding sequence.}

\item{Strand}{Strand of the coding sequence.}

\item{SourceString}{DNAString from which to pull the gene.}

\item{Mode}{Argument for future use to allow multi-exon genes to be extracted in other formats.}

}

\description{
Allow extraction of individual genes into a DNAStringSet.
}

\examples{
ExtractCDSs()
}

\keyword{Gene}
\keyword{Extract}
\keyword{GeneCall}