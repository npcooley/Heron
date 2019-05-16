\name{HSSPScore}
\alias{HSSPScore}
\title{A function for determining the HSSP score for a given alignment}
\usage{
HSSPScore(PID, MatchLength)
}
\arguments{

\item{PID}{Percent Identity calculated by BLAST.}

\item{MatchLength}{Number of Matches between two aligned sequences}
}
\description{Determine HSSP Score, a modification of PID that takes match length into account.}
\examples{
HSSPScore()
}
\keyword{Alignment}
\keyword{Identity}