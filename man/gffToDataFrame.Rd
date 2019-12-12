\name{gffToDataFrame}
\alias{gffToDataFrame}
\title{Convert gff files directly to a dataframe}
\usage{
gffToDataFrame(gff, additionalAttrs = NULL, additionalTypes = NULL, rawTableOnly = FALSE)
}

\arguments{

\item{gff}{A location or address for a gff3 file.}

\item{additionalAttrs}{Additional attribute fields to include in pull.}

\item{additionalTypes}{Additional line types to parse further.}

\item{rawTableOnly}{Logical specifying whether to parse GFF file into more succinct format.}

}

\description{
Convert a gff file into a raw dataframe, or one that is parsed down to specific entry types.
}

\examples{
gffToDataFrame()
}

\keyword{GeneCalls}
\keyword{gff}