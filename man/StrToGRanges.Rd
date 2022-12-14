% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StrToGRanges.r
\name{StrToGRanges}
\alias{StrToGRanges}
\title{Convert String in GRanges object}
\usage{
StrToGRanges(x.chr_vec)
}
\arguments{
\item{x.chr_vec}{<character>: strings to convert on GRanges.}
}
\value{
a GRanges object.
}
\description{
StrToGRanges
Convert ranges describe with string (i.e seqname:start-end:strand) in GRanges object.
}
\examples{
StrToGRanges("chr1:1-100:+")
StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))
}
