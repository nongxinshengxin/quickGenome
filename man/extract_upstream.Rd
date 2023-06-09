% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_upstream.R
\name{extract_upstream}
\alias{extract_upstream}
\title{Extraction of upstream sequences of genes.}
\usage{
extract_upstream(fasta, gff, outputfile, range = 2000)
}
\arguments{
\item{fasta}{A character vector (of arbitrary length when reading, of length 1 when writing) containing the path(s) to the fasta file(s) to read. Reading files in gzip format (which usually have the '.gz' extension) is supported.

Note that special values like "" or "|cmd" (typically supported by other I/O functions in R) are not supported here.

Note the format of file must be "fasta" or fasta file in gzip.}

\item{gff}{A character vector (of arbitrary length when reading, of length 1 when writing) containing the path(s) to the gff file(s) to read.

Note the format of file must be "gff".}

\item{outputfile}{A character vector (of arbitrary length when reading, of length 1 when writing) containing the path(s) to the fasta file(s) to write.}

\item{range}{The length of upstream sequences. Default is 2000}
}
\description{
Extraction of upstream sequences of genes.
}
\references{
\enumerate{
\item Biostrings
\url{https://doi.org/doi:10.18129/B9.bioc.Biostrings}
}

Please, when using quickGenome, cite us using the reference:
\url{https://github.com/nongxinshengxin/quickGenome}
}
\author{
Cy Gao
}
