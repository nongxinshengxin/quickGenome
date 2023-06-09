% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_exonNum.R
\name{calculate_exonNum}
\alias{calculate_exonNum}
\title{The number of gene exons was counted according to gff/gtf.}
\usage{
calculate_exonNum(filepath, format = "gff", collapse = FALSE)
}
\arguments{
\item{filepath}{Input GFF3 or GTF file. Can be a path to a file, or an URL, or a connection object, or a GFF3File or GTFFile object.}
\item{format}{Format of the input file. Accepted values are: "gff" (the default), or "gtf".}
\item{collapse}{TRUE or FALSE (the default). If TRUE, genes with more than 10 exons will be collapsed together and counted.}
}
\value{
Eventually, return a data frame describing the number of genes with different numbers of exons.
}
\description{
The number of gene exons was counted according to gff/gtf.
}
\references{
\enumerate{
\item GenomicFeatures
\url{https://doi.org/doi:10.18129/B9.bioc.GenomicFeatures}
}

Please, when using quickGenome, cite us using the reference:
\url{https://github.com/nongxinshengxin/quickGenome}
}
\author{
Hm Wang
}
