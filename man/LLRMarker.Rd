% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markerREC.R
\name{LLRMarker}
\alias{LLRMarker}
\title{Calculate marker specificity}
\usage{
LLRMarker(mat, clusterIDs, binmethod = "min", expR = 0.1, cpu = 1, marker = 20)
}
\arguments{
\item{mat}{expression matrix}

\item{clusterIDs}{cluster IDs for each cell.}

\item{binmethod}{Method used to binarize the matrix}

\item{expR}{The proportion of feature expressed cells in a cluster.}

\item{cpu}{How many cores will be used for detection of markers}

\item{marker}{The number of markers returned.}
}
\value{
data.frame with marker features of each cluster.
}
\description{
Calculate marker specificity
}
\author{
Qiong Zhang
}
