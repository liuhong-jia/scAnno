% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random.R
\name{get_random_score}
\alias{get_random_score}
\title{get_random_score}
\usage{
get_random_score(
  obj.seu,
  ref.obj,
  ref.markers,
  sig.mat,
  max.names,
  second.names
)
}
\arguments{
\item{obj.seu}{Seurat object, which need to be annotated.}

\item{ref.obj}{Reference Seurat object.}
}
\value{
Cell Type Score from Permutation Test.
}
\description{
get_random_score
}
