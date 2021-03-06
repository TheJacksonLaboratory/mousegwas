% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gemma_tools.R
\name{average_strain}
\alias{average_strain}
\title{Return the beta of each strain. Strains are similar if identical in SNPs}
\usage{
average_strain(
  strains_genomes,
  phenotypes,
  covars,
  downsample,
  sex,
  strain = NULL
)
}
\arguments{
\item{strains_genomes}{A table of mice genotypes}

\item{phenotypes}{A table of mice phenotypes}

\item{covars}{Covariates to include in lm}

\item{downsample}{Maximal number of representatives. If 0 use average}

\item{sex}{A vector which assign sex to each individual, to be crossed with strain}

\item{strain}{A vector with strain name. If given use it instead of matching.}
}
\value{
A list with $genotypes and $phenotypes
}
\description{
Return the beta of each strain. Strains are similar if identical in SNPs
}
