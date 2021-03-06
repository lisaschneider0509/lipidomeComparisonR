% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipidome_comparison_pca.R
\name{plot_contrib_to_pc}
\alias{plot_contrib_to_pc}
\title{Plot contribution to pca}
\usage{
plot_contrib_to_pc(
  pca_element,
  title = "Contribution of variables to the principal components",
  out_path = "none"
)
}
\arguments{
\item{pca_element}{Object produced by FactoMineR::PCA.}

\item{title}{string. Main title of the plot. Default = "Contribution of variables to the principal components"}

\item{out_path}{string. Path to save plot as png. If "none" (default), plot is printed to device.}
}
\description{
`plot_contrib_to_pc` plots the contributions of each variable to each principal component in a grid.
}
