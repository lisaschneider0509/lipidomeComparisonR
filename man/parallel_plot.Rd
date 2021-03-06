% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipidome_comparison_visualization.R
\name{parallel_plot}
\alias{parallel_plot}
\title{Parallel coordinates plot}
\usage{
parallel_plot(
  input_df,
  group_vector,
  out_path = "none",
  title = "Parallel Plot",
  x_title = "",
  y_title = "",
  legend_title = "",
  col_names = colnames(dplyr::select_if(input_df, is.numeric)),
  scale = "globalminmax"
)
}
\arguments{
\item{input_df}{data frame.}

\item{group_vector}{vector. Vector (part of the data frame), to sort by.}

\item{out_path}{string. Path to save parallel plot to png. If out_path is empty the parallel plot is printed to the device.}

\item{title}{string. Main Title. Default = "Parallel plot"}

\item{x_title}{string. Title of x-axis. Default: none}

\item{y_title}{string. Title of y-axis. Default: none}

\item{legend_title}{string. Title of legend. Default: none}

\item{col_names}{vector. Column names of the data frame and lables on the x-axis. Default: colnames(input_df)}

\item{scale}{string. Method used to scale the variable. Default = "globalminmax". Options: "std", "robust", "uniminmax", "globalminmax", "center", "centerObs". For more information on the options see help(ggparcoord).}
}
\description{
`parallel_plot` prints a paralell coordinates plot using a data frame
}
\details{
This function takes a data frame with at least one factor variable and displays it in a pralell coordinates plot, where the different groups are color coded. This plot works best for <= 10 parameters.
}
\examples{
parallel_plot(iris, iris$Species)
\dontrun{
dir.create(paste(getwd(), "/examples", sep = ""), showWarnings = FALSE)
dir <- paste(getwd(), "/examples/iris", sep = "")
parallel_plot(iris, 5, out_path = dir)}
}
