% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipidome_comparison_visualization.R
\name{spider_chart}
\alias{spider_chart}
\title{Spider Chart}
\usage{
spider_chart(
  minimized_df,
  title = "Spider chart",
  legend_lab = rownames(minimized_df),
  out_path = "none"
)
}
\arguments{
\item{minimized_df}{data frame. A data frame that contains only one value per group and variable. Most often this will be a data frame of means calculated from another data frame.}

\item{title}{string. Main title of the chart. Default = "Spider chart"}

\item{legend_lab}{vector. Lables. Default: rownames(minimized_df)}

\item{out_path}{string. Path to save spider chart to png. If out_path is empty, the spider chart is printed to the device.}
}
\description{
`spider chart` takes a minimized data frame (one value per group) and prints a spider chart
}
\details{
{This function takes a data frame of with one value per group (i.e. calculate mean groupwise).
This minimized data frame is used to draw a spider chart (also radar chart odr network plot).
The ideal numer of parameters for a spider chart is <= 10. Also the shape of the graph depends on
the order of parameters. If the data has large differences in size, normalizing or scaling the data is necessary.
Doesn't work with non normal data.}
}
\examples{
minimized_iris <- aggregate(dplyr::select_if(iris, is.numeric),
                            by = list(iris$Species),
                            FUN = mean)
rownames(minimized_iris) <- minimized_iris$Group.1
spider_chart(minimized_iris, title = "Spider chart of iris", legend_lab = c(1, 2, 3))
}
