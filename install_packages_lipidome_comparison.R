### Install packages
#' Install packages for lipidome comparison workflow
#' @description `install_packages_lipidome_comparison` installs all packages for the lipidome comparison workflow
#' @example 
#' \dontrun
#' install_packages_lipidome_comparison()
install_packages_lipidome_comparison <- function(){
  ## general
  install.packages("dplyr") # select part of data
  install.packages("stringr") # string manipulation
  install.packages("data.table") # better transpose function than r-base
  install.packages("tibble") # data frame manipulation
  ## data transforations
  install.packages("imputeLCMD") # various imputation procedures, including left censored imputation
  install.packages("impute") # dependecy for imputeLCMD
  ## graphics
  install.packages("ggplot2") # plots
  install.packages("GGally") # paralell plot
  install.packages("fmsb") # spider chart
  install.packages("scales") # color scales, for example opacity with alpha
  install.packages("viridis") # color blind and printer save color schemes
  install.packages("ggrepel") # avoids overlapping lables in ggplot
  install.packages("plotly") # interactive ggplots
  install.packages("htmlwidgets") # save plotly-plots as html
  install.packages("ggpubr") # arrange ggplots in grid
  ## PCA
  install.packages("FactoMineR")
  install.packages("ggfortify") # improves autoplot-function in ggplot, for example biplot with ggplot
  install.packages("factoextra") # various graphs for PCA
  install.packages("corrplot")
  ## clustering
  install.packages("dendextend") # find the best clustering parameters
  install.packages("gplots") # nicer heatmaps than ggplot
  install.packages("heatmaply") # interactive hatmaps
  ## Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  BiocManager::install("lipidr")
  
  # install.packages("LipidMS")
  # install.packages("LipidMSdata")
}

