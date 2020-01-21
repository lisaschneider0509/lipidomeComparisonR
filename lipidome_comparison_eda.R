### Install packages
# install.packages("ggplot2")
# install.packages("stringr")
install.packages("tidyr")
install.packages("data.table")
install.packages("textshape")
install.packages("tibble")

### load packages
library(gridExtra)
library(stringr)
library(ggplot2)
library(tidyr)
library(data.table)
library(textshape)
library(tibble)

theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

# set variables
input_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/prefix_data.csv"

# load & transform data
lipid_data <- read.csv(input_path, sep = ",", dec = ".", header = TRUE) #read data
t_lipid_data <- transpose(as.data.table(lipid_data), make.names = 1)
t_lipid_data <- as.data.frame(t_lipid_data)
row.names(t_lipid_data) <- colnames(lipid_data[-1])


test_data <- read.csv(test_path, sep = ",", dec = ".", header = TRUE) #read data
t_test_data <- transpose(as.data.table(test_data), make.names = 1)
t_test_data <- as.data.frame(t_test_data)
row.names(t_test_data) <- colnames(test_data[-1])

working_data <- t_test_data



### summary biological & technical replicates
sample_ID <- list(substr(row.names(working_data), 1, 3)) # letter 1-3 identify the sample
# biol_replicate <- list(substr(row.names(working_data), 1, 5)) # letter 1-5 identify the biological replicates

meta_info <- read.table(text = row.names(working_data), sep = "_")
meta_info$V2 <- paste(meta_info$V1, meta_info$V2, sep = "_")
meta_info$V3 <- NULL
working_data <- tibble::add_column(working_data, treatment = meta_info$V1, .before = 1)
working_data <- tibble::add_column(working_data, biol_replicate = meta_info$V2, .after = 1)

means_biol <- as.data.frame(aggregate(working_data[-(1:2)], by=list(working_data$treatment), FUN=mean)) # grouping by sample-ID calculates the mean over all biological replicates
means_tech <- as.data.frame(aggregate(working_data[-(1:2)], by=biol_replicate, FUN=mean)) # grouping by biological replicate calculates the mean over all technical replicates

plot_qq <- function(sample_df,plot_type){
  print(plot_type)
  plotlist <- list()
  for (i in 2:ncol(sample_df[,1: ncol(sample_df)])){
    if (plot_type=="qqplot"){
      p1 <- ggplot(sample_df, aes(sample=sample_df[, i]))
      p1 <- p1 + labs(title = colnames(means_tech)[i])
      p1 <- p1 +stat_qq()+stat_qq_line()
    }
    else if(plot_type=="boxplot"){
          p1 <- ggplot(means_tech, aes(x=means_tech[[1]], y=means_tech[,i])) 
    p1 <- p1 + geom_boxplot()+labs(title = colnames(means_tech)[i])
    }
    else if(plot_type=="histogram"){
      p1 <- ggplot(means_tech, aes(means_tech[,i])) 
      p1 <- p1 + geom_histogram()+labs(title = colnames(means_tech)[i])
    }
    plotlist[[i]] <- p1
  }
  
  n_plots <- length(plotlist)
  nCol <- floor(sqrt(n_plots))
  grid.arrange(grobs = plotlist[2:ncol(sample_df)], ncol = nCol)
}
plot_qq(working_data, "boxplot")




