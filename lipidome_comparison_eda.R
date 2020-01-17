# Install packages
# install.packages("ggplot2")

# load packages
library(gridExtra)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "top")
)

# set variables
input_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"
test_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/test.csv"

# load data
lipid_data <- read.csv(input_path, header = TRUE, sep = ",", dec = ".", row.names = 1) #read data
lipid_data <- t(lipid_data)

test_data <- read.csv(test_path, header = TRUE, sep = ",", dec = ".", row.names = 1) #read data
test_data <- as.data.frame(t(test_data))

working_data <- lipid_data

# summary biological & technical replicates
sample_ID <- list(substr(row.names(working_data), 1, 3)) # letter 1-3 identify the sample
biol_replicate <- list(substr(row.names(working_data), 1, 5)) # letter 1-5 identify the biological replicates

means_biol <- as.data.frame(aggregate(working_data, by=sample_ID, FUN=mean)) # grouping by sample-ID calculates the mean over all biological replicates
means_tech <- as.data.frame(aggregate(working_data, by=biol_replicate, FUN=mean)) # grouping by biological replicate calculates the mean over all technical replicates

plotlist <- list()
for (i in 2:ncol(means_tech[,1: ncol(means_tech)])){  
  p1 <- ggplot(means_tech, aes(sample=means_tech[, i]))+stat_qq()+stat_qq_line()+labs(title = colnames(means_tech)[i])
  plotlist[[i]] <- p1
}

n_plots <- length(plotlist)
nCol <- floor(sqrt(n_plots))
grid.arrange(grobs = plotlist[2:ncol(means_tech)], ncol = nCol)

