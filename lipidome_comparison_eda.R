# Install packages

# load packages

# set variables
data_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/Probe-Datensatz_lisa.csv"

# load data
input_data <- read.csv(data_path, header = TRUE, sep = ",", dec = ".", row.names = 1) #read data
transposed_data <- t(input_data)

# summary biol. replicates
sample_id <- row.names(transposed_data)
lipids <- colnames(transposed_data)

x <- sample_id[1]


for(j in 1:nrow(transposed_data)){
  if(substr(sample_id[j], 1, 3)==substr(sample_id[j+1], 1, 3))
  print(1)
}

# summary technical replicates
summary(c(1, 2, 3))
