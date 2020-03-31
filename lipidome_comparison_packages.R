### Install packages
## general
# install.packages("tibble")
# install.packages("stringr")
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("devtools")
# install.packages(DT) # pretty print df (interactive)

## graphs
# install.packages("ggplot2")
# install.packages("scales")
# install.packages("viridis")
# install.packages("ggrepel")

## correlation plot 
# install.packages("psych")

## for spider chart
# install.packages("fmsb")

## for paralell plot
# install.packages("GGally")

# install.packages("hrbrthemes")

## for PCA
# install.packages("ggfortify")
# install.packages("factoextra")

## clustering 
# install.packages("plotly") # interactive heatmap
# install.packages("gplots")
# install.packages("heatmaply")
# install.packages("dendextend") # find the best clustering parameters

## hypothesis testing
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano', 'limma')
BiocManager::install("airway")
BiocManager::install('DESeq2')
BiocManager::install("ALL")
BiocManager::install("a4Base")

## preprocessing
# install.packages("imputeLCMD") # impute missing values (for left censored missing values)
