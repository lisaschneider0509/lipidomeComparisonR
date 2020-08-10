library(mixOmics)
### load packages
{library(dplyr) # select part of data
  library(stringr) # count separators
  library(data.table) # transpose data frame
  # library(crmn)
  # library(impute)
  # library(imputeLCMD)
  library(tibble) # data frame manipulation
  
  library(tidyverse)
  library(ggplot2)#, # plots
  library(viridis) # colorblind save color schemes
  # library(GGally) # paralell plot
  # library(fmsb) # spider chart
  # library(scales) # scale opacity of filling (alpha)
  # library(ggpubr) # multiple plots on one page
  # library(ggmosaic)
  
  # library(ggrepel)
  # library(factoextra)
  # library(ggfortify) # biplot with ggplot
  # library(corrplot)
  # library(FactoMineR)
  
  library(heatmaply) # interactive heatmap
  library(gplots) # heatmap
  library(plotly) # interactive ggplots
  library(htmlwidgets) # save plotly-plots as html
  library(dendextend)
  
  source("R/lipidome_comparison_dataTransformaions.R")
  source("R/lipidome_comparison_EDA.R")
  source("R/lipidome_comparison_pca.R")
  source("R/lipidome_comparison_clustering.R")
  source("R/lipidome_comparison_hypothesis_testing.R")}

  source("R/lipidome_comparison_plsDA.R")

############### set variables #############################
{
  project <- "meat"
  
  working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
  setwd(working_directory)
  
  data_dir <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data"
  lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
  annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"
  data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/LipidomeComparison/data/", project, "_data_matrix.csv", sep = "")
  
  plot_path <- paste(working_directory, "/plots", sep = "")
}

############## import lipid data #############################

lipid_list <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipid_list) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)


set.seed(1)
train_ind <- sample(seq_len(nrow(lipid_list)), size = 0.75*nrow(lipid_list))
train_data <- lipid_list#[train_ind, ]
test_data <- lipid_list[-train_ind, ]


X <- select_if(train_data, is.numeric)
Y <- train_data$Group
summary(Y)
dim(X); length(Y)

XT <- select_if(test_data, is.numeric)
YT <- test_data$group

meat_pls <- splsda(X, Y, keepX = c(50,50)) # 1 Run the method
x <- plot_pls_scores(meat_pls, title = "PLS-DA of meat - training data")
background <- background.predict(meat_pls, comp.predicted=2, dist = "max.dist")
plotIndiv(meat_pls, comp = 1:2, group = Y,
          ind.names = FALSE, 
          title = "Maximum distance",
          legend = TRUE,  
          background = background, 
          col = viridis(n = 3), 
          style = "ggplot2" , 
          pch = c(15:17), 
          cex = 1.5
          )


### spls-da ###
list.keepX <- c(5:10, seq(5, 50, 2))
set.seed(2543) # for reproducibility here,
# to speed up the computational time, consider the cpu argument
# take ~ 4 min to run
tune_splsda <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 2, 
                           progressBar = FALSE, dist = 'max.dist',
                           test.keepX = list.keepX, nrepeat = 50, cpus = 2) #nrepeat 50-100 for better estimate

choice.ncomp <- tune_splsda$choice.ncomp$ncomp; choice.ncomp
choice.keepX <- tune_splsda$choice.keepX[1:choice.ncomp]; choice.keepX
plot(tune_splsda, col = viridis(n = 5, alpha = 0.5))
new_splsda <- splsda(X, Y, ncomp = choice.ncomp, keepX = choice.keepX)
perf_splsda2 <- perf(new_splsda, validation = "Mfold", folds = 2, 
                    progressBar = FALSE, auc = TRUE, nrepeat = 100) 

{par(mfrow = c(1, 1), cex.main = 1.2, family = "AvantGarde", col = "grey40", col.lab = "grey40", font.main = 1)
  plot(perf_splsda2, col = viridis(n = 3))
title(main = "Performance of tuned PLS-DA")}
selectVar(new_splsda, comp = 1)$value
auc.plsda <- auroc(new_splsda, line.col = viridis(n = 3, alpha = 0.5))


tuned_pls <- plot_pls_scores(new_splsda, 
                title = "Tuned PLS-DA", 
                xlab = "component 1", 
                ylab = "component 2")
ggsave(paste(plot_path, "tuned_pls-da.png", sep ="/" ), 
       tuned_pls,
       device = "png", 
       width = 10, 
       height = 5)

loadings1 <- plotLoadings(new_splsda, 
             comp = 1, 
             title = 'Loadings on comp 1', 
             contrib = 'max', 
             method = 'mean', 
             legend.color = viridis(n = 3), 
             size.title = 1)

loadings2 <- plotLoadings(new_splsda, 
             comp = 2, 
             title = 'Loadings on comp 2', 
             contrib = 'max',
             method = 'mean', 
             legend.color = viridis(n = 3), 
             size.title = 1)



lp1 <- ggplot(loadings1, aes(reorder(rownames(loadings1), importance, sum), importance, color = GroupContrib, fill = GroupContrib))+
  geom_col() +
  coord_flip() +
  scale_color_viridis_d(begin = 0.5) +
  scale_fill_viridis_d(begin = 0.5) +
  labs(title = "Loadings of component 1", x = "feature") +
  theme(legend.title = element_blank()) 

lp2 <- ggplot(loadings2, aes(reorder(rownames(loadings2), importance, sum), importance, color = GroupContrib, fill = GroupContrib))+
    geom_col() +
    coord_flip() +
    scale_color_viridis_d(begin = 1) + 
    scale_fill_viridis_d(begin = 1) +
    labs(title = "Loadings of component 2", x = NULL) +
    theme(legend.title = element_blank())

predict_meat <- predict(new_splsda, XT, YT, ncomp = 2)
prediction <- predict_meat$class$max.dist[,2]
confusion.mat <- get.confusion_matrix(truth = YT, predicted = prediction)
get.BER(confusion.mat)


pls_sub1 <- selectVar(new_splsda, comp = 1)$name
pls_sub2 <- selectVar(new_splsda, comp = 2)$name
pls1 <- subset(lipid_list, select = c("group", pls_sub1))
pls1 <- cbind(pls1, subset(lipid_list, select = c(pls_sub2)))


pp <- parallel_plot(pls1, pls1$group)
loadings <- ggarrange(lp1, lp2, legend = "right", 
                      labels = c("(A)", "(B)"), 
                      font.label = list(size = 10, 
                                        color = "grey40", 
                                        face = "plain", 
                                        family = "AvantGarde"), ncol = 2, common.legend = TRUE)
ggsave(paste(plot_path, "pls_da_loadings.png", sep ="/" ), 
       loadings,
       device = "png", 
       width = 10, 
       height = 4)

ggsave(paste(plot_path, "pp_pls_da.pbg", sep ="/" ), 
       pp,
       device = "png", 
       width = 10, 
       height = 4)

