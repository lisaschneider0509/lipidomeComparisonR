### load packages
library(tidyverse)
library(viridis) # colorblind save color schemes
library(mixOmics)
library(ggpubr)

source("R/lipidome_comparison_EDA.R")
source("R/lipidome_comparison_plsDA.R")

# set ggplot theme
my_theme <- theme_set(
  theme_minimal() +
    theme(plot.title = element_text(size=12, hjust = 0.5, family="AvantGarde"),
          plot.subtitle = element_text(size = 8, hjust = 0.5, family = "AvantGarde", colour = "grey40"),
          axis.text.x = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.text.y = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          axis.title.x = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          axis.title.y = element_text(size = 10, hjust = 0.5, colour = "grey40", family="AvantGarde"),
          legend.text = element_text(size = 8, colour = "grey40", family="AvantGarde"),
          legend.title = element_text(size = 10, colour = "grey40", family="AvantGarde"))
)

############### set variables #############################
project <- "meat"

working_directory <- "/home/lisa/FH/Masterarbeit/LipidomeComparison"
setwd(working_directory)

data_dir <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data"
lipid_list_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_fish_final_raw.csv"
annotation_path <- "/home/lisa/FH/Masterarbeit/LipidomeComparison/data/meat_annotation.csv"
data_matrix_path <- paste("/home/lisa/FH/Masterarbeit/LipidomeComparison/data/", project, "_data_matrix.csv", sep = "")

plot_path <- paste(working_directory, "/plots", sep = "")


############## import lipid data #############################

lipid_list <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipid_list) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)

############ separate numeric data from group data ################
X <- dplyr::select_if(lipid_list, is.numeric)
Y <- lipid_list$Group
summary(Y)
dim(X); length(Y)

########### SPLS-DA without tuning #####
meat_pls <- mixOmics::splsda(X, Y, keepX = c(50,50)) # SPLS-DA without training
meat_pls_wo_training <- plot_pls_scores(meat_pls, title = "PLS-DA of meat")
background <- mixOmics::background.predict(meat_pls, comp.predicted=2, dist = "max.dist")
mixOmics::plotIndiv(meat_pls, comp = 1:2, group = Y,
          ind.names = FALSE,
          title = "Maximum distance",
          legend = TRUE,
          background = background,
          col = viridis(n = 3),
          style = "ggplot2" ,
          pch = c(15:17),
          cex = 1.5
          )


### tuning the SPLS model ###
set.seed(2543) # for reproducibility 
list.keepX <- c(5:10, seq(5, 50, 2))
## may take a couple of minuites to run depending on the data 
tune_splsda <- mixOmics::tune.splsda(X, Y, 
                           ncomp = 5, 
                           validation = 'Mfold', 
                           folds = 2,            # <= number of observations in the smallest group
                           progressBar = FALSE,  
                           dist = 'max.dist',
                           test.keepX = list.keepX, 
                           nrepeat = 50, # repeats of crossvalidation
                           cpus = 2) # faster computation when using the cpu argument

choice.ncomp <- tune_splsda$choice.ncomp$ncomp; choice.ncomp
choice.keepX <- tune_splsda$choice.keepX[1:choice.ncomp]; choice.keepX
plot(tune_splsda, col = viridis(n = 5, alpha = 0.5))

### tuned SPLS-DA model ###
new_splsda <- mixOmics::splsda(X, Y, ncomp = choice.ncomp, keepX = choice.keepX)

### tuned spls-da performance ###
perf_new_splsda <- mixOmics::perf(new_splsda, validation = "Mfold", folds = 2, 
                    progressBar = FALSE, auc = TRUE, nrepeat = 100) 

{par(mfrow = c(1, 1), 
     cex.main = 1.2, 
     family = "AvantGarde", 
     col = "grey40", 
     col.lab = "grey40", 
     font.main = 1)
  plot(perf_new_splsda, col = viridis(n = 3))
  title(main = "Performance of tuned PLS-DA")}

## display selected variables
mixOmics::selectVar(new_splsda, comp = 1)$value
mixOmics::selectVar(new_splsda, comp = 2)$value

## roc curve
mixOmics::auc.plsda <- auroc(new_splsda, line.col = viridis(n = 3, alpha = 0.5))

## plot samples in the first two components of PLS-DA
tuned_pls_scores <- plot_pls_scores(new_splsda, 
                title = "Tuned PLS-DA", 
                xlab = "component 1", 
                ylab = "component 2")
ggsave(paste(plot_path, "tuned_pls-da.png", sep ="/" ), 
       tuned_pls_scores,
       device = "png", 
       width = 10, 
       height = 5)

## loading weights of selected variables
loadings1 <- mixOmics::plotLoadings(new_splsda, 
             comp = 1, 
             title = 'Loadings on comp 1', 
             contrib = 'max', 
             method = 'mean', 
             legend.color = viridis(n = 3), 
             size.title = 1)

loadings2 <- mixOmics::plotLoadings(new_splsda, 
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

loadings <- ggpubr::ggarrange(lp1, lp2, legend = "right", 
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

## visualize the abundances of the selected variables
pls_sub1 <- mixOmics::selectVar(new_splsda, comp = 1)$name
pls_sub2 <- mixOmics::selectVar(new_splsda, comp = 2)$name
pls1 <- subset(lipid_list, select = c("Group", pls_sub1))
pls_all <- cbind(pls1, subset(lipid_list, select = c(pls_sub2)))

parallel_plot <- parallel_plot(pls_all, pls_all$Group)

ggsave(paste(plot_path, "pp_pls_da.pbg", sep ="/" ), 
       parallel_plot,
       device = "png", 
       width = 10, 
       height = 4)

