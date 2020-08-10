library(mixOmics)
source("lipidome_comparison_plsDA.R")
set.seed(1)
train_ind <- sample(seq_len(nrow(beef)), size = 0.75*nrow(beef))
train_data <- beef[train_ind, ]
test_data <- beef[-train_ind, ]


X <- select_if(train_data, is.numeric)
Y <- train_data$nutrition
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
          col = viridis(n = 2), 
          style = "ggplot2" , 
          pch = c(16:17), 
          cex = 1.5
)


### spls-da ###
list.keepX <- c(1:10, seq(10, 300, 5))
set.seed(2543) # for reproducibility here,
# to speed up the computational time, consider the cpu argument
# take ~ 4 min to run
tune_splsda <- tune.splsda(X, Y, ncomp = 2, validation = 'Mfold', folds = 2, 
                           progressBar = FALSE, dist = 'max.dist',
                           test.keepX = list.keepX, nrepeat = 10, cpus = 2) #nrepeat 50-100 for better estimate

choice.ncomp <- tune_splsda$choice.ncomp$ncomp + 1; choice.ncomp
choice.keepX <- tune_splsda$choice.keepX[1:choice.ncomp]; choice.keepX
plot(tune_splsda, col = viridis(n = 2, alpha = 0.5))
new_splsda <- splsda(X, Y, ncomp = choice.ncomp, keepX = choice.keepX)
perf_splsda2 <- perf(new_splsda, validation = "Mfold", folds = 2, 
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) 

{par(mfrow = c(1, 1), cex.main = 1.2, family = "AvantGarde", col = "grey40", col.lab = "grey40", font.main = 1)
  plot(perf_splsda2, col = viridis(n = 3))
  title(main = "Performance of tuned PLS-DA")}
selectVar(new_splsda, comp = 1)$value
auc.plsda <- auroc(new_splsda, line.col = viridis(n = 1, alpha = 0.5))


tuned_pls <- plot_pls_scores(new_splsda, 
                             title = "Tuned PLS-DA", 
                             xlab = "component 1", 
                             ylab = "component 2")
ggsave(paste(plot_path, "tuned_pls-da.png", sep ="/" ), 
       tuned_pls,
       device = "png", 
       width = 7, 
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



lp1 <- ggplot(loadings1, aes(reorder(rownames(loadings1), -abs(importance), sum), importance, color = GroupContrib, fill = GroupContrib))+
  geom_col() +
  coord_flip() +
  scale_color_viridis_d(guide = F) +
  scale_fill_viridis_d(guide = F) +
  labs(title = "Loadings of component 1", x = "feature") +
  theme(legend.title = element_blank(), legend.position = "none") 

lp2 <- ggplot(loadings2, aes(reorder(rownames(loadings2), -abs(importance), sum), importance, color = GroupContrib, fill = GroupContrib))+
  geom_col() +
  coord_flip() +
  scale_color_viridis_d(begin = 1, guide = F) + 
  scale_fill_viridis_d(begin = 1, guide = F) +
  labs(title = "Loadings of component 2", x = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

predict_meat <- predict(new_splsda, XT, YT, ncomp = 2)
prediction <- predict_meat$class$max.dist[,2]
confusion.mat <- get.confusion_matrix(truth = YT, predicted = prediction)
get.BER(confusion.mat)

caret::createFolds()

pls_sub1 <- selectVar(new_splsda, comp = 1)$name
pls_sub2 <- selectVar(new_splsda, comp = 2)$name
pls1 <- subset(beef, select = c("nutrition", pls_sub1))
pls1 <- cbind(pls1, subset(beef, select = c(pls_sub2)))


pp <- parallel_plot(pls1, pls1$nutrition)
loadings <- ggarrange(lp1, lp2, pp, legend = "right", 
                      labels = c("(A)", "(B)", "(C)"), 
                      widths = c(0.8, 0.8, 1), 
                      font.label = list(size = 10, 
                                        color = "grey40", 
                                        face = "plain", 
                                        family = "AvantGarde"), ncol = 3)
ggsave(paste(plot_path, "pls_da_loadings.png", sep ="/" ), 
       loadings,
       device = "png", 
       width = 10, 
       height = 4)
