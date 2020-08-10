library(lipidr)
library(dplyr)
library(fgsea)
source("lipidome_comparison_hypothesis_testing.R")


lipidr_list <- read.csv(file = paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), row.names = 1)
colnames(lipidr_list) <- base::unlist(read.csv(paste(data_dir, "/", project, "_renamed_data.csv", sep = ""), header = FALSE)[1, -1], use.names = FALSE)
lipidr_list <- subset(lipidr_list, subset = lipidr_list$Group == "beef")

lipidr_matrix <- t(select_if(lipidr_list, is.numeric))
colnames(lipidr_matrix) <- lipidr_list$Replicate

annotation_df <- subset(lipidr_list, select = c("Replicate", "Sample"))
annotation_data <- read.csv(annotation_path, sep = ",", dec = ".", header = TRUE)
map <- data.frame(Sample = annotation_data$Sample, 
                  Nutrition = annotation_data$Nutrition)
annotation_df <- left_join(annotation_df, map, by = "Sample")
names(annotation_df) <- c("Sample", "ID", "Nutrition")

lipidr_experiment <- as_lipidomics_experiment(lipidr_matrix)
lipidr_experiment <- add_sample_annotation(lipidr_experiment, annotation_df)

plot_samples(lipidr_experiment, 
             type = "tic", 
             log = TRUE)
plot_molecules(lipidr_experiment, 
               "sd")
plot_lipidclass(lipidr_experiment, "boxplot")
lipidr_norm <- normalize_pqn(lipidr_experiment, 
              measure = "Area", 
              log = TRUE)
plot_samples(lipidr_norm, "boxplot")

lipidr_beef <- lipidr_norm[, lipidr_norm$Group == "beef"]

pca_results <- mva(lipidr_experiment, 
                   measure = "Area", 
                   method = "PCA")
plot_mva(pca_results, color_by = "Group", components = c(1, 2)) +
  theme_minimal() +
  scale_color_viridis_d() + 
  scale_fill_viridis_d()

oplsda_results <- mva(lipidr_beef, 
                      method = "OPLS-DA", 
                      group_col = "Nutrition", 
                      groups = c("grazing", "stall"))
plot_mva(oplsda_results, color_by="Nutrition") +
  theme_minimal() +
  scale_color_viridis_d() + 
  scale_fill_viridis_d()
lipidr_experiment

de_results = de_analysis(
  data=lipidr_beef, 
  group_col = "Nutrition",
  stall - grazing,
  measure="Area"
)
head(de_results)
plot_results_volcano(de_results, show.labels = T) +
  scale_color_viridis_d()

vp <- volcano_plot(volcano_df = de_results, foldchange_col = de_results$logFC, significance_col = de_results$adj.P.Val)
ggsave(filename = paste(plot_path, "beef_volcano.png", sep = "/"), 
       plot = vp, 
       device = "png", 
       width = 7, 
       height = 5
)

x <- subset(de_results$Molecule, de_results$adj.P.Val <= 0.05)


enrich_results <- lsea(de_results, rank.by = "logFC")
sig_lipi <- significant_lipidsets(enrich_results)
plot_class_enrichment(de_results, significant_lipidsets(enrich_results))

sig_lipi <- sub(".*_", "", sig_lipi$`stall - grazing`)
lc <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = Class, y = logFC, fill = is.element(de_results$Class, sig_lipi))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by lipid classes", y = "log2FC", x = "lipid class")

tc <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = as.factor(total_cl), y = logFC, fill = is.element(de_results$total_cl, sig_lipi))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by total chain length", y = "log2FC", x = "total chain length")

tus <- ggplot(data = de_results) + 
  geom_boxplot(aes(x = as.factor(total_cs), y = logFC, fill = is.element(de_results$total_cs, sig_lipi))) +
  scale_fill_manual(name = 'Significant', values = setNames(c(viridis(n = 1, begin = 0.5),'grey90'),c(T, F))) +
  labs(title = "Foldchange by unsaturated bonds", y = "log2FC", x = "total unsaturated bonds")


mylist <- list(lc, tc, tus)

lsea <- ggarrange(plotlist = mylist, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", 
                  labels = c("(A)", "(B)", "(C)"), 
                  font.label = list(size = 10, 
                                    color = "grey40", 
                                    face = "plain", 
                                    family = "AvantGarde"))
ggsave(filename = paste(plot_path, "lsea.png", sep = "/"), plot = lsea, device = "png", 
       width = 15, height = 5)
