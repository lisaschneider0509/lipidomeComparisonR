################## Functions for PLS-DA #################################
#' PLS-DA scores plot
#' 
#' `plot_pls_scores` plots the loadings between the first two principal components. 
#' @param pls_da_obj Object produced by MixOmics::splsda or MixOmics::plsda. 
#' @param title string. Main title of the plot. Default = "PLS-DA"
#' @param colour logical. TRUE ... color scale by PC1. FALSE (default) ... all black. 
#' @param top_loadings integer. Number of loadings to label (the top n loadings of both PCs are labelled).
#' @param out_path string. Path to save plot as png. If "none" (default), plot is only printed to device. 
#' @example 
#' data(srbct)
#' srbc_pls <- plsda(srbct$gene, srbct$class)
#' plot_pls_scores(srbc_pls)
#' plot_pls_scores(srbc_pls, colour = FALSE)
plot_pls_scores <- function(pls_da_obj,
                            colour = TRUE,
                            top_loadings = 10,
                            title = "PLS-DA",
                            xlab = "PLS-DA 1",
                            ylab = "PLS-DA 2"){
  X <- as.data.frame(pls_da_obj$variates$X)
  
  pc1_max <- max(abs(X$comp1))
  pc2_max <- max(abs(X$comp2))
  
  if(colour == TRUE){
    loadings_plot <- pls_plot <- ggplot(X, aes(x=comp1, 
                                               y=comp2, 
                                               color = pls_da_obj$Y, 
                                               fill = pls_da_obj$Y, 
                                               shape = pls_da_obj$Y)) + 
      geom_point(size=1)  +
      stat_ellipse(type = "norm", geom = "polygon")
  }
  else{ loadings_plot <- pls_plot <- ggplot(X, aes(x=comp1, 
                                                   y=comp2)) +
    geom_point(size=1)}
  
  pls_plot <- pls_plot + 
    geom_point(size= 1.5)  +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               colour = "grey60") +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               colour = "grey60") +
    labs(title = title, x = xlab, y = ylab) +
    
    scale_color_viridis_d() +
    scale_fill_viridis_d(alpha = 0.1) +
    # coord_cartesian(xlim = c(-pc1_max, pc1_max), ylim = c(-pc2_max, pc2_max)) +
    # xlim(-pc1_max, pc1_max) +
    # ylim(-pc2_max, pc2_max) +
    theme(legend.title = element_blank())
  my_theme
  
  pls_plot
}
