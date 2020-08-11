################## Functions for PLS-DA #################################
#' PLS-DA scores plot
#'
#' `plot_pls_scores` plots the loadings between the first two principal components.
#' @param pls_da_obj Object produced by MixOmics::splsda or MixOmics::plsda.
#' @param title string. Main title of the plot. Default = "PLS-DA"
#' @param colour logical. TRUE ... color scale by PC1. FALSE (default) ... all black.
#' @param top_loadings integer. Number of loadings to label (the top n loadings of both PCs are labelled).
#' @param xlab string. Title of x-axis. Default: "PLS-DA 1"
#' @param ylab string. Title of y-axis. Default: "PLS-DA 2"
#' @export
#' @example
#' {srbct <- data(srbct)
#' srbc_pls <- mixOmics::plsda(srbct$gene, srbct$class)
#' plot_pls_scores(srbc_pls)
#' plot_pls_scores(srbc_pls, colour = FALSE)}
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
    pls_plot <- ggplot2::ggplot(X, ggplot2::aes(x=comp1,
                                            y=comp2,
                                            color = pls_da_obj$Y,
                                            fill = pls_da_obj$Y,
                                            shape = pls_da_obj$Y)) +
      ggplot2::geom_point(size=1)  +
      ggplot2::stat_ellipse(type = "norm", geom = "polygon")
  }
  else{ pls_plot <- ggplot2::ggplot(X, ggplot2::aes(x=comp1,
                                                y=comp2)) +
    ggplot2::geom_point(size=1)}

  pls_plot <- pls_plot +
    ggplot2::geom_point(size= 1.5)  +
    ggplot2::geom_hline(yintercept = 0,
               linetype = "dashed",
               colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0,
               linetype = "dashed",
               colour = "grey60") +
    ggplot2::labs(title = title, x = xlab, y = ylab) +

    viridis::scale_color_viridis(discrete = TRUE) +
    viridis::scale_fill_viridis(alpha = 0.1, discrete = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=12, hjust = 0.5),
          axis.text.x = ggplot2::element_text(size = 8),
          axis.title.x = ggplot2::element_text(size = 10, hjust = 0.5),
          axis.title.y = ggplot2::element_text(size = 10, hjust = 0.5),
          legend.text = ggplot2::element_text(size = 8),
          legend.title = ggplot2::element_blank())

  pls_plot
}
