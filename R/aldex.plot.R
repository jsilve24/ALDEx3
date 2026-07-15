#' Simple plots for an ALDEx3  object
#'
#' Provides volcano, effect, MA and waterfall plots from an ALDEx3  object.
#'
#' This method plots combinations of adjusted p-values from `object$p.val.adj`,
#' posterior estimates, standard errors and log abundance values
#' averaged across Monte Carlo samples.
#' The result is returned as a single plot of the desired type. Effect plots
#' use the ALDEx2-inspired diagnostics returned by `aldex.effect`. The contrast
#' must identify exactly one binary model coefficient. Effect diagnostics and
#' MA plots require `logComp` and `logScale`, which are not returned when
#' `aldex` streams large jobs.
#'
#' @title Plot Method for ALDEx3 objects with pairwise comparisons
#' @param object An object of class \code{aldex}
#' @param plot type of plot (default='volcano')
#' @param contrast the name of the comparison to plot, must be provided
#' @param threshold FDR significance threshold (default=0.05)
#' @param min.diff (default=0.5) used for MA (display), and waterfall (cutoff) plots
#' @param cohen Cohen's d threshold for effect plot (slope)
#' @param sig.col color for significant features
#' @param water.show number features to show in waterfall plot
#' @param water.col vector of colors for waterfall plot
#' @param water.names logical, show names of features for waterfall plot 
#' @param ... additional plot parameters (cex)
#' @return the desired plot
#' @author Greg Gloor
#' @importFrom grDevices rgb
#' @importFrom graphics abline barplot mtext plot.new points text
#' @export
aldex.plot <-function(object, plot=c("volcano", "effect", "MA", "water", "3d"), 
  contrast=NULL, threshold=0.05, min.diff=0.5, cohen=0.5, sig.col=rgb(1,0,0,0.5),
  water.show=5, water.col=c("red", "blue"), water.names=TRUE,  ... ){
    
    plot=match.arg(plot)
    if( is.null(contrast) ) stop("\nPlease enter the name of the contrast to plot")

    contrast.idx <- .aldex_resolve_contrast(object, substitute(contrast))
    parameter <- rownames(object$X)[contrast.idx]
    contrast.labels <- .aldex_contrast_labels(object, contrast, contrast.idx)

    sum.output <- summary(object)
    sum.output <- sum.output[sum.output$parameter == parameter,, drop=FALSE]
    if (nrow(sum.output) == 0) stop("summary output does not contain contrast: ", parameter)

    sig <- sum.output$p.val.adj < threshold
    
    if(plot=="volcano"){
      p.val1 <- sum.output$p.val.adj > 0
      p.val0 <- sum.output$p.val.adj == 0
      min.p <- min(sum.output$p.val.adj[p.val1], na.rm=TRUE)
      if (!is.finite(min.p)) min.p <- .Machine$double.xmin
      sum.output$p.val.adj[p.val0] <- min.p/10
      y.val <- -(log10(sum.output$p.val.adj))
      plot(sum.output$estimate, y.val, pch=19, col=rgb(0,0,0,0.3),
        xlab="estimate", ylab="-log10(p.adjust)", ...)
      abline(h=-log10(threshold), lty=2)
      mtext(contrast.labels[1], side=1, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(contrast.labels[2], side=1, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(sum.output$estimate[sig], y.val[sig], pch=19, col=sig.col, ...)
    }else if(plot=="effect"){
      effect.output <- aldex.effect(object, contrast)
      plot(effect.output$pooled.SD, effect.output$estimate,
      	pch=19, col=rgb(0,0,0,0.3),xlab="std dev", ylab='estimate', ...)
      mtext(contrast.labels[1], side=2, line = 2, at = min(effect.output$estimate),
            col = "grey", cex = 0.8)
      mtext(contrast.labels[2], side=2, line = 2, at = max(effect.output$estimate),
            col = "grey", cex = 0.8)
    
      points(effect.output$pooled.SD[sig], effect.output$estimate[sig],
      	pch=19, col=sig.col, ... )
	  abline(0,cohen, lty=2, col='grey') 
	  abline(0,-cohen, lty=2, col='grey')
    }else if(plot=="MA"){
      if(!all(c("logComp", "logScale") %in% names(object))) {
        stop("\nlogComp/logScale slots not found\ntry reducing streamsize or nsample when running aldex()\nnsample=32 is a reasonable minimum")
      }
      logW <- sweep(object$logComp, c(2,3), object$logScale, FUN = `+`)
      vals <- rowMeans(logW, dims=1)
      plot(vals,sum.output$estimate, pch=19, col=rgb(0,0,0,0.3),
        xlab='log abundance', ylab='estimate', ... )
      mtext(contrast.labels[1], side=2, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(contrast.labels[2], side=2, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(vals[sig], sum.output$estimate[sig], pch=19, col=sig.col, ... )
      abline(h=min.diff, lty=2, col='grey')
      abline(h=-min.diff, lty=2, col='grey')
    }else if(plot=="water"){
	  pos <- which(sum.output$p.val.adj < threshold & sum.output$estimate > min.diff)    
	  neg <- which(sum.output$p.val.adj < threshold & sum.output$estimate < -min.diff)  

      if (length(pos) == 0 && length(neg) == 0) {
        plot.new()
        text(0.5, 0.5, "No features pass the plotting thresholds")
        return(invisible(NULL))
      }

	  pos.p <-  pos[order(sum.output[pos, "p.val.adj"], sum.output[pos, "estimate"])]
	  neg.p <-  neg[order(sum.output[neg, "p.val.adj"], sum.output[neg, "estimate"])]
	  
      display.p <- min(water.show, length(pos.p))
      display.n <- min(water.show, length(neg.p))
	  pos.p <- if (display.p > 0) pos.p[seq_len(display.p)] else integer(0)
	  neg.p <- if (display.n > 0) neg.p[seq_len(display.n)] else integer(0)

	  pos.e <-  pos.p[order(sum.output[pos.p, "estimate"], decreasing=TRUE)]
	  neg.e <-  rev(neg.p[order(sum.output[neg.p, "estimate"])])
      selected <- c(pos.e, neg.e)
      if (length(selected) == 0) {
        plot.new()
        text(0.5, 0.5, "No features pass the plotting thresholds")
        return(invisible(NULL))
      }
	  
	  min.estimate <- min(sum.output[selected, "estimate"])
	  
	  plot.offset <- barplot(sum.output[selected, "estimate"], horiz=T, las=2, col=c(rep(water.col[1], display.p), rep(water.col[2],display.n)), xlab="estimate", ...)
	  if(water.names==TRUE) {
	    text(c(rep(min.estimate,display.p), rep(0.5,display.n)), plot.offset[,1], 
	      adj=0, sum.output[selected, "entity"])
	  }
    }
  
}

.aldex_contrast_labels <- function(object, contrast, contrast.idx) {
  x <- object$X[contrast.idx,]
  labels <- c("0", "1")
  if (!is.null(object$data) && is.character(contrast) && contrast %in% colnames(object$data)) {
    labels[1] <- as.character(object$data[[contrast]][which(x == 0)[1]])
    labels[2] <- as.character(object$data[[contrast]][which(x == 1)[1]])
  }
  labels
}
