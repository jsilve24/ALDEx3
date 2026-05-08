#' Simple plots for an ALDEx3  object
#'
#' Provides volcano, effect, MA and waterfall plots from an ALDEx3  object.
#'
#' This method plots combinations of adjusted p-values from `object$p.val.adj`,
#' posterior estimates, standard errors and log abundance values
#' averaged across Monte Carlo samples.
#' The result is returned as a single plots of the desired type
#' Note-calls the aldex `summary` and `cohensd` functions internally
#' Note-the contrast must be supplied
#' Note-the data for MA plots may not be available for very large datasets
#' Note-if there are many tied values the waterfall function will show the first n features
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
#' @importFrom graphics abline barplot mtext points text
#' @export
aldex.plot <-function(object, plot=c("volcano", "effect", "MA", "water"), 
  contrast=NULL, threshold=0.05, min.diff=0.5, cohen=0.5, sig.col=rgb(1,0,0,0.5),
  water.show=5, water.col=c("red", "blue"), water.names=TRUE,  ... ){
    
  	# this allows partial matching because I'm lazy; defaults to volcano
    plot=match.arg(plot)
    if( is.null(contrast) ) stop("\nPlease enter the name of the contrast to plot")
    
  	# check names in the passed object 
  	# throws an error if the needed slot is not populated
    # only needed for MA plot, so test is inside that if statement  
  	# if(!("logScale" %in% names(object))) stop("\nlogScale slot not found\ntry reducing nsample when running aldex()\nnsamples=32 is a reasonable minimum")

  	#contrast = names(object[[9]])
  	#contrast = names(object$data)
 	
  	nms0 <- as.character(object$data[,contrast])[as.numeric(which(object$X[2,] == 0)[1])]
  	nms1 <- as.character(object$data[,contrast])[as.numeric(which(object$X[2,] == 1)[1])]
  	
  	# 0 should be the negative direction for the estimate

    # call to summary() 
    # this will need to be updated for glm with multiple contrasts
    sum.output <- summary(object)
    nsamples <- length(object$data[,1])
    # get the sig features
    sig <- sum.output$p.val.adj < threshold
    
    if(plot=="volcano"){
      # replace 0 with min pval/10
      p.val1 <- sum.output$p.val.adj > 0
      p.val0 <- sum.output$p.val.adj == 0
      min.p <- min(sum.output$p.val.adj[p.val1])
      sum.output$p.val.adj[p.val0] <- min.p/10
      y.val <- -(log10(sum.output$p.val.adj))
      plot(sum.output$estimate, y.val, pch=19, col=rgb(0,0,0,0.3),
        xlab="estimate", ylab="-log10(p.adjust)", ...)
      abline(h=-log10(threshold), lty=2)
      mtext(nms0, side=1, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(nms1, side=1, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(sum.output$estimate[sig], y.val[sig], pch=19, col=sig.col, ...)
    }else if(plot=="effect"){
      cohens <- cohensd(object, contrast)
      plot(cohens$pooled.SD, sum.output$estimate,
      	pch=19, col=rgb(0,0,0,0.3),xlab="std dev", ylab='estimate', ...)
      mtext(nms0, side=2, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(nms1, side=2, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
    
      points(cohens$pooled.SD[sig], sum.output$estimate[sig],
      	pch=19, col=sig.col, ... )
	  abline(0,cohen, lty=2, col='grey') 
	  abline(0,-cohen, lty=2, col='grey')
    }else if(plot=="MA"){
    ####
    # currently logComp and logScale can be too large to print out
    # in this case there is an error and no plot
    # example, yeast dataset with nsample > 128(ish)
    # there is a known hack in aldex()
    ####
    if(!("logScale" %in% names(object))) stop("\nlogScale slot not found\ntry reducing nsample when running aldex()\nnsamples=32 is a reasonable minimum")
      vals <- vector()
      for(i in 1:length(object$logComp[,1,1])){vals[i]= mean(object$logComp[i,,])}
      plot(vals,sum.output$estimate, pch=19, col=rgb(0,0,0,0.3),
        xlab='log abundance', ylab='estimate', ... )
      mtext(nms0, side=2, line = 2, at = min(sum.output$estimate),
            col = "grey", cex = 0.8)
      mtext(nms1, side=2, line = 2, at = max(sum.output$estimate),
            col = "grey", cex = 0.8)
      points(vals[sig], sum.output$estimate[sig], pch=19, col=sig.col, ... )
      abline(h=min.diff, lty=2, col='grey')
      abline(h=-min.diff, lty=2, col='grey')
    }else if(plot=="water"){
	  # positions below threshold
	  pos <- which(sum.output$p.val.adj < threshold & sum.output$estimate > min.diff)    
	  neg <- which(sum.output$p.val.adj < threshold & sum.output$estimate < -min.diff)  

	  # order from lowest to highest adjusted p value and then by estimate
	  pos.p <-  pos[order(sum.output[pos, "p.val.adj"], sum.output[pos, "estimate"], decreasing=c(T,F))]
	  neg.p <-  neg[order(sum.output[neg, "p.val.adj"], sum.output[neg, "estimate"])]
	  
	  # select n values to display
	  display = water.show
	  display.p <- display
	  display.n <- display
      # if there are a more 0 values than n to display order them by estimate
	  if(length(which(sum.output$p.val.adj == min(sum.output$p.val.adj)) > display)){
	  	
	  }

	  if(length(pos.p) <= display) display.p <- length(pos.p)
	  if(length(neg.p) <= display) display.n <- length(neg.p)
	  
	  pos.p <- pos.p[1:display.p]
	  neg.p <- neg.p[1:display.n]

	  # get L2FC values and order from largest to smallest 
	  pos.e <-  pos.p[order(sum.output[pos.p, "estimate"], decreasing=TRUE)]
	  neg.e <-  rev(neg.p[order(sum.output[neg.p, "estimate"])])
	  
	  min.estimate <- min(sum.output[c(pos.e, neg.e), "estimate"])
	  
	  plot.offset <- barplot(sum.output[c(pos.e, neg.e), "estimate"], horiz=T, las=2, col=c(rep(water.col[1], display.p), rep(water.col[2],display.n)), xlab="estimate", ...)
	  if(water.names==TRUE) {
	    text(c(rep(min.estimate,display.p), rep(0.5,display.n)), plot.offset[,1], 
	      adj=0, sum.output[c(pos.e, neg.e), "entity"])
	  }
    }
  
}