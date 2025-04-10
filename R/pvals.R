##' Calculate p-values adjusting for changes in sign as described by Nixon
##' et al. (2024) in Beyond Normalization: Incorporating Scale Uncertainty
##' in Microbiome and Gene Expression Analysis (internal only)
##'
##' @param p.lower A P x D x S matrix for P covariates, D taxa/genes, and S
##'   monte carlo samples representing the lower tail p.values
##' @param p.lower A P x D x S matrix for P covariates, D taxa/genes, and S
##'   monte carlo samples representing the upper tail p.values
##' @param p.adjust.method An adjutment method for p.adjust
##' @return A list with P x D matrices with the non-adjusted and adjusted
##'   p-values.
##' @author Justin Silverman, Kyle McGovern
aldex.pvals <- function(p.lower, p.upper, p.adjust.method) {

  p.lower <- array(pmin(1, 2 * p.lower), dim = dim(p.lower))
  p.upper <- array(pmin(1, 2 * p.upper), dim = dim(p.upper))

  p.lower.adj <- apply(p.lower, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.upper.adj <- apply(p.upper, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })

  p.lower.mean <- apply(p.lower, c(2,1), mean)
  p.upper.mean <- apply(p.upper, c(2,1), mean)
  rm(p.lower, p.upper)

  p.res <- c()
  for(col_i in 1:ncol(p.lower.mean)) {
    tmp_mat <- cbind(p.lower.mean[,col_i],
                     p.upper.mean[,col_i])
    p.res <- rbind(p.res, apply(tmp_mat, 1, min))
  }
  rm(p.lower.mean, p.upper.mean)

  p.lower.mean.adj <- apply(p.lower.adj, c(1,2), mean)
  p.upper.mean.adj <- apply(p.upper.adj, c(1,2), mean)
  p.adj.res <- c()
  for(col_i in 1:ncol(p.lower.mean.adj)) {
    tmp_mat <- cbind(p.lower.mean.adj[,col_i],
                     p.upper.mean.adj[,col_i])
    p.adj.res <- rbind(p.adj.res, apply(tmp_mat, 1, min))
  }
  rm(p.lower.mean.adj, p.upper.mean.adj)
  return(list(p.res=p.res, p.adj.res=p.adj.res))
}
