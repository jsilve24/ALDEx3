##' Combine output from aldex streams (internal only)
##'
##' Takes a list l, where every element is itself a list of 3D arrays. Each
##' array should have matching first and second dimentions. This function
##' combines those arrays along the third dimension.
##' @param l a list, where every element is itself a list of 3D arrays.
##' @return a list of 3D arrays
##' @author Justin Silverman, Kyle McGovern
combine.streams <- function(out) {
  out.names <- names(out[[1]])
  out <- lapply(out.names, function(name) {
    ## TODO this is a bit hacky
    if(length(dim(out[[1]][[name]]))>2) {
      do.call(abind::abind, c(lapply(out, function(x) x[[name]]), list(along = 3)))
    } else {
      do.call(cbind, c(lapply(out, function(x) x[[name]])))
    }

  })
  names(out) <- out.names
  out
}
