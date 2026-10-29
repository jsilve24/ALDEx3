##' Combine output from aldex streams (internal only)
##'
##' Takes a list l, where every element is itself a list of 3D arrays. Each
##' array should have matching first and second dimentions. This function
##' combines those arrays along the third dimension.
##' @param l a list, where every element is itself a list of 3D arrays.
##' @return a list of 3D arrays
##' @author Justin Silverman
##' @keywords internal
combine.streams <- function(l) {
   fields <- names(l[[1]])
   stats::setNames(lapply(fields, function(field) {
     abind::abind(lapply(l, function(chunk) chunk[[field]]), along=3)
   }), fields)
}



##' Test object contains elements
##'
##' @param obj object (list type) 
##' @param names names of elements required to be in the list 
##' @author Justin Silverman
##' @keywords internal
req <- function(obj, names){
  present <- sapply(obj[names], is.null)
  if(any(present)){
    stop("object does not contain required components:", names[present])
  }
}
