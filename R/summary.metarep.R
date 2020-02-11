summary.metarep <- function(object, ...) {
  
  meta:::chkclass(object, "metarep")
  # chkclass(object, "metarep")
  
  res <- summary.meta(object)
  ##
  
  res$r.value <- round( object$r.value , digits = 4 )
  if(!is.null(object$u_R) ) res$u.increased <- object$u_R
  if(!is.null(object$u_L) ) res$u.decreased <- object$u_L
  
  ##
  class(res) <- c("summary.metarep", class(res))
  ##
  res
}
