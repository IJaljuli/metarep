summary.metarepl <- function(object, ...) {
  
  meta:::chkclass(object, "metarepl")
  # chkclass(object, "metarepl")
  
  res <- summary.meta(object)
  ##
  
  res$r.value <- object$r.value
  if(!is.null(object$u_R) ) res$u.increased <- object$u_R
  if(!is.null(object$u_L) ) res$u.decreased <- object$u_L
  
  ##
  class(res) <- c("summary.metarepl", class(res))
  ##
  res
}
