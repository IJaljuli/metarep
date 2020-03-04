#' Print meta-analysis with replicability-analysis results
#' 
#' @description
#' Print method for objects of class \code{metarep}.
#' 
#' 
#' @param x An object of class \code{metarep}
#'
#' @param ... 
#' @import meta
#' @export
#' @examples
#' n.i.1 <- c( 20, 208, 24, 190, 58, 36, 51)
#' a.i <- c( 2,79,0,98,15,34,9) 
#' n.i.2 <- c( 20, 119, 22, 185, 29, 51, 47)
#' c.i <- c(9,106,14,98,12,49,9) 
#' m1 <- meta::metabin( event.e = a.i,n.e = n.i.1,event.c = c.i,n.c = n.i.2,
#'                studlab = paste0('Study ' , 1:7) , sm = 'OR' ,
#'                comb.fixed = FALSE, comb.random = TRUE )
#' mr1 <- metarep(  m1 , u = 2, common.effect = FALSE , t = 0.05 , 
#'                alternative = 'two-sided', report.u.max = TRUE) 
#' print(mr1, digits = 2)
print.summary.metarep <- function(x, ...) {
  
  meta:::chkclass(x, "summary.metarep")
  # chkclass(x, "summary.metarep")
  
  print.summary.meta(x)
  cat(paste0("- replicability analysis (r-value = ",
             round( x$r.value, digits = 4 ) , ")\n"))
  
  if( is.null( x$u.decreased)&(!is.null( x$u.increased)) ){
    cat(paste0("- out of ", x$k, " studies, at least: ",
               x$u.increased, " with increased effect\n"))
  }
  if( is.null( x$u.increased)&(!is.null( x$u.decreased)) ){
    cat(paste0("- out of ", x$k, " studies, at least: ",
               x$u.decreased, " with decreased effect\n"))
  }
  
  if( (!is.null( x$u.increased)) & (!is.null( x$u.decreased)) ){
    cat(paste0("- out of ", x$k, " studies, at least: ",
               x$u.increased,  " with increased effect and "  ,
               x$u.decreased , " with decreased effect.\n"))
  }
  
  ##
  invisible(NULL)
}
