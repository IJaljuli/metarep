print.summary.metarep <- function(x, ...) {
  
  meta:::chkclass(x, "summary.metarep")
  # chkclass(x, "summary.metarep")
  
  print.summary.meta(x)
  cat(paste0("- replicability analysis (r-value = ",
             x$r.value, ")\n"))
  
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
