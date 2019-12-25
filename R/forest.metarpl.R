forest.metarepl <- function(x, ...) {
  
  # meta:::chkclass(x, "metarepl")
  chkclass(x, "metarepl")
  
  u_max_text <- NULL
  
  if( is.null( x$u.decreased) ){
    u_max_text <-paste0("- out of ", x$k, " studies, at least: ",
                        x$u.increased, " with increased effect\n")
  }
  if( is.null( x$u.increased) ){
    u_max_text <-paste0("- out of ", x$k, " studies, at least: ",
                        x$u.decreased, " with decreased effect\n")
  }
  
  if( (!is.null( x$u.increased)) & (!is.null( x$u.decreased)) ){
    u_max_text <- paste0("- out of ", x$k, " studies, at least: ",
                         x$u.increased,  " with increased effect and "  ,
                         x$u.decreased , " with decreased effect.\n")
  }
  
  forest.meta(x,
              text.addline1 = paste0("Replicability analysis (r-value = ",
                                     x$r.value, ")"),
              text.addline2 = u_max_text ,
              ...)
  ##
  invisible(NULL)
}
