forest.metarepl <- function(x, ...) {
  
  # meta:::chkclass(x, "metarepl")
  chkclass(x, "metarepl")
  
  u_max_text <- NULL
  
  rvalue.text <- paste0("Replicability analysis (r-value = ",
                        round(x$r.value,digits = 4) , ")")
  
  if( is.null( x$u_L) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_R, " with increased effect")
  }
  if( is.null( x$u_R) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_L, " with decreased effect")
  }
  
  if( (!is.null( x$u_R)) & (!is.null( x$u_L)) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_R,  " with increased effect and "  ,
               x$u_L , " with decreased effect.")
  }
  
  
  forest.meta(x,
              text.addline1 = rvalue.text ,
              text.addline2 = u_max_text ,
              ...)
  ##
  invisible(NULL)
}
