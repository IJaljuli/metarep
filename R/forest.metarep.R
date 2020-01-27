forest.metarep <- function(x, ...) {
  
  meta:::chkclass(x, "metarep")
  # chkclass(x, "metarep")
  
  u_max_text <- NULL
  
    rvalue.text <- meta:::formatPT(x$r.value , digits = 4)
    if( rvalue.text == '1.0000'  ){ 
      rvalue.text  <- '1' }
    
    rvalue.text <- paste0("Replicability analysis (r-value = ",
                          rvalue.text, ")")
  
  if( is.null( x$u_L)&(!is.null( x$u_R)) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_R, " with increased effect")
  }
  if( is.null( x$u_R)&(!is.null( x$u_L)) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_L, " with decreased effect")
  }
  
  if( (!is.null( x$u_R)) & (!is.null( x$u_L)) ){
    u_max_text <-paste0("Out of ", x$k, " studies, at least: ",
               x$u_R,  " with increased effect and "  ,
               x$u_L , " with decreased effect.")
  }
  
  
  if(!is.null(u_max_text)){
    forest.meta(x,
                text.addline1 = rvalue.text ,
                text.addline2 = u_max_text ,
                ...)
  }else{
    forest.meta(x,
                text.addline1 = rvalue.text ,
                ...)
  }
  ##
  invisible(NULL)
}
