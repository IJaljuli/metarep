metarepl <- function(x, u = 2, common.effect = FALSE , t = NULL , alternative = 'two-sided') {
  
  meta:::chkclass(x, "meta")
  performing.truncated <- F
  
  if( is.null(t) ){
    if (!common.effect){
      stop("Error: Must specify truncation threshold t. 
         For replicability analysis with common-effect assumption, set common.effect = TRUE ")
    }
  }else{
    if ( t == 1  ) {
      common.effect <- FALSE
      message( "Performing Replicability analysis via original Pearson's test" )
    }
    
    if ( common.effect & ( t < 1) ) {
      common.effect <- FALSE
      message( paste0("Performing Replicability analysis via truncated Pearson's test, 
                      with truncation threshold t = ", t ) )
    }
    
  }
  
  ##
  ## Do replicability analysis
  ##
  
  rvalue <- rvalue.less <- rvalue.greater <- NULL 
  if ( alternative != 'two-sided' ){
    rvalue.results <- metaRvalue.onesided.U( x = x , u = u , alternative = alternative  , 
                                             comb.fixed = common.effect,
                                             comb.random = !common.effect,
                                             do.truncated.umax = !is.null(t), 
                                             alpha.tilde = ifelse(is.null(t) , 1 , t ) )
    rvalue <- rvalue.results$pvalue.onesided
    side <- alternative
  }else{
    rvalue.results.less <-metaRvalue.onesided.U( x = x , u = u , alternative = 'less', 
                                                 comb.fixed = common.effect,
                                                 comb.random = !common.effect,
                                                 do.truncated.umax = !is.null(t), 
                                                 alpha.tilde = ifelse(is.null(t) , 1 , t ) )
    rvalue.less <- rvalue.results.less$pvalue.onesided
    
    rvalue.results.greater <- metaRvalue.onesided.U( x = x , u = u , alternative = 'greater', 
                                                     comb.fixed = common.effect,
                                                     comb.random = !common.effect,
                                                     do.truncated.umax = !is.null(t), 
                                                     alpha.tilde = ifelse(is.null(t) , 1 , t ) )
    
    rvalue.greater <- rvalue.results.greater$pvalue.onesided
    
    if ( rvalue.less < rvalue.greater ){
      rvalue.results <- rvalue.results.less
      side <- 'less'
    }else{
      rvalue.results <- rvalue.results.greater
      side <- 'greater'
    }
    
    rvalue = min(1 , 2*rvalue.results$pvalue.onesided )
  }
  
  ## 
  ## Replicability analysis results
  ##
  res <- x
  res$r.value <- rvalue
  res$Side <- side
  res$worst.case.studies <- (rvalue.results$worst.case)$studlab
  
  Umax_right  <- Umax_left <- Umax <- Side <-  NULL

    if( common.effect ){
      
      Umax = find_umax( x , alternative = alternative, confidence = 0.95 ,do.truncated.umax = F ,
                        comb.fixed = T , comb.random = F )
      
      res$u <-  Umax$u_max 
      res$side <-  Side
      
    }else{
      
      if( alternative != 'less' ){
        Umax_right = find_umax( x , alternative = 'greater', confidence = 0.975 ,do.truncated.umax = !is.null(t),
                                comb.fixed = common.effect , comb.random = !common.effect , 
                                alpha.tilde = ifelse(is.null(t) , 1 , t ) )
        res$u_R <- Umax_right$u_max
      }
      
      if( alternative !=  'greater' ){
        Umax_right = find_umax( x , alternative = 'less', confidence = 0.975 ,do.truncated.umax = !is.null(t),
                                comb.fixed = common.effect , comb.random = !common.effect , 
                                alpha.tilde = ifelse(is.null(t) , 1 , t ) )
      
        res$u_L <- Umax_left$u_max
      }
      
    }
  
 
  
  ##
  class(res) <- c("metarepl", class(res))
  res
  }


summary.metarepl <- function(object, ...) {
  
  meta:::chkclass(object, "metarepl")
  
  res <- summary.meta(object)
  ##
  res$r.value <- object$r.value
  res$n.increased <- object$n.increased
  ##
  class(res) <- c("summary.metarepl", class(res))
  ##
  res
}


print.summary.metarepl <- function(x, ...) {
  
  meta:::chkclass(x, "summary.metarepl")
  
  print.summary.meta(x)
  cat(paste0("- replicability analysis (r-value = ",
             x$r.value, ")\n"))
  cat(paste0("- out of ", x$k, " studies, at least: ",
             x$n.increased, " with increased effect\n"))
  ##
  invisible(NULL)
}


forest.metarepl <- function(x, ...) {
  
  meta:::chkclass(x, "metarepl")
  
  forest.meta(x,
              text.addline1 = paste0("Replicability analysis (r-value = ",
                                     x$r.value, ")"),
              text.addline2 = paste0("Out of ", x$k, " studies, at least: ",
                                     x$n.increased, " with increased effect"),
              ...)
  ##
  invisible(NULL)
}


# library(meta)
# m1 <- metagen(1:5, rep(0.25, 5))
# mr1 <- metarepl(m1)
# 
# m1
# mr1
# summary(mr1)
# forest(mr1)
# 
# class(m1)
# class(mr1)
# class(summary(mr1))
