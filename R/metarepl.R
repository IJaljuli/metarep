metarepl <- function(x, u = 2, common.effect = FALSE , t = NULL , alternative = 'two-sided') {
  
  # meta:::chkclass(x, "meta")
  chkclass(x, "meta")
  performing.truncated <- F
  
  if( is.null(t) ){
    if (!common.effect){
      stop("Error: Must specify truncation threshold t <= 1 .
           For replicability-analysis with common-effect assumption, set common.effect = TRUE ")
    }
    }else{
      if ( t == 1  ) {
        if(common.effect){ t <- NULL }else{
          message( "Performing Replicability analysis via original-Pearson's test" )
        }
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
  
  # compute r-value
  rvalue <- rvalue.less <- rvalue.greater <- NULL 
  if ( alternative != 'two-sided' ){
    rvalue.results <- metaRvalue.onesided.U( x = x , u = u , alternative = alternative  , 
                                             comb.fixed = common.effect,
                                             comb.random = !common.effect,
                                             do.truncated.umax = ifelse(is.null(t) , F , t < 1 ), 
                                             alpha.tilde = ifelse(is.null(t) , 1 , t ) )
    rvalue <- rvalue.results$pvalue.onesided
    side <- alternative
  }else{
    rvalue.results.less <- metaRvalue.onesided.U( x = x , u = u , alternative = 'less', 
                                                  comb.fixed = common.effect,
                                                  comb.random = !common.effect,
                                                  do.truncated.umax = ifelse(is.null(t) , F , t < 1 ), 
                                                  alpha.tilde = ifelse(is.null(t) , 1 , t ) )
    rvalue.less <- rvalue.results.less$pvalue.onesided
    
    rvalue.results.greater <- metaRvalue.onesided.U( x = x , u = u , alternative = 'greater', 
                                                     comb.fixed = common.effect,
                                                     comb.random = !common.effect,
                                                     do.truncated.umax = ifelse(is.null(t) , F , t < 1 ), 
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
  
  # find u_max
  Umax_right  <- Umax_left <- Umax <- Side <-  NULL
  
  if( common.effect ){
    
    Umax = find_umax( x , alternative = alternative, confidence = 0.95 ,
                      do.truncated.umax = F ,
                      comb.fixed = T , comb.random = F )
    
    if( Umax$side == 'less '){
      res$u_L <-  Umax$u_max 
    }else{
      res$u_R <-  Umax$u_max 
    }
    
    
  }else{
    
    if( alternative != 'less' ){
      Umax_right = find_umax( x , alternative = 'greater',
                              confidence = 1 - 0.05/(1+alternative == 'two-sided') ,
                              do.truncated.umax = ifelse(is.null(t) , F , t < 1 ), 
                              comb.fixed = common.effect , comb.random = !common.effect , 
                              alpha.tilde = ifelse(is.null(t) , 1 , t ) )
      res$u_R <- Umax_right$u_max
    }
    
    if( alternative !=  'greater' ){
      Umax_right = find_umax( x , alternative = 'less', 
                              confidence = 1 - 0.05/(1+alternative == 'two-sided') ,
                              do.truncated.umax = ifelse(is.null(t) , F , t < 1 ), 
                              comb.fixed = common.effect , comb.random = !common.effect , 
                              alpha.tilde = ifelse(is.null(t) , 1 , t ) )
      
      res$u_L <- Umax_left$u_max
    }
    
  }
  
  
  ## 
  ## Replicability analysis results
  ##
  res <- x
  res$r.value <- rvalue
  res$side <- side
  res$worst.case.studies <- (rvalue.results$worst.case)$studlab
  if( common.effect ){
    res$repl.method <- 'Common-effect' 
  }else{
    res$repl.method <- ifelse( t == 1,  'Truncated Pearson' , 'Pearson')
  }
  
  
  ##
  class(res) <- c("metarepl", class(res))
  res
  
}
