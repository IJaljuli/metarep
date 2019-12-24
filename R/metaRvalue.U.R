metaRvalue.U <- function (x, u = 2 ,comb.fixed = F , comb.random = T , 
                          alternative = 'two-sided', do.truncated.umax = T , alpha.tilde = 0.05 ){
  # function (x, u = 2 ,comb.fixed = x$comb.fixed , comb.random = x$comb.random, 
  #           alternative = 'two-sided', do.truncated.umax = F , alpha.tilde = 0.5 ){
  if (do.truncated.umax & (alpha.tilde == 1) ) {
    message(' truncation at 0.5 is performed')
    alpha.tilde <- .5
  }
  
  if ( comb.fixed == comb.random )  {
    comb.random <- T ; comb.fixed <- F ; do.truncated.umax <- F 
    warning('Random-effect replicability analysis is performed, without truncation.')
  }
  
  if ( (!comb.random) & (do.truncated.umax) )  {
    comb.random <- T ; comb.fixed <- F ; do.truncated.umax <- T ; alpha.tilde <- 0.5
    warning('Random-effect replicability analysis is performed, with truncation at 0.5.')
  }
  
  rvalue <- rvalue.less <- rvalue.greater <- NULL 
    if (alternative != 'two-sided' ){
      rvalue.results <- metaRvalue.onesided.U( x , u , alternative = alternative  , 
                                             comb.fixed = comb.fixed ,comb.random = comb.random ,
                                             do.truncated.umax = do.truncated.umax , 
                                             alpha.tilde = alpha.tilde)
      rvalue <- rvalue.results$pvalue.onesided
      side <- alternative
    }else{
      rvalue.results.less <- metaRvalue.onesided.U( x , u , alternative = 'less'  , 
                                                    comb.fixed = comb.fixed ,comb.random = comb.random ,
                                                    do.truncated.umax = do.truncated.umax , 
                                                    alpha.tilde = alpha.tilde)
      rvalue.less <- rvalue.results.less$pvalue.onesided

      rvalue.results.greater <- metaRvalue.onesided.U( x , u , alternative = 'greater'  , 
                                                       comb.fixed = comb.fixed ,comb.random = comb.random ,
                                                       do.truncated.umax = do.truncated.umax, 
                                                       alpha.tilde = alpha.tilde )
      
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
    return(list (worst.case = rvalue.results$worst.case ,
                 rvalue = rvalue , Side = side) )

}