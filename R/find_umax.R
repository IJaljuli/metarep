find_umax <- function(x , comb.fixed = F , comb.random = T ,
                      alternative = 'two-sided' ,
                      do.truncated.umax = T , alpha.tilde = .05, confidence = 0.95 ){
  
  # function(x , comb.fixed = x$comb.fixed , comb.random = x$comb.random ,
  #          alternative = 'two-sided' ,
  #          do.truncated.umax = F , alpha.tilde = .5, confidence = 0.95 ){
    
  
  na.pvs <- all(is.na(x$pval))
  na.zvs <- (sum(!is.na(x$zval))==0 )
  if ( na.zvs & na.pvs ){
    warning('Please supply valid p-values or zvalues.')
    return( list ( u_max = NULL , worst.case = x , Side = NULL , rvalue = NULL )[c(2,1,3,4)] )
  }
  if ( comb.fixed == comb.random )  {
    comb.random <- T ; comb.fixed <- F ; do.truncated.umax <- T 
    warning('Random-effect replicability analysis is performed at 0.5.')
  }
  
  alpha = (1 - confidence)
  # do.random.u_max = Fisher's \ karl pearson . do.truncated.umax  = truncated fisher
  # Or else - alternative = 'greater' , 'less'
  chkclass(x, "meta")
  #if((!x$comb.fixed )& (! do.random.u_max)) stop('Fixed meta-analysis must be reported')
  twoSided <- (alternative == 'two-sided')
  
  nstudlab <- sum(!is.na(x$pval))
  nstudlab <- ifelse( na.pvs , sum(!is.na(x$zval)) , nstudlab )
  
  pv.greater <- ifelse(comb.fixed , x$zval.fixed , x$zval.random ) 
  
  if(!is.na(pv.greater)){ 
    pv.greater <- pnorm( pv.greater , lower.tail = F )
  }else{
    pv.greater <- ifelse(comb.fixed , x$pval.fixed , x$pval.random ) 
    TE.sign <- ifelse(comb.fixed , x$TE.fixed , x$TE.random ) 
    pv.greater <- ifelse(TE.sign > 0 , pv.greater/2 , 1-pv.greater/2) 
  }
  
  pv.less <- 1 - pv.greater
  
  rvl <- rvg <- 1
  ul <- ug <- u_max <-  0
  meta_ug <- meta_ul <- meta_ul_last_sig <- meta_ug_last_sig <- NULL
  
  # perform random replicability analysis  
  if( comb.random ){
    meta_ul_prev <-  meta_ug_prev <- NULL
    if ( alternative != 'greater' ){
      u1 <-  1  
      meta_ul <-  metaRvalue.onesided.U(x,u = u1 ,comb.fixed = F , comb.random = T,
                                      alternative = 'less',
                                      do.truncated.umax = do.truncated.umax  ,
                                      alpha.tilde = alpha.tilde )
      rvl <- meta_ul$pvalue.onesided
      while( (u1 < nstudlab ) & ( rvl <=  alpha / (1 + twoSided )) ){
        meta_ul_prev <- meta_ul
        u1 <- u1 + 1
        meta_ul <- metaRvalue.onesided.U(x,u = u1 ,comb.fixed = F , comb.random = T,
                                         alternative = 'less',
                                         do.truncated.umax = do.truncated.umax  ,
                                         alpha.tilde = alpha.tilde )
        rvl <- meta_ul$pvalue.onesided
      }
      
      ul <- u1
      
      if ( (rvl >  alpha / (1 + twoSided )) ){
        if(is.null(meta_ul_prev)){
          ul <- 0
        }else{
          rvl <- meta_ul_prev$pvalue.onesided
          meta_ul <- meta_ul_prev
          ul <- u1-1 
        }

      }
      
      u_max <- ul ;  names( u_max) <- 'u^L'
      rvalue <- rvl
      worst.case.meta <- meta_ul 
      side = 'less'
      rep.text <- paste0('out of ' , nstudlab , ' studies ', ul ,
                         ' with decreased effect.')
    }
    
    if ( alternative != 'less' ){
      u1 <-  1  
      meta_ug = metaRvalue.onesided.U(x,u = u1 ,comb.fixed = F , comb.random = T,
                                      alternative = 'greater',
                                      do.truncated.umax = do.truncated.umax  ,
                                      alpha.tilde = alpha.tilde )
      rvg <- meta_ug$pvalue.onesided
      while((u1 < nstudlab )&( rvg <=  alpha / (1 + twoSided )) ){
        meta_ug_prev <- meta_ug
        u1 <- u1 + 1
        meta_ug <- metaRvalue.onesided.U(x,u = u1 ,comb.fixed = F , comb.random = T,
                                         alternative = 'greater',
                                         do.truncated.umax = do.truncated.umax  ,
                                         alpha.tilde = alpha.tilde )
        rvg <- meta_ug$pvalue.onesided
      }
      ug <- u1
      if (rvg >  alpha / (1 + twoSided ) ){ 
        if ( is.null(meta_ug_prev) ){
          ug <- 0
          }else{
          rvg <- meta_ug_prev$pvalue.onesided
          meta_ug <- meta_ug_prev
          ug <- u1-1 
        }
      }
      u_max <- ug ;  names( u_max) <- 'u^R'
      rvalue <- rvg
      worst.case.meta <- meta_ug 
      side <- 'greater'
      rep.text <- paste0('out of ' , nstudlab , ' studies, ' , ug , ' with increased effect.')
    }
    
    if(is.null(rvl)) { rvl <- 1 ; ul <- 0 }
    if(is.null(rvg)) { rvg <- 1 ; ug <- 0 }
    
    names(rvalue) <- 'r^R'
    
    if ( alternative == 'two-sided' ){
      if( (ul > ug) | ((ul == ug)&(rvl<rvg)) ){
        u_max <- ul
        worst.case.meta <- meta_ul 
        side <- 'less'
        names(rvalue) <- 'r^L'
        names( u_max) <- 'u^L'
      }
      u_max <- c(u_max , ul , ug )
      names(u_max) <- c('u_max' , 'u^L', 'u^R')
      
      rvalue <- 2*min( c( rvl, rvg, 0.5 ))
      rvalue <- c(  rvalue , rvl , rvg )
      names( rvalue ) <- c( 'rvalue' , 'r^L' , 'r^R')
      
      rep.text <- paste0('out of ' , nstudlab , ' studies: ', ul ,
                         ' with decreased effect, and ', ug , ' with increased effect.')
      
      worst.case.meta$pvalue.onesided <- min( c( 0.5 , rvl , rvg) )*2 
    }

    names(side) <- 'Direction of the stronger signal'
    return(list(worst.case =  worst.case.meta$worst.case,
                Side = side , u_max = u_max , rvalue = rvalue ,
                Replicability_Analysis = rep.text))
  }
  
  # if reaches here, means that fixed-effect replicability analysis is performed. 
  
  
  
  ul <- 0 ; meta_ul <- NULL
  if ( alternative != 'greater' ){
    u1 <- 1 ; u2 <- nstudlab
    final_ul <- NULL
    # eleminate the radical cases: 
    meta_ul_last_sig <- meta_ul <- 
      metaRvalue.onesided.U(x,u = 1 ,comb.fixed = T , comb.random = F,
                            alternative = 'less',  do.truncated.umax = F ,
                            alpha.tilde = 1 )
    if( pv.less >  alpha / (1 + twoSided )) {
      meta_ul_last_sig <- NULL
      final_ul <- 
        list(u_max = 0 , worst.case =  meta_ul$worst.case,
             Side = 'less' , rvalue = meta_ul$pvalue.onesided )
      
    }
    
    # u1 <- u1 + 1
    
    # model with one study only, replicability at all: 
    meta_ul <-  metaRvalue.onesided.U(x,u = u2 ,comb.fixed = T , comb.random = F,
                                      alternative = 'less',
                                      do.truncated.umax = F ,
                                      alpha.tilde = 1)
    rvl <- meta_ul$pvalue.onesided
    
    if(rvl <=  alpha / (1 + twoSided )) {
      final_ul <- 
        list(u_max = u2 , worst.case =  meta_ul$worst.case, Side = 'less' , rvalue = rvl )
    }
    
    # u2 <- u2 - 1 
    # if 1 < u_max < n.studlab , search for it by devide-and-concur:
    if ( is.null(final_ul)){
      u_mid <- ceiling( (u1 + u2)/2 )
      
      while ( u_mid != u2 ){
        meta_ul <- metaRvalue.onesided.U(x,u = u_mid ,comb.fixed = T , comb.random = F,
                                         alternative = 'less',
                                         do.truncated.umax = F ,
                                         alpha.tilde = 1)
        
        if ( meta_ul$pvalue.onesided < alpha / (1 + twoSided ) ){
          u1 <- u_mid 
          meta_ul_last_sig <- meta_ul
        }else{
          u2 <- u_mid 
        }
        u_mid <- ceiling( (u1 + u2)/2 )
      }
      ul <- u1
      meta_ul <- meta_ul_last_sig 
      rvl <- meta_ul$pvalue.onesided
      side <- 'less'
      
      final_ul <- 
        list(u_max = ul , worst.case =  meta_ul$worst.case, Side = 'less' , rvalue = rvl )
    }
    
    
    if( alternative == 'less'){
      return(final_ul[c(2,1,3,4)])
    }
    
  }  
  
  ug <- 0 ; meta_ug <- NULL
  if ( alternative != 'less' ){
    u1 <- 1 ; u2 <- nstudlab
    final_ug <- NULL
    
    # eleminate the radical cases: 
    meta_ug <- meta_ug_last_sig <-
      metaRvalue.onesided.U(x,u = 1 ,comb.fixed = T , comb.random = F,
                            alternative = 'greater', do.truncated.umax = F ,
                            alpha.tilde = 1 )
    
    if( pv.greater >  alpha / (1 + twoSided )) {
      meta_ug_last_sig <- NULL
      final_ug <-
        list(u_max = 0 , worst.case =  meta_ug$worst.case,
             Side = 'greater' , rvalue = meta_ug$pvalue.onesided )
    }
    
    # the original model
    
    # u1 <- u1 + 1
    
    # model with one study only, replicability at all: 
    meta_ug = metaRvalue.onesided.U(x,u = u2 ,comb.fixed = T , comb.random = F,
                                    alternative = 'greater',
                                    do.truncated.umax = F ,
                                    alpha.tilde = 1 )
    rvg <- meta_ug$pvalue.onesided
    
    if(rvg <=  alpha / (1 + twoSided )) {
      final_ug <-
        list(u_max = u2 , worst.case =  meta_ug$worst.case, Side = 'greater' , rvalue = rvg )
      
    }
    
    # u2 <- u2 - 1 
    
    # if 1 < u_max < n.studlab , search for it by devide-and-concur:
    if (is.null(final_ug)) {
      u_mid <- ceiling( (u1 + u2)/2 )
      
      while ( u_mid != u2 ){
        meta_ug <- metaRvalue.onesided.U(x,u = u_mid ,comb.fixed = T , comb.random = F,
                                         alternative = 'greater',
                                         do.truncated.umax = F ,
                                         alpha.tilde = 1 )
        
        if ( meta_ug$pvalue.onesided < alpha / (1 + twoSided ) ){
          u1 <- u_mid 
          meta_ug_last_sig <- meta_ug 
        }else{
          u2 <- u_mid 
        }
        u_mid <- ceiling( (u1 + u2)/2 )
      }
      ug <- u1
      meta_ug <- meta_ug_last_sig 
      rvg <- meta_ug$pvalue.onesided
      final_ug <-
        list(u_max = ug , worst.case =  meta_ug$worst.case, Side = 'greater' , rvalue = rvg )
    }
    final <- final_ug
    side = 'greater'
    if( alternative == 'greater'){
      return(final_ug[c(2,1,3,4)])
    }
    
  }
  
  if(is.null(rvl)) { rvl <- 1 ; ul <- 0 }
  if(is.null(rvg)) { rvg <- 1 ; ug <- 0 }
  
  if ( alternative == 'two-sided' ){
    if( (ul > ug) | ((ul == ug)&(rvl<rvg)) ){
      final <- final_ul
      side <- 'less'
    }
    final$rvalue <- min( min(rvl,rvg)*2 ,1 )
    final <- final[c(2,1,3,4)]
    return(final)
  }
  
}