metaRvalue.onesided.U <- function (x,u = 2 , comb.fixed = F , comb.random = T ,
                                   alternative = 'less',
                                   do.truncated.umax = T ,
                                   alpha.tilde = .05 ){
  # function (x,u = 2 , comb.fixed = x$comb.fixed , comb.random = x$comb.random ,
  #           alternative = 'less',
  #           do.truncated.umax = F ,
  #           alpha.tilde = NULL ){
    # Returns the worst case meta model. 
  chkclass(x, "meta")
  metaInf <- inherits(x,'metainf')
  x.original = x
  if (do.truncated.umax & (alpha.tilde == 1) ) {
    message(' truncation at 0.05 is performed')
    alpha.tilde <- .05
  }
  if( is.null(x$zval) | all( is.na(x$zval) ) ){
    stop('Please supply model with z-values')
  }
  # x = update.meta(object = x , subset = (1:length(x$studlab)) [ !is.na(x$pval)] )
  x.fixed.metainf <- NULL
  x.fixed.w <- NULL
  if(metaInf){ x <- extract.meta.modle (x) }
  nstudlab <- length(x$studlab)
  if (( u > nstudlab ) | (u < 1) ){
    stop ( 'invalid tuning parameter u ' )
  }
  
  
  if ( (! comb.fixed )&( ! comb.random ) ) stop('Desired replicability model must be supplied')
  if(!(alternative %in% c('less','greater'))) stop('Supply informative alternative ( "less" or "greater" )')

  
  
  
    
  ## fixed-effects replicability analysis. 
  worst.case.fixed <- studies_subsets <- NULL
  if ( comb.fixed ){
    
    if( u == 1 ){
      return(list(worst.case = x ,
                  Side = alternative,
                  pvalue.onesided = pnorm(x$zval.fixed , lower.tail = ( alternative == 'less'))  ))
    }
    
    studies_subsets <- combn(x = 1:nstudlab , m = nstudlab-u+1 )
    studies_subsets <- rbind(studies_subsets , rep(NA , ncol(studies_subsets)))
    rownames(studies_subsets) <- c( paste0('s' , 1:(nstudlab-u+1)) , 'rvalue' )
    for ( k in 1:ncol(studies_subsets)){
      worst.studies.fisher = x$studlab[ studies_subsets[-c(nrow(studies_subsets)) ,k] ]
      worst.studies.fisher = which( x$data$.studlab %in% worst.studies.fisher)
      x.sub <- update.meta(x,subset =  worst.studies.fisher )
      z.val <- x.sub$zval.fixed
      studies_subsets['rvalue' , k] <-  x.sub$pval.fixed / 2
      if ( ((z.val < 0) & (alternative == 'greater')) | ( (z.val > 0) & (alternative == 'less') ) ){
        studies_subsets['rvalue' , k] = 1 - studies_subsets['rvalue' , k]  
      }
    }
    # derive the line of worst case study 
    k = max(which.max(studies_subsets['rvalue',])) # studies column pointer. 
  
    worst.case.studies = x$studlab[ studies_subsets[-c(nrow(studies_subsets)) ,k] ]
    worst.case.studies = which( x$data$.studlab %in% worst.case.studies )
    
    worst.case <- update.meta(x ,subset = worst.case.studies )
    return(list(worst.case = worst.case ,
                Side = alternative,
                pvalue.onesided = studies_subsets['rvalue' , k] ))
  }
  # so far so good 
  
  ## random-effects replicability analysis:
  
  zval.all <- x$zval[!is.na(x$zval)]
  pvs.all <- pnorm( zval.all ,lower.tail = (alternative == 'less'))
  
  if( u == 1 ){
    pvo <- truncatedPearson( pvs.all ,alpha.tilde =  alpha.tilde)
    return(list(worst.case = x ,
                Side = alternative,
                pvalue.onesided =  pvo$rvalue ))
  }
  if( u == length(pvs.all) ){
    pvo <- max(pvs.all)
    if ( pvo < alpha.tilde){
      return(list(worst.case = NULL ,
                  Side = alternative,
                  pvalue.onesided =  pvo ))
    }else{
      return(list(worst.case = NULL ,
                  Side = alternative,
                  pvalue.onesided =  1 ))
    }
    
  }
  
    if(alternative == 'greater'){
      wsf <- which( (nstudlab+1-rank(zval.all)) >= u )
    }else{
      wsf <- which( rank(zval.all) >= u )
    }
    
  # }else{
  #   pvs.all <- x$pval/2 ; pvs.all <- pvs.all[!is.na(pvs.all)]
  #   if(alternative == 'greater'){
  #     pvs.all[zval.all<0] = 1 - pvs.all[zval.all<0]
  #   }else{
  #     pvs.all[zval.all>0] = 1 - pvs.all[zval.all>0]
  #   }
  # }
  # wsf <- which( rank(pvs.all) >= u )
  
  worst.studies.fisher <- x$studlab[ wsf ]
  worst.studies.fisher <- which( x$data$.studlab %in% worst.studies.fisher)
  worst.case <- update.meta(x ,subset = worst.studies.fisher )
  worst.pvs.fisher <- truncatedPearson( pvs.all[ wsf ] ,alpha.tilde =  alpha.tilde)
  
  return(list(worst.case = worst.case ,
              Side = alternative,
              pvalue.onesided = worst.pvs.fisher$rvalue ))
}
