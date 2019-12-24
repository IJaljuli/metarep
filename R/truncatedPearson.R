truncatedPearson <- function( p , alpha.tilde = 1 ){
  p <- p[!is.na(p)]
  if( alpha.tilde < 1){
    L = length(p)
    w = prod( p[ p <= alpha.tilde ] )
    db = dbinom(x = 1:L , size = L ,prob = alpha.tilde)
    pg = pgamma( -log( w / (alpha.tilde^(1:L) ) ) , shape = 1:L ,lower.tail = F  )
    
    TP.pvalue <- ifelse( sum( p <= alpha.tilde ) >0 ,  sum( db * pg) , 1 )
    
    return( list( chisq = NULL , df = NULL , rvalue = TP.pvalue , validp = p ) )
  }
  
  if ( length(p)<=1 ){
    stop( 'Error: Meta-analysis must include at least 2 p-values' )
  }
  
  p [  which( p < 10^-40 ) ] <- 10^-39
  
  output <- metap::sumlog(p) 
  return( list( chisq = output$chisq , df = output$df ,
                rvalue = output$p , validp = output$validp ) )
  
}