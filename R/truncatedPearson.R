#' Truncated-Pearsons' test 
#' @description Apply Truncated-Pearsons' test or ordinary Pearsons' test on one-sided p-values.
#' @param p one-sided p-values of the individual studies for testing one-sided alternative based on z-test. 
#' @param alpha.tilde truncartion threshold for truncated-Pearson test. Use alpha.tilde = 1 for ordinary Pearsons' test for combining p-values.
#'
#' @return 
#' \itemize{
#'   \item{chisq: }{Pearson test statistic }
#'   \item{df: }{degrees of freedom of truncated-Pearson statistic }
#'   \item{rvalue: }{p-value of the test }
#'   \item{validp: }{p-values used in the test. }
#' }
#' 
#' @export
#'
#' @examples 
#' truncatedPearson( p = c( 0.001 , 0.01 , 0.1  ) , alpha.tilde = 1 )
#' truncatedPearson( p = c( 0.001 , 0.01 , 0.1  ) , alpha.tilde = 0.05 )
truncatedPearson <- function( p , alpha.tilde = 1 ){
  if ( length(p)<=1 ){
    stop( 'Error: Meta-analysis must include at least 2 p-values' )
  }
  p <- p[!is.na(p)]
  if ( length(p)<=1 ){
    stop( 'Error: Meta-analysis must include at least 2 p-values' )
  }
  
  if( alpha.tilde < 1){
    L = length(p)
    w = prod( p[ p <= alpha.tilde ] )
    db = dbinom(x = 1:L , size = L ,prob = alpha.tilde)
    pg = pgamma( -log( w / (alpha.tilde^(1:L) ) ) , shape = 1:L ,lower.tail = F  )
    
    TP.pvalue <- ifelse( sum( p <= alpha.tilde ) >0 ,  sum( db * pg) , 1 )
    
    return( list( p.value = TP.pvalue , validp = p ) )
  }
  
  p [  which( p < 10^-40 ) ] <- 10^-39
  
  output <- metap::sumlog(p) 
  return( list( p.value = output$p , validp = output$validp,
                chisq = output$chisq , df = output$df ) )
  
}