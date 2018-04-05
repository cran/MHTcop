#' Approximate adjustment factor for "False Discovery Rate Control under Archimedean Copula"
#'
#' Approximates the adjustement factor described in Bodnar, Dickhaus (2014) using a Monte Carlo integration.
#'
#' @param cop The copula
#' @param m The number of hypotheses
#' @param m0 A reasonable lower bound for the number of true hypotheses
#' @param alpha the desired significance level
#' @param num.reps The number of Monte Carlo simulations to use
#' @param calc.var Additionaly calculate the variance (mostly for the examples)
#'
#' @export
#' @keywords internal
#' @references T. Bodnar and T. Dickhaus (2014). False discovery rate control under Archimedean copula. \emph{Electronic Journal of Statistics} Volume 8, Number 2 (2014), 2207-2241.

ac_fdr.calc_delta <- function(cop,m,m0,alpha=0.05,num.reps=1e5,calc.var=FALSE) {
  if(m >= 50) warning("ac_fdr.calc_delta: m is large and therefore the result is likely to be inaccurate (and possibly even invalid)!");
  if(m0>m) stop( "ac_fdr.calc_delta: m0 cannot be larger than m (at most all hypotheses are true)")
  if(m0 < 1) stop(m0>=1,"ac_fdr.calc_delta: m0 needs to be larger than zero!")
  qk<-seq(1,m,1)*alpha/m  # critical values

  theta <- cop@theta
  iPsi <- function(u) cop@iPsi(u,theta)
  iPsi.qk <- iPsi(qk)

  #calculate delta
  zst_kf<-function(s){log(1+1/s)/(iPsi.qk[s]-iPsi.qk[s+1])}
  zst_k<-sapply(1:(m0-1),zst_kf) # calculation of z_k^*
  zst_k.max <- max(zst_k)
  Gk_f<-function(s){
    a<-exp(-zst_k[s]*iPsi.qk[(m-m0+1):m])
    a[1:s]<-0
    return(bolshev.rec.vec(cbind(a[m0:2]))[1])
    #return(bolshev.rec(a[2:m0]))
  }
  Gk<-sapply(1:(m0-1),Gk_f)
  sampleAndCalc <- function(n) {
    Z<-sample.Z(cop,n)
    Z<-Z[Z>=zst_k.max]
    a <- g_z <- matrix(0,nrow=m0,ncol=length(Z))
    for(s in 1:m0) {
      a[m0-s+1,] <- exp(-Z*iPsi.qk[s+m-m0])
      g_z[s,] <- a[m0-s+1,]/qk[s+m-m0]
    }
    GkZ.diff <- bolshev.rec.vec(a[1:(m0-1),])-Gk
    g_z.diff <- g_z[2:m0,]-g_z[1:(m0-1),]
    s <- sum(g_z.diff*GkZ.diff/num.reps)
    if(!calc.var) s
    else c(s,sum(g_z.diff*GkZ.diff^2/num.reps))
  }
  n <- 1e4
  res <- replicate(ceiling(num.reps/n),sampleAndCalc(n))
  if(calc.var) {
    res <- rowSums(res)
    c(1-res[1],res[2]-res[1]^2)
  }
  else 1 - sum(res)
}
