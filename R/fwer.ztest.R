#' Copula-based multiple z-test which controlls the FWER
#'
#' Perform a multiple (two-sided) z-test controlling the family-wise error rate (FWER) using the procedure described in Stange, Bodnar, Dickhaus (2015).
#'
#' Let \eqn{X_1,\cdots,X_n} denote an i.i.d. sample with values in \eqn{{\rm I\!R}^m}. Furthermore let \eqn{\mu_j={\rm I\!E}[X_{1,j}]} be the component-wise
#' expectations. Then the multiple (two-sided) z-test simultaneously tests the hypotheses \eqn{H_{0,j}: \mu_j = \mu_j^*} versus the corresponding alternatives \eqn{H_{1,j}: \mu_j\not=\mu_j^*}.
#'
#' For usage examples and figure reproduction see \code{vignette('fwer-ztest',package='MHTcop')}.
#'
#' Note: If the parameter \code{sigma} is passed it needs to be a consistent estimate of the covariance matrix of \eqn{X_1}.
#'
#' @param sample The observed sample
#' @param mu The mean \eqn{\mu^*}
#' @param sigma The estimated covariance matrix (the copula parameter). If it is omitted it will be estimated from an AR(1) model
#' @param sigLevel The desired significance level
#' @references J. Stange, T. Bodnar and T. Dickhaus (2015). Uncertainty quantification for the family-wise error rate in multivariate copula models. \emph{AStA Advances in Statistical Analysis} 99.3 (2015): 281-310.
#' @export
#' @return list l, where\itemize{
#'                 \item{l$statistic contains the values of the test statistics,}
#'                 \item{l$critvalues are the calibrated critical values,}
#'                 \item{l$test contains the test decisions,}
#'                 \item{l$etahat is estimated parameter of the Gumbel copula}}
fwer.ztest <- function(sample,mu,sigma=NULL,sigLevel=0.05) {
  statistic <- function (X, mu) {
    M<-colMeans(X)
    n<-nrow(X)
    T<-sqrt(n)*abs(M-mu)
    return(T)
  }

  estimate <- function(X)
  {
    helper1<-function(x) {y<-x[-1]; z<-x[-length(x)]; return(sum(y*z))}
    helper2<-function(x,I=(1:length(x))) sum((x[I])^2)
    n<-nrow(X)
    m<-ncol(X)
    N<-n*(m-1)
    M<-colMeans(X)
    CX<-t(t(X)-M)
    SP<-sum(apply(CX,1,helper1))/N
    SQ1<-sum(apply(CX,1,function(x) helper2(x)))/N
    SQ2<-sum(apply(CX,1,function(x) helper2(x,I=(2:(m-1)))))/N
    solutions<-polyroot(c(-SP,SQ1+SQ2-1,-SP,1))
    rhohat<-Re(solutions[1])
  }

  Sigma_GaussianAR1<-function(m,rho) {
    exponents<-abs(outer(1:m,1:m,"-")) #exponents for AR(1) matrix
    Sigma<-rho^exponents               #computes the AR(1) matrix
    return(Sigma)
  }

  critical_values <- function(mu,sigma,sigLevel=0.05) {
    m <- length(mu)
    cval<-list()
    cval$Bonferroni <- stats::qnorm(1-sigLevel/(2*m))
    alphasidak<-1-(1-sigLevel)^(1/m)
    cval$Sidak <- stats::qnorm(1-alphasidak/2)
    cval$Empirical <- mvtnorm::qmvnorm(1-sigLevel,interval=c(stats::qnorm(1-sigLevel),cval$Sidak),mean=mu,sigma=sigma,tail="both.tails")$quantile
    return(cval)
  }

  stat            <- statistic(sample,mu)
  rhohat          <- estimate(sample)
  if(is.null(sigma)) {
    sigmahat      <- Sigma_GaussianAR1(length(mu),rhohat)
  } else {
    sigmahat      <- sigma
  }
  cval            <- critical_values(mu,sigmahat,sigLevel)
  t               <- lapply(cval,function(c) stat>c)

  result <- list(statistic=stat,critvalues=cval,test=t,sigmahat=sigmahat)
}
