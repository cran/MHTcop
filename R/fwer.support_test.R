#' Copula-based multiple support test which controlls the FWER
#'
#' Perform a multiple support test controlling the family-wise error rate (FWER) using the procedure described in Stange, Bodnar, Dickhaus (2015).
#'
#' The test is performed assuming an i.i.d. sample \eqn{X_1,\cdots,X_n} which has the stochastic representation
#' \deqn{X_{i,j}=\vartheta_j Z_j} where \eqn{Z_j} takes values in \eqn{[0,1]} and which is distributed according to a Gumbel copula with Beta margins. The test simultaneously tests the hypotheses \eqn{H_{0,j}: \vartheta_j \le \vartheta_j^*} versus the corresponding alternatives \eqn{H_{1,j}: \vartheta_j>\vartheta_j^*}.
#'
#' For usage examples and figure reproduction see \code{vignette('fwer-support-test',package='MHTcop')}.
#'
#' Note: If the copula is only in the domain of attraction of the Gumbel copula (but not a Gumbel copula) then it is necessary to pass the
#' number of boot strap repetitions `boot.reps` as an additional parameter since the non-bootstrapped parameter estimate would not be consistent.
#'
#' @param sample The observed sample (a matrix whose columsn are the observations)
#' @param theta The hypothesized scale \code{theta=c}\eqn{(\vartheta_1^*,\cdots,\vartheta_m^*}\code{)}
#' @param alpha First shape parameter of the Beta margins
#' @param beta Second shape parameter of the Beta margins
#' @param boot.reps number of bootstrap repetitions for estimating the parameter \eqn{\eta} of the Gumbel copula.
#'                  If this parameter is NULL then \eqn{\eta} is estimated from Kendalls tau and no bootstrap is performed.
#' @param sigLevel The desired significance level
#' @references J. Stange, T. Bodnar and T. Dickhaus (2015). Uncertainty quantification for the family-wise error rate in multivariate copula models. \emph{AStA Advances in Statistical Analysis} 99.3 (2015): 281-310.
#' @export
#' @return list l, where\itemize{
#'                 \item{l$statistic contains the values of the test statistics,}
#'                 \item{l$critvalues are the calibrated critical values,}
#'                 \item{l$test contains the test decisions,}
#'                 \item{l$etahat is estimated parameter of the Gumbel copula}}
fwer.support_test <- function(sample,theta,alpha=3,beta=4,boot.reps=NULL,sigLevel=0.05) {
  statistic <- function (X,theta) {
    T<-apply(X,2,max)/theta
    return(T);
  }

  estimate <- function(X)
  {
    Tau<-stats::cor(X,method="kendall")
    tauhat<-mean(Tau[upper.tri(Tau)])
    etahat<-1/(1-tauhat)
    etahat
  }

  estimate_bootstrap <- function(X, num.reps = 400)
  {
    num.reps <- as.numeric(num.reps)
    n<-nrow(X)
    m<-ncol(X)
    k<-sqrt(n)
    bootstrap_statistic<-replicate(num.reps,apply(X[sample(n,k,replace=TRUE),],2,max)) #returns m rows and B columns
    #each column is a resampled statistic
    #scale invariance of kendall's tau does not require further manipulations
    #but observations have to be transposed for call to cor()
    Tau<-stats::cor(t(bootstrap_statistic),method="kendall")
    tauhat<-mean(Tau[upper.tri(Tau)])
    etahat<-1/(1-tauhat)
    etahat
  }

  #helper function to invert from global to (equicoordinate) local significance level for given Gumbel copula parameter eta
  Invert_Gumbel<-function(y,m,eta) exp(m^(-1/eta)*log(y))

  critical_values <- function(n,theta,alpha,beta,eta,sigLevel=0.05) {
    m <- length(theta)
    cval<-list()
    cval$Bonferroni <- stats::qbeta((1-sigLevel/m)^(1/n),alpha,beta)
    cval$Sidak <- stats::qbeta((1-sigLevel)^(1/(m*n)),alpha,beta)
    UU<-Invert_Gumbel(1-sigLevel,m,eta)
    cval$Empirical <- stats::qbeta(UU^(1/n),alpha,beta)
    return(cval)
  }

  stat            <- statistic(sample,theta)
  if(!is.null(boot.reps)) {
    etahat        <- estimate_bootstrap(sample,boot.reps)
  } else {
    etahat        <- estimate(sample)
  }
  cval            <- critical_values(nrow(sample),theta,alpha,beta,etahat,sigLevel)
  t               <- lapply(cval,function(c) stat>c)

  result <- list(statistic=stat,critvalues=cval,test=t,etahat=etahat)
}
