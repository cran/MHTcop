#' Perform a FDR controlling test on marginal p-values that are distributed according to an Archmidean copula
#'
#' Performs a test on marginal p-values according to the procedure described in Bodnar, Dickhaus (2014). See the
#' vignette \code{vignette('fdr-test',package='MHTcop')} for a detailed explanation of the example.
#'
#' @param p The vector of marginal p-values
#' @param cop The dependency model for the p-values (for example
#'   copula::copClayton)
#' @param m0Lower A lower bound on the number of true null hypotheses (i.e.
#'   m0Lower is a reasonable lower bound for the number of true null
#'   hypotheses), \eqn{1 \le m0Lower \le length(p)}
#' @param alpha The desired FDR level
#' @param num.reps The number of samples to draw for the Monte-Carlo integration
#'   (default = 1e5)
#'
#' @return The adjusted p-values \code{p.adjusted} such that performing the test by
#'   rejecting the i-th hypothesis if and only if \code{p.adjusted[i]} \eqn{\le} \code{alpha} is a test at FDR
#'   level \code{alpha}
#' @export
#' @references T. Bodnar and T. Dickhaus (2014). False discovery rate control under Archimedean copula. \emph{Electronic Journal of Statistics} Volume 8, Number 2 (2014), 2207-2241.
#'
#' @examples
#' #(Using p-values generated from the model (16))
#' library(copula)
#' set.seed(1)
#' m <- 20
#' m0 <- 0.8*m
#' p_values <- rCopula(1,onacopulaL(copClayton,list(1,1:20)))
#' mu<-runif(m-m0, min=-1, max=-1/2)
#' p_values[1,(m0+1):m]<-pnorm(sqrt(m)*mu+qnorm(p_values[(m0+1):m]),0,1)
#' ac_fdr.test(p_values,setTheta(copClayton,1),m0,0.05,1e4)$test

ac_fdr.test <- function(p,cop,m0Lower,alpha=0.05,num.reps=1e5) {
  m <- length(p)

  if(m0Lower >= m || m0Lower <= 1) stop(paste("m0Lower should be larger than 1 and less than ", m))

  delta <- sapply(m0Lower:m,function(m0)ac_fdr.calc_delta(cop,m,m0,alpha=alpha,num.reps=num.reps))
  delta <- max(delta);

  list(test = stats::p.adjust(p,'BH') <= alpha/delta, delta = delta)
}
