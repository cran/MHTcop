#' Generate a sample from the inverse Laplace-Stieltjes transform of a copula's generator
#'
#' \code{sample.Z} generates a sample of size \code{n} from the inverse Laplace-Stieltjes transform of the generator of the copula \code{cop}. For further details see \url{https://doi.org/10.1016/j.csda.2008.05.019} (especially table 1).
#'
#' @param cop The copula
#' @param n The desired sample size
#'
sample.Z <- function(cop,n) {
  switch(cop@name,
         Clayton = {
           stats::rgamma(n, shape=1/cop@theta, scale = 1)
         },
         Gumbel = {
           #QRM::rstable(num.reps,1/cop@theta) * (cos(pi/(2 * cop@theta)))^(cop@theta)
           stabledist::rstable(n,1/cop@theta,1,cos(pi/(2 * cop@theta))^cop@theta,0,1)
         },
         AMH = {
           #inversion sampling
           ceiling(log(1-stats::runif(n),base=cop@theta))
         },
         Frank = {
           sample.discrete(function(k) (1-exp(-cop@theta))^k / (k*cop@theta),n)
         },
         Joe = {
           sample.discrete(function(k) (-1)^(k+1)*choose(1/cop@theta,k),n)
         },
         {
           stop(paste("Simulation of the mixing variable of the mixture representation not implemented for copula type",
                      sQuote(cop@name)))
         })
}

#' Generate a sample from a discrete distribution
#'
#' \code{sample.discrete} generates a sample of size \code{n} given its density function \code{df}
#'
#' @param df The density function - It is assumed that the support is a subset of the natural numbers
#' @param n The desired sample size

sample.discrete <- function(df,n) {
  x <- stats::runif(n)
  idx <- seq_along(x)
  n <- 1
  curr <- 0
  while(length(idx) > 0) {
    df.eval <- df(n)
    if(df.eval < .Machine$double.eps) {
      x[idx] <- n
      break;
    }
    curr <- curr + df.eval
    idx.remove <- which(x[idx] <= curr)
    x[idx[idx.remove]] <- n
    n <- n+1
    if(length(idx.remove)>0) idx <- idx[-idx.remove]
  }
  as.integer(x)
}
