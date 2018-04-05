#' Evaluate the inverse Laplace-Stieltjes transform of a copula's generator
#'
#' \code{cdf.Z} evaluates the inverse Laplace-Stieltjes transform of the generator of the copula \code{cop} at \code{z}. Not: The evaluated mapping is a distribution function.
#'
#' @param cop The copula
#' @param z Argument to the inverse Laplace-Stieltjes transform of the copula's generator
#'
#' @export
#' @keywords internal
#'
cdf.Z <- function(cop,z) {
  switch(cop@name,
         Clayton = {
           stats::pgamma(z, shape=1/cop@theta, scale = 1)
         },
         Gumbel = {
           res <- tryCatch(stabledist::pstable(z,1/cop@theta,1,cos(pi/(2 * cop@theta))^cop@theta,0,1), warning=function(w) 0)
           if(res <= 5*1e-3) res <- sum(sample.Z(cop,1e6)<=z)/1e6
           res
         },
         AMH = {
           1-cop@theta^floor(z)
         },
         #Frank = {

         #},
         #Joe = {

         #},
         {
           stop(paste("Simulation of the mixing variable of the mixture representation not implemented for copula type",
                      sQuote(cop@name)))
         })
}
