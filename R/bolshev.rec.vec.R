#' Distribution function of the order statistics of i.i.d. uniform random variables
#'
#' \code{bolshev.rec.vec} is a vectorized and unrolled implementation of the Bolshev recursion described in Shorack, Wellner (1986)
#' which can be utilized to calculate probabilities for order statistics of i.i.d. uniform random variables.
#'
#' Denote by \eqn{U_1,\cdots,U_n} n i.i.d. uniform random variables on \eqn{[0,1]}. Denote by \eqn{U_{1:n},\cdots,U_{n:n}} their order statistics.
#' Then the return value \code{p} contains the probabilities \deqn{p[i,j] = P\left(\bigcap\limits_{k=i}^n\left\{m[n-k+1,j] \le U_{k:n}\right\}\right)}{p[i,j] = P(\forall k=i,\cdots,n: m[n-k+1,j] \le U_{k:n})}
#' @param m matrix whose columns are p-values sorted in descending order
#' @return matrix p containing the calculated probabilities
#' @references G. R. Shorack and J. A. Wellner (1986). Empirical Processes with Applications to Statistics
#'
#' @export
#' @examples
#' bolshev.rec.vec(cbind(rev(c(0.7,0.8,0.9))))
#' #result: c(0.016, 0.079, 0.271)
#' #monte carlo simulation
#' sim <- function(v) mean(replicate(1e4,all(v <= sort(runif(3)))))
#' set.seed(0)
#' c(sim(c(0.7,0.8,0.9)),sim(c(0,0.8,0.9)),sim(c(0,0,0.9)))
#' #similar result: c(0.0176, 0.0799, 0.2709)

bolshev.rec.vec <- function (m)
{
  dim.row <- nrow(m)
  dim.col <- ncol(m)
  s <- 0:(dim.row-1)
  summands <- m
  summands[-1,] <- 0
  for (k in 2:dim.row) {
    Fk <- 1 - .colSums(summands,dim.row,dim.col,na.rm=TRUE)
    summands <- (k / (k-s))  * summands * m
    summands[k,] <- k * Fk * m[k,]
  }
  ret <- 1 - matrixStats::colCumsums(summands)
  ret[dim.row:1,]
}
