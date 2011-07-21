\name{pmultinom}
\Rdversion{1.1}
\alias{pmultinom}
\title{The Cumulative Distribution Function for the Multinomial Distribution}
\description{Compute cumulative probability for the multinomial distribution.}
\usage{
 pmultinom(x, size, prob)
}
\arguments{
  \item{x}{Vector of length \code{k-1} of non-negative integers 
   where \code{k} is the number of classes.}

  \item{size}{The sample size}

  \item{prob}{Numeric non-negative vector of length \code{k-1} 
    specifying the probability for \code{k-1} of the classes.}
}
\details{The multinomial distribution is used for sampling \emph{with} replacement, or
  if the lot is large compared to the sample size. If a sample
  of size  \code{n} is drawn from a lot whose \code{k} classes have probabilities
  \eqn{p[1],..., p[k-1], p[k]}, let \eqn{X[1], ..., X[k-1], X[k]} denote the 
   number of observations drawn from each of the \code{k} classes. Then 
  the cumulative probability \code{pmultinom(x, size=n, prob=p)} is given by   

  \deqn{\begin{array}{ll}
    &P(X[1] <= x[1],\dots, X[k-1] <= x[k-1])\\\\
    &\quad = \sum\limits_{y[1]=0}^{x[1]} \dots \sum\limits_{y[k-1]=0}^{x[k-1]}
     P(X[1] = y[1], \dots, X[k-1] = y[k-1])\\\\
    &\quad =\sum\limits_{y[1]=0}^{x[1]} \dots \sum\limits_{y[k-1]=0}^{x[k-1]}
     \frac{n!}{y[1]! \dots y[k-1]!\left(n-\sum\limits_{i=1}^{k-1}y[i]\right)!}
    p[1]^{y[1]} \dots p[k-1]^{y[k-1]} \\\\
    & \quad \quad \quad \quad \quad \quad 
     *\left(1-\sum\limits_{i=1}^{k-1}p[i]\right)^
      {\left(n-\sum\limits_{i=1}^{k-1}y[i]\right)}
  \end{array}
  }{%
    P(X[1] <= x[1], ..., X[k-1] <= x[k-1])= sum(...( sum(
      P(X[1] = y[1], ..., X[k-1] = y[k-1]), y[1]=0...x[1]) ... )y[k-1]=0 ... x[k-1])}
  \deqn{}{
    = sum(...( sum(
\{n!/(y[1]! ... y[k-1]![n - sum(y)]!)\}*(p[1]^y[1]...p[k-1]^y[k-1]*[1- sum(p)]^[n- sum(y)]), 
   y[1]=0...x[1]) ... )y[k-1]=0 ... x[k-1])}

  
where the sum is over all values of \code{y} such that 
\eqn{y[1]+y[2]+...+y[k-1] <=n }.   
}

\value{
\code{pmultinom} gives the value of the cumulative distribution function.
 Invalid arguments will stop running.
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Johnson, N. L., Kotz, S., and Balakrishnan, N. (1997).
\emph{Discrete Multivariate Distributions}, Wiley, New York. 
}

\author{ Aaron Childs and Yalin Chen}

\section{Source}{\code{pmultinom} is computed using recursive algorithms for 
  the Dirichlet \bold{J} function given in \cr

  Sobel, M., Uppuluri, V. R. R., and Frankowski, K. (1977). Dirichlet distributions  
  type-1. In \emph{Selected Tables in Mathematical Statistics} \bold{4},
  American Mathematical Society, Providence, Rhode Island.

  Sobel, M. and Frankowski, K. (2004). Extensions of Dirichlet integrals: their 
  computation and probability applications. In \emph{Gupta, A.K. and Nadarajah, S. (eds) 
  Handbook of Beta Distribution and its applications,} 319-360, Marcel Dekker, New York.   
  }

\seealso{
\code{\link{pmultihyper}}, \code{\link{pnmultinom}}, \code{\link{pnmultihyper}}
}

\examples{
X <- c(2,3)
n <- 20
pr <- c(0.12, 0.15)
cdf <- pmultinom(x = X, size = n, prob = pr)
cdf

}
\keyword{methods}
