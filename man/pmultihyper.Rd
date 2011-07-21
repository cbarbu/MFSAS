\name{pmultihyper}
\Rdversion{1.1}
\alias{pmultihyper}

\title{The Cumulative Distribution Function for the Multivariate Hypergeometric Distribution}

\description{Compute cumulative probability for the multivariate 
 hypergeometric distribution.}

\usage{pmultihyper(x, n, M, N)}


\arguments{
  \item{x}{Vector of length \code{k-1} of non-negative integers where \code{k} is the number of
   classes.}

  \item{n}{The sample size.}

  \item{M}{Vector of length \code{k-1} containing the total number of objects 
   in each of \code{k-1} of the classes.}
 
  \item{N}{The size of the lot from which the sample is drawn.}
}

\details{The multivariate hypergeometric distribution is used for sampling 
 from a finite lot \emph{without} replacement. If a sample
  of size  \code{n} is drawn from a lot of size \code{N} which 
  has \eqn{M[i]} objects of type \eqn{i} (for \eqn{i= 1,2, ..., k}),
  let \eqn{X[i]} be the number of objects of type \eqn{i} in the sample 
   (for \eqn{i= 1,2, ..., k}). Then
  the cumulative probability \code{pmultihyper(x, n, M, N)} is given by,  \cr

  \deqn{\begin{array}{ll}
   &P(X[1] <= x[1],\dots, X[k-1] <= x[k-1]) \\\\
  &=\sum\limits_{y[1]=0}^{x[1]}\dots\sum\limits_{y[k-1]=0}^{x[k-1]}
  P(X[1]= y[1],\dots, X[k-1] = y[k-1])\\\\
  &=\sum\limits_{y[1]=0}^{x[1]}\dots\sum\limits_{y[k-1]=0}^{x[k-1]}
   \frac{{M[1] \choose y[1]} \dots 
  {M[k-1]\choose y[k-1]}
  {N-\sum\limits_{i=1}^{k-1}M[i]\choose n-\sum\limits_{i=1}^{k-1}y[i]}}{{N \choose n}}
  \end{array}
  }{%
   P(X[1] <= x[1], ... , X[k-1] <= x[k-1]) = 
   sum(...(sum(P(X[1] = y[1],..., X[k-1] = y[k-1]), 
  y[1]=0 ... x[1])...)y[k-1]=0 ... x[k-1])}
 \deqn{}{ = sum(...(sum(choose(M[1],y[1]), ... ,choose(M[k-1], y[k-1])
  choose(N-sum(M),n-sum(y))/choose(N,n), y[1]=0 ... x[1])
  ...)y[k-1]=0 ... x[k-1]) 
        } 

where the sum is over all values of \code{y} such that 
 \eqn{y[1]+y[2]+...+y[k-1] <= n} and 
\eqn{n - (y[1]+y[2]+...+y[k-1]) <= N - (M[1]+M[2]+...+M[k-1])}.
}

\value{
 \code{pmultihyper} gives the value of the cumulative distribution function.
Invalid arguments will stop running.
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Johnson, N. L., Kotz, S., and Balakrishnan, N. (1997).
\emph{Discrete Multivariate Distributions}, Wiley, New York.
}

\author{Aaron Childs and Yalin Chen}

%\note{Here the Dirichlet \bold{HJ} function is used to calculate the lower tail probability 
% for the multivariate hypergeometric distribution. See reference.}

\section{Source}{\code{pmultihyper} is computed using recursive algorithms for the Dirichlet \bold{HJ}  
  function given in \cr

  Sobel, M. and Frankowski, K. (1994). Hypergeometric analogues of multinomial type-1
  Dirichlet problems. \emph{Congressus Numerantium} \bold{101}, 65-82.
} 

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pmultinom}}, \code{\link{pnmultinom}}, \code{\link{pnmultihyper}}.}

\examples{
X <- c(1,3,4)
n <- 15
M <- c(8, 10, 14)
N <- 50
pr <- pmultihyper(X, n, M, N)
pr
}
\keyword{methods}

