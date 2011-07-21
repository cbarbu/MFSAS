\name{pnmultihyper}
\Rdversion{1.1}
\alias{pnmultihyper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Cumulative Distribution Function for the Negative Multivariate Hypergeometric Distribution}

\description{Compute cumulative probability for the negative
  multivariate  hypergeometric distribution.
}

\usage{
pnmultihyper(x, m, M, N)
}
 
\arguments{
  \item{x}{Vector of length \code{k-1} for the failures where \code{k} is 
    the number of classes.}


  \item{m}{The target number of successful trials.
           Must be a strictly positive integer.}

  \item{M}{Vector of length \code{k-1} containing the total number of each type of 
     failure in the lot.}

  \item{N}{Total lot size from which the sample is drawn.}
}
\details{The negative multivariate hypergeometric distribution is used for sequential
  sampling from a finite lot \emph{without} replacement. Suppose that the 
  lot of size \code{N} has \code{k-1} different types of failures represented 
  \code{M[1], ..., M[k-1]} times. Let \code{X[1], ..., X[k-1]} denote the number 
  of failures of each type that are selected in a sequence of trials 
  before a target number \code{m} of successes is reached. Then
  \code{pnmultihyper(x, m, M, N)} is the cumulative probability:
 
  \deqn{\begin{array}{ll}
    &P(X[1] <= x[1], \dots, X[k-1] <= x[k-1])\\\\
    &= \sum\limits_{y[1]=0}^{x[1]} \dots \sum\limits_{y[k-1]=0}^{x[k-1]} 
    P(X[1] = y[1] , \dots, X[k-1] = y[k-1])\\\\
    &=\sum\limits_{y[1]=0}^{x[1]} \dots \sum\limits_{y[k-1]=0}^{x[k-1]} 
    \frac{{M[1] \choose y[1]}\cdots {M[k-1] \choose y[k-1]}
  {N-\sum\limits_{i=1}^{k-1}M[i] \choose m-1}}{{N \choose m-1+
   \sum\limits_{i=1}^{k-1}y[i]}}
   \frac{N - m +1 -\sum\limits_{i=1}^{k-1}M[i]}{N - m +1 -\sum\limits_{i=1}^{k-1}y[i]}
  \end{array}
  }{%
    P(X[1] <= x[1], ..., X[k-1] <= x[k-1]) = sum(...(sum(P(X[1] = y[1] , 
   ..., X[k-1] = y[k-1]),y[1]=0 ...x[1])...)y[k]=0 ...x[k])
  }
   \deqn{}{%
    = sum(...(sum([choose(M[1],y[1]), ... ,choose(M[k-1], y[k-1])
     choose(N-sum(M),m-1)/choose(N,(sum(y)+m-1))]}
    \deqn{}{%
     *[N-m+1-sum(M)]/[N-m+1-sum(y)],
      y[1]=0 ... x[1])...)y[k-1]=0 ... x[k-1])}

where the sum is over all values of \code{y} such that 
 \eqn{y[1]+y[2]+...+y[k-1] <= N-m}.
}

\value{
\code{pnmultihyper} gives the value of the cumulative distribution function.
 Invalid \code{M} or \code{N} will stop running with a warning.
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Johnson, N. L., Kotz, S., and Balakrishnan, N. (1997).
\emph{Discrete Multivariate Distributions}, Wiley, New York. 
}

\author{Aaron Childs and Yalin Chen}


\section{Source}{\code{pnmultihyper} is computed using recursive algorithms 
  for the Dirichlet \bold{HD} function given in \cr

  Childs, A. (2010). Vector extensions of the Dirichlet HC and HD
  functions, with applications to the sharing problem. 
  \emph{Methodology and Computing in Applied Probability}
   \bold{12,} 91-109.

  Sobel, M. and Frankowski, K. (1995). Hypergeometric analogues of multinomial type-2
    problems via Dirichlet methodology. \emph{Congressus Numerantium} \bold{106}, 171-191. 
} 



\seealso{
\code{\link{pmultinom}}, \code{\link{pmultihyper}}, \code{\link{pnmultinom}}.
}
\examples{
 X <- c(2,3,4,1)
 m <- 5 
 M <-c(5,7,8,3)
 N <-30
cdf <- pnmultihyper(X, m, M, N)
cdf
}
\keyword{methods}
