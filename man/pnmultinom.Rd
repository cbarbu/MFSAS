\name{pnmultinom}
\Rdversion{1.1}
\alias{pnmultinom}

\title{The Cumulative Distribution Function for the Negative Multinomial Distribution}

\description{Compute cumulative probability for the negative multinomial distribution.}

\usage{pnmultinom(x, m, prob)}

\arguments{
  \item{x}{Vector of length \code{k-1} for the failures where k is the number of 
    classes.}


  \item{m}{The target number of successful trials.
    Must be a strictly positive integer.}

  \item{prob}{Numeric non-negative vector of length \code{k-1} specifying the 
     probability for the \code{k-1} classes of  failures.} 
}

\details{The negative multinomial distribution is used for sequential sampling 
  \emph{with} replacement. Suppose that the lot has \code{k-1} different types of  
 failures, with corresponding probabilities 
 \code{p[1], ..., p[k-1]} in each trial.  Let \code{X[1], ..., X[k-1]} denote 
  the number of failures of each type that are selected in a sequence of trials 
 before a target number \code{m} of successes is reached. Then 
 \code{pnmultinom(x, m, prob = p)} is the cumulative probability:

\deqn{\begin{array}{ll}
 &P(X[1] <= x[1], \dots, X[k-1] <= x[k-1]) \\\\
 & = \sum\limits_{y[1]=0}^{x[1]} \dots 
\sum\limits_{y[k-1]=0}^{x[k-1]} P(X[1] = y[1] , \dots, X[k-1] = y[k-1] )\\\\
&=\sum\limits_{y[1]=0}^{x[1]} \dots \sum\limits_{y[k-1]=0}^{x[k-1]} 
\frac{\left(\sum\limits_{i=1}^{k-1}y[i] + m -1\right)!}{(m-1)!y[1]! \dots  y[k-1]!}
p[1]^{y[1]} \dots p[k-1]^{y[k-1]}\left(1-\sum\limits_{i=1}^{k-1}p[i]\right)^m
\end{array}
}{%
P(X[1] <= x[1], ..., X[k-1] <= x[k-1]) = sum(...(sum(P(X[1] = y[1] , 
..., X[k-1] = y[k-1],  y[1]=0 ...x[1])...) y[k]=0 ...x[k])}
\deqn{}{
=sum(...(sum((sum(y) + m - 1)!/[(m-1)!y[1]!... y[k]!] p[1]^y[1]
  ... p[k-1]^y[k-1](1-sum(p)^m), y[1]=0 ...x[1])...)y[k-1]=0 ...x[k-1])}

}

\value{
 \code{pnmultinom} gives the value of the cumulative distribution function.
  Invalid \code{m} or \code{prob} will stop running with a warning.
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Johnson, N. L., Kotz, S., and Balakrishnan, N. (1997).
\emph{Discrete Multivariate Distributions}, Wiley, New York. 
}

\author{Aaron Childs and Yalin Chen}

%\note{Here the Dirichlet \bold{D} function is used to calculate the lower tail 
 %probability for the negative multinomial distribution. See reference.}

\section{Source}{\code{pnmultinom} is computed using recursive algorithms for the Dirichlet \bold{D} 
  function given in \cr

  Sobel, M., Uppuluri, V.R.R, and Frankowski, K. (1985). Dirichlet integrals of type-2 and their application.  
    In \emph{Selected Tables in Mathematical Statistics} \bold{9},  
    American Mathematical Society, Providence, Rhode Island. 
   
   Sobel, M. and Frankowski, K. (2004). Extensions of Dirichlet integrals: their 
  computation and probability applications. In \emph{Gupta, A. K. and Nadarajah, S. (eds) 
  Handbook of Beta Distribution and its applications,} 319-360, Marcel Dekker, New York.   
} 

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pmultinom}}, \code{\link{pmultihyper}}, \code{\link{pnmultihyper}}
}
\examples{
 X <- c(4,5,6)
 m <-3
 pr <-c(0.10,0.15,0.18)
 pnmultinom(x = X, m = m, prob = pr)
}
\keyword{methods}
