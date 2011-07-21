\name{find.multi.plan}
\Rdversion{1.1}
\alias{find.multi.plan}

\title{Utility Function for Finding Multilevel Sampling Plans}

\description{Find the \code{k}-level sampling plan with the 
   smallest (expected) sample size such that the specified 
   \emph{Producer's Risk Point (PRP)} and \emph{Consumer's Risk 
   Point (CRP)} are met.}

\usage{
 find.multi.plan(PRP, CRP, N=100, type=c("multinomial", "hypergeom"), 
 stype=c("fixed", "sequential"))
 }

\arguments{
  \item{PRP}{The Producer's Risk Point in the form of a two part
   numeric vector of the form \code{(pd, pa)}. The first part
   \code{pd}, a vector of length \code{k-1}, specifies the quality 
   level at which to evaluate the plan. 
   The second part \code{pa}, indicates the \emph{minimum} 
   probability of acceptance to be achieved by the plan.}

  \item{CRP}{The Consumer's Risk Point in the form of a two part
   numeric vector of the form \code{(pd, pa)}. The first part
  \code{pd}, a vector of length \code{k-1}, specifies the quality level 
   at which to evaluate the plan. The second part \code{pa}, 
   indicates the \emph{maximum} probability of acceptance 
   to be achieved by the plan.}

  \item{N}{The lot size from which the sample is drawn; 
    applicable for \code{type="hypergeom"}.}

  \item{type}{The distribution on which the sampling plan is based.
    Possible values are \code{"multinomial"} and \code{"hypergeom"}.
    The default is \code{"multinomial"}.}

  \item{stype}{The type of sampling. Possible values are \code{"fixed"} 
   and \code{"sequential"}. The default is \code{"fixed"}.}
}

\details{
   Typical usages are:
   \preformatted{
     find.multi.plan(PRP, CRP)
     find.multi.plan(PRP, CRP, stype ="sequential")
     find.multi.plan(PRP, CRP, N, type="hypergeom")
     find.multi.plan(PRP, CRP, N, type="hypergeom", stype="sequential")
    }
 In the first form, the default \code{type} \code{"multinomial"} and the
   default \code{stype} \code{"fixed"} are used. 

  The second form is based on the negative multinomial distribution.
    
 The third form uses a default \code{stype} of \code{"fixed"}
  and is based on the multivariate hypergeometric distribution. 
 
 The fourth form is based on the negative multivariate hypergeometric distribution.

In both third and fourth cases, the lot size \code{N} needs to be 
  specified, and \code{pd} in \code{PRP} and \code{CRP} is the vector of
  the proportions of defectives in the lot. Since \code{pd*N} gives a vector
  containing the actual numbers of each type of defective in the lot, 
  all of its entries must be integers.       
}

\value{
      The values returned are
  \item{n}{The smallest possible sample size for \code{stype="fixed"}.}
 
  \item{rn}{Vector of length \code{k-1} of rejection numbers for 
    fixed sampling, or cell quotas for the defective items in sequential sampling.}

  \item{m}{The smallest possible number of good items for \code{stype="sequential"}.}

  \item{p.PRP}{The actual probability of acceptance at the producer's quality level
     for the sampling plan.}

  \item{p.CRP}{The actual probability of acceptance at the consumer's quality level 
     for the sampling plan.}

  \item{ASNp}{The average sampling number 
    at the producer's quality level for \code{stype="sequential"}.}

  \item{ASNc}{The average sampling number 
    at the consumer's quality level for \code{stype="sequential"}.}
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Schilling, E. G. and Neubauer, D. V. (2009). \emph{Acceptance Sampling in 
   Quality Control}, Second Edition, CRC Press, New York.

}

\author{Aaron Childs and Yalin Chen}

\section{Source}{For sequential sampling, the average sampling numbers \code{ASNp} and \code{ASNc} 
  are computed using algorithms for the Dirichlet \code{D} function 
 (for \code{type="multinomial"}) or \bold{HD} function (for \code{type="hypergeom"}),
  together with equation (5.30) in 
  Sobel, Uppuluri,  and Frankowski (1985) (for \code{type="multinomial"}), or  
  equation (5.3) in Childs (2010) (for \code{type="hypergeom"}).\cr

  Childs, A. (2010). Vector extensions of the Dirichlet HC and HD
   functions, with applications to the sharing problem. 
   \emph{Methodology and Computing in Applied Probability}
   \bold{12}, 91 - 109.

  Sobel, M., Uppuluri, V. R. R., and Frankowski, K. (1985). 
   Dirichlet integrals of type-2 and their application.  
   In \emph{Selected Tables in Mathematical Statistics} \bold{9},  
    American Mathematical Society, Providence, Rhode Island. 

  } 

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Ocmult-class}}, \code{\link{Ocmult}}, \code{\link{assess.multi}}.}

\examples{
 find.multi.plan(PRP=c(0.03, 0.05,  0.8), CRP=c(0.15, 0.16, 0.1), stype="seq")

 find.multi.plan(PRP=c(0.06, 0.04, 0.06, 0.8), CRP=c(0.14, 0.16, 0.2, 0.1), 
  N= 100,type = "h", stype="seq")
}

\keyword{methods}
