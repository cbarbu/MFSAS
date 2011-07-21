\name{assess.multi}
\Rdversion{1.1}
\alias{assess.multi}

\title{Utility Function for Assessing Multilevel Sampling Plans}

\description{Assess whether the \code{k}-level fixed or sequential sampling plan can meet 
 the specified \emph{Producer's Risk Point (PRP)} and/or \emph{Consumer's Risk Point (CRP)}.}

\usage{
assess.multi(rn, n = 30, m, N = 100, PRP, CRP, type=c("multinomial", 
 "hypergeom"), stype=c("fixed", "sequential"), print = TRUE)
}

\arguments{
   \item{rn}{A vector of length \code{k-1} of rejection numbers for fixed sampling, 
     or cell quotas for the defective items in sequential sampling.}

  \item{n}{Sample size; applicable for \code{stype="fixed"}.}

  \item{m}{The cell quota for good items; applicable for \code{stype="sequential"}.}

  \item{N}{The lot size from which the sample is drawn; 
    applicable for \code{type="hypergeom"}.}

  \item{PRP}{The Producer's Risk Point in the form of a two part
    numeric vector \code{(pd, pa)}. The first part
    \code{pd}, a vector of length \code{k-1}, specifies the quality level 
    at which to evaluate the plan. The second part \code{pa}, 
    indicates the \emph{minimum} probability of acceptance to be 
    achieved by the plan.}

  \item{CRP}{The Consumer's Risk Point in the form of a two part
    numeric vector \code{(pd, pa)}. The first part
    \code{pd}, a vector of length \code{k-1}, specifies the quality level 
    at which to evaluate the plan. The second part \code{pa}, 
    indicates the \emph{maximum} probability of acceptance 
    to be achieved by the plan.}

  \item{type}{The type of distribution on which the sampling plan is based.
   Possible values are \code{"multinomial"} and \code{"hypergeom"}. The default 
   is \code{"multinomial"}.}

  \item{stype}{The type of sampling. Possible values are \code{"fixed"} 
   and \code{"sequential"}. The default is \code{"fixed"}.}

  \item{print}{Logical, indicating whether or not a summary of the
    assessment should be printed. }
}

\details{
   Typical usages are:
   \preformatted{
    assess.multi(rn, n, PRP, CRP)
    assess.multi(rn, m, PRP, CRP, stype ="sequential")
    assess.multi(rn, n, N, PRP, CRP, type="hypergeom")
    assess.multi(rn, m, N, PRP, CRP, type = "hypergeom", stype="sequential")
    }
 In the first form, the default \code{type} \code{"multinom"} and the default 
   \code{stype}  \code{"fixed"} are used. 

  The second form is based on the negative multinomial distribution.
  
 The third form uses a default \code{stype} of \code{"fixed"}
  and is based on the multivariate hypergeometric distribution.
 
 The fourth form is based on the negative multivariate hypergeometric distribution.

 The cell quota \code{m} for the good items  must be provided in both second and fourth forms.

In both third and fourth cases, the lot size \code{N} needs to be 
  specified, and \code{pd} in \code{PRP} and \code{CRP} is the vector of
  the proportions of defectives in the lot. Since \code{pd*N} gives a vector
  containing the actual numbers of each type of defective in the lot, 
  all of its entries must be integers.       
}

\value{
The function will return the result of whether the plan meets the acceptance requirement(s), 
 along with the actual acceptance probability achieved by the sampling plan. In the case of 
 sequential sampling, average sampling numbers \code{ASNp} and \code{ASNc} for the quality
 levels in \code{PRP} and \code{CRP} are also returned.
}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Schilling, E. G. and Neubauer, D. V. (2009). \emph{Acceptance Sampling in 
   Quality Control}, Second Edition, CRC Press, New York.

}



\author{Aaron Childs and Yalin Chen}


\section{Source}{For sequential sampling, the average sampling number \code{ASN} 
  is computed using  algorithms for the Dirichlet \code{D} function 
 (for \code{type="multinomial"}) or \bold{HD} function (for \code{type="hypergeom"}),
  together with equation (5.30) in 
  Sobel, Uppuluri, and Frankowski (1985) (for \code{type="multinomial"}), or  
  equation (5.3) in Childs (2010) (for \code{type="hypergeom"}).\cr

  Childs, A. (2010). Vector extensions of the Dirichlet HC and HD
   functions, with applications to the sharing problem. 
   \emph{Methodology and Computing in Applied Probability}
   \bold{12}, 91-109.

  Sobel, M., Uppuluri, V.R.R, and Frankowski, K. (1985), 
   Dirichlet integrals of type-2 and their application.  
   In \emph{Selected Tables in Mathematical Statistics} \bold{9},  
    American Mathematical Society, Providence, Rhode Island. 
} 


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{find.multi.plan}}, \code{\link{Ocmult-class}}
}
\examples{
assess.multi(n=30, rn=c(2,2,3), PRP = c(0.05,0.06, 0.08, 0.95), 
 CRP = c(0.15,0.18, 0.20, 0.075))
assess.multi(rn=c(7,8), m=5, PRP = c(0.1,0.05, 0.95), CRP = c(0.2,0.15, 0.075),
 type="multinomial", stype="seq")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
