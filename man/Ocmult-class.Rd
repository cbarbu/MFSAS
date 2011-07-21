\name{Ocmult-class}
\Rdversion{1.1}
\docType{class}
\alias{Ocmult-class}
\alias{Ocmult.multinomial-class}
\alias{Ocmult.hypergeom-class}

\title{Class "Ocmult" }
\description{"Operating Characteristic" function of the class \code{"Ocmult"}  
  provides methods for creating, plotting and printing \code{k}-level acceptance sampling 
  plans based on the Multinomial, Negative Multinomial (\code{"Ocmult.multinomial"}),  
  Multivariate Hypergeometric and Negative Multivariate Hypergeometric 
  (\code{"Ocmult.hypergeom"}) distributions.}

\section{Objects from the Class}{The \code{"Ocmult"} class is a virtual
  class: No objects may be created from it.

  However, objects from the derived classes \code{Ocmult.multinomial} 
  and \code{Ocmult.hypergeom} can be created using the create function
  \code{\link{Ocmult}}.
  }

\section{Slots}{
	 \describe{
    \item{\code{n}:}{Object of class \code{"numeric"}.  
      The sample size; applicable for \code{stype="fixed"}.}

    \item{\code{m}:}{Object of class \code{"numeric"}.  
      The cell quota for good items; applicable for \code{stype="sequential"}.}

    \item{\code{rn}:}{Object of class \code{"numeric"}.  A vector of length \code{k-1} of 
      rejection numbers for fixed sampling, or cell quotas for the defective items in 
      sequential sampling.}

    \item{\code{pd}:}{Object of class \code{"matrix"}.  A matrix whose rows are vectors 
      containing the proportions of each type of defective.}

    \item{\code{pa}:}{Object of class \code{"numeric"}.  A numeric vector contains the 
      probabilities of acceptance according to the proportion of defectives in the 
      rows of \code{pd}.}

    \item{\code{ASN}:}{Object of class \code{"numeric"}.  Only for sequential sampling.
      A numeric vector containing average sampling numbers according 
      to the proportion of defectives in the rows of \code{pd}.}

    \item{\code{stype}:}{Object of class \code{"character"}.  The type of sampling. 
      Possible values are \code{"fixed"} and \code{"sequential"}.}

    \item{\code{type}:}{Object of class \code{"character"}.  The type of distribution on which the 
       plans are based. Possible values are \code{"multinomial"} and \code{"hypergeom"}.}

    \item{\code{N}:}{Object of class \code{"numeric"}.  Only for class
      \code{"Ocmult.hypergeom"}. A number giving the lot size from
       which the sample is drawn.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "Ocmult")}: Show the details of 
      the sampling plan. }

    \item{summary}{\code{signature(object = "Ocmult")}: Summarize the sampling
      plan. Optional argument \code{"full"} (defaults to \code{FALSE}) will
      show the details at all quality values \code{"pd"} supplied when the
      object was created.}

    \item{plot}{\code{signature(x = "Ocmult.multinomial", y = "missing")},\cr
      \code{signature(x = "numeric", y = "Ocmult.multinomial")},\cr
      \code{signature(x = "Ocmult.hypergeom", y = "missing")}: \cr
      \code{signature(x = "numeric", y = "Ocmult.hypergeom")}: \cr
       Plot the OC curve for 2-level sampling plans.}

     \item{persp}{\code{signature(x = "Ocmult.multinomial")}, \cr
       \code{signature(x = "Ocmult.hypergeom")}: \cr
       Plot the OC surface for 3-level sampling plans.}

    \item{contour}{\code{signature(x = "Ocmult.multinomial")}, \cr
      \code{signature(x = "Ocmult.hypergeom")}: \cr
      Plot the OC contour curve for 3-level sampling plans.}
  }
}
\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.

Schilling, E. G. and Neubauer, D. V. (2009). \emph{Acceptance Sampling in 
   Quality Control}, Second Edition, CRC Press, New York.

}

\author{Aaron Childs and Yalin Chen}
%\note{ ~~further notes~~ }

 
\seealso{\code{\link{Ocmult}}, \code{\link{find.multi.plan}}, 
          \code{\link{assess.multi}}.} 

\keyword{classes}
