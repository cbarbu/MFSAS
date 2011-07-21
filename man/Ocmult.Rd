\name{Ocmult}
\Rdversion{1.1}
\alias{Ocmult}

\alias{show,Ocmult-method}
\alias{summary,Ocmult-method}
%\alias{initialize,Ocmult.multinomial-method}
%\alias{initialize,Ocmult.hypergeom-method}
\alias{plot,Ocmult.multinomial,missing-method}
\alias{plot,numeric,Ocmult.multinomial-method}
\alias{plot,Ocmult.hypergeom,missing-method}
\alias{plot,numeric,Ocmult.hypergeom-method}
\alias{contour,Ocmult.multinomial-method}
\alias{contour,Ocmult.hypergeom-method}
\alias{persp,Ocmult.multinomial-method}
\alias{persp,Ocmult.hypergeom-method}


\title{Operating Characteristics of Multilevel Acceptance Sampling Plans}

\description{Creating new objects from the \code{"Ocmult"} classes.}

\usage{Ocmult(rn, type=c("multinomial", "hypergeom"), 
       stype=c("fixed", "sequential"), ...)
 }

\arguments{
    \item{rn}{A vector of length \code{k-1} of rejection numbers for \code{k}-level 
     (\code{k} different types of items in the lot) fixed sampling, 
     or a vector of length \code{k-1} of cell quotas for 
     the defective items in \code{k}-level sequential sampling.}
    
    \item{type}{The type of distribution on which the plan is based. 
     Possible values are \code{"multinomial"} and \code{"hypergeom"}.
     The default is \code{"multinomial"}.}

    \item{stype}{The type of sampling. Possible values are \code{"fixed"} 
     and \code{"sequential"}. The default is \code{"fixed"}.}
    
    \item{\dots}{Additional arguments to be passed to the class
      generating function for each type. See Details for options.}  
  }

\details{
   Typical usages are:
   \preformatted{
     Ocmult(rn, n)
     Ocmult(rn, n, pd)
     Ocmult(rn, n, N, pd, type="hypergeom")
     Ocmult(rn, m, stype="sequential", pd)
     Ocmult(rn, m, N, pd, type="hypergeom", stype="sequential")
   }
  In the first and second forms, the default \code{type} \code{"multinomial"} 
   and the default \code{stype} \code{"fixed"} are used. The OC function is 
   calculated based on the proportion of defectives \code{pd}, whose default 
   values are used in the first form (and depend on the length of \code{rn}).  


 The third form is the OC function based on the multivariate hypergeometric
  distribution. In this case, the lot size \code{N} needs to be 
  specified, and \code{pd} is a matrix whose rows are vectors containing 
  the proportions of each type of defective. Since \code{pd*N} is a 
  matrix containing the actual numbers of each type of defective in the 
  lot, all of its entries must be integers. If \code{N} is
  not specified, it takes a default value of \code{N=100}.
  
 The fourth form uses a default \code{type} of \code{"multinomial"}. 
  Its OC function is based on the negative multinomial distribution, 
  hence the cell quota \code{m} for good items must be specified.  
   
 In the fifth form, the OC function is based on the negative multivariate 
  hypergeometric distribution. The cell quota \code{m} for good items and 
  the lot size \code{N} need to be specified. 
}

\value{
An object from \code{Ocmult-class} returns the class
  \code{Ocmult.multinomial} or \code{Ocmult.hypergeom}.\cr
There is a logic argument \code{detail} in the function \code{summary}.
  If \code{detail=TRUE}, all of the information for the object is shown.
  For sequential sampling the average sampling number (\code{ASN}) is also
  provided. The default value for this argument is \code{detail=FALSE}.

}

\references{Childs, A., and Chen, Y. (2011). Multilevel Fixed and Sequential
  Acceptance Sampling: The R Package MFSAS. \emph{Journal of Statistical
  Software}, \bold{43}(6), 1-20. URL http://www.jstatsoft.org/v43/i06/.}

\author{Aaron Childs and Yalin Chen}


\seealso{\code{\link{Ocmult-class}}, \code{\link{find.multi.plan}}, 
          \code{\link{assess.multi}}.} 
	
\examples{
px <- as.matrix(expand.grid(seq(0,0.5, 0.1),seq(0,0.5, 0.1)))
p.multinom <- Ocmult(n=30,rn=c(3,4), pd=px)
summary(p.multinom ,detail=TRUE )

p.multih<- Ocmult(c(3,4),n=15,  N=100, type="h")
summary(p.multih,detail=TRUE)
persp(p.multih)

p.nmultinom <- Ocmult(c(3,4), m=5, stype="s")
p.nmultinom 
summary(p.nmultinom ,detail=TRUE)
persp(p.nmultinom)
}
\keyword{methods}
\keyword{classes}

