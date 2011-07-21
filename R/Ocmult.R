#############################################
# Class definitions                         #
#                                           #
# Create a virtual class                    #
#                                           #
#############################################
setClass("Ocmult",
         representation(n="numeric",        # An integer - the sample size for fixed sampling
                        m="numeric",        # An integer - cell quota of good items for sequential sampling
                        rn="numeric",       # A vector of rejection numbers for fixed sampling or 
                                            # cell quotas of defectives for sequential sampling
                        pd="matrix",        # The proportions of each type of defective in the 
                                            # lot (by row)
                        pa="numeric",       # Probability of acceptance
                        asn="numeric",      # Average sampling number
                        stype="character",  # Type of sampling
                        type="character",   # Distribution type
                        "VIRTUAL"),

         validity=function(object){
                        if(any(is.na(object@rn))) 
                          return("The 'rn' vector is not allowed to contain missing values.")

                      # Check that the rejection numbers are reasonable         
                        if(any(object@rn < 1))
                        return("rejection number(s) 'rn' must be greater than 0.")
                    
                      if(any(is.na(object@pd))) 
                          return("'pd' is not allowed to contain missing values.")

                          l=ncol(object@pd)

                       # Check that the rows of pd have the same length as rn
                        if(length(object@rn)!=l)
                        return("The number of columns in 'pd' must be the same as the length of 'rn'.")

                       # Check that the proportions of defectives are reasonable
                          if(any(object@pd < 0))
                          return("The entries in 'pd' must be greater than or equal to 0.")

                          if(any(apply(object@pd,1,sum) > 1))
                          return("The row sums of 'pd' must be less than or equal to 1.")

                    return(TRUE)
                  }                     
     )

#############################################
# Create two classes                        #
#                                           #
#############################################
setClass("Ocmult.multinomial", 
         representation("Ocmult"                        
                       ), 
         prototype=list("Ocmult",
                         stype="fixed",  
                         type="multinomial",
                         n=30,
                         pd=as.matrix(seq(0, 0.1, by=0.01),nrow=1)
                        ),
         contains="Ocmult",
         validity=function(object){
                     
                      ## fixed sampling plan
                         if(object@stype=="fixed"){
                            if(is.na(object@n) )
                            return("The 'n' is not allowed to contain missing value.")

                         if(length(object@n)!= 1) 
                          return("The length of sample size 'n' should be equal to 1")

                      # Check that the sample size is reasonable
                        if(object@n <= 0) 
                        return("Sample size 'n' should be greater than 0.")

                      # Check that the rejection numbers are reasonable         
                        if(any(object@rn-1 > object@n))
                        return("Each rejection number in the vector 'rn' must be less than 'n+1'.")

                         }
                       # Check that the length of acceptance is reasonable for 
                      ## sequential sampling plan

                            if(object@stype=="sequential"){
                            if(is.na(object@m) )
                            return("The 'm' is not allowed to contain missing value.")

                         if(length(object@m)!= 1) 
                          return("The length of cell quota of good item 'm' should be equal to 1")

                      # Check that the cell quota is reasonable
                        if(object@m <= 0) 
                        return("cell quota 'm' should be greater than 0.")
                       }
                       return(TRUE)
                     }
        )

setClass("Ocmult.hypergeom", 
         representation("Ocmult", 
                        N="numeric" # An integer - the lot size
                        ),  
         prototype = list("Ocmult",
                           stype="fixed",  
                           type="hypergeom", 
                           N=100,
                           n=30,
                           pd=as.matrix(seq(0, 0.1, by=0.01),nrow=1)
                          ),
         contains="Ocmult",
         validity=function(object){
                         if(is.na(object@N)) return("N contains NA. ")

                       # Check that the lot size is reasonable   
                         if(length(object@N)>1) return("N must be of length 1.")                 
                         if(object@N < 1) return("N is less than 1. ")

                         if(length(object@rn) >= object@N)
                         return("The length of 'rn' must be less than 'N'.")

                      # Check that the numbers of defectives for each type are integers
                        if(any((object@N*object@pd < round(object@N*object@pd
                            )-1e-6)|(object@N*object@pd > round(object@N*object@pd
                            )+1e-6)))
                        return("N times pd must be integer numbers.")

                     ## fixed sampling plan
                        if(object@stype=="fixed"){
                          if(is.na(object@n) )
                          return("The 'n' is not allowed to contain missing value.")

                         if(length(object@n)!= 1) 
                          return("The length of sample size 'n' should be equal to 1")

                       # Check that the sample size is reasonable
                         if( object@n <= 0) 
                         return("Sample size 'n' should be greater than 0.")

                       # Check that the rejection numbers are reasonable         
                         if(any(object@rn-1 > object@n))
                         return("Any rejection number in 'rn' must be less than 'n+1'.")

                       # Check that the lot size and sample size are reasonable
                         if(object@N <= object@n) 
                         return("N must be greater than n.")
                        }

                     # Check that the length of acceptance is reasonable for 
                    ## sequential sampling plan

                         if(object@stype=="sequential")
                           {
                         if(is.na(object@m))
                          return("The 'm' is not allowed to contain missing value.")

                         if(length(object@m)!= 1) 
                          return("The length of cell quota 'm' should be equal to 1")

                       # Check that the cell quota of good item is reasonable
                         if( object@m <= 0) 
                         return("Cell quota 'm' should be greater than 0.")

                       # Check that the lot size and cell quota are reasonable
                         if(object@N <= object@m) 
                         return("N must be no less than m.")
                             
                           }         
                      
                     return(TRUE)
                     }
      )
##############################################################################################
##################################
# Creation of the object         #
#                                #
################################## 
Ocmult <- function(rn, type=c("multinomial", "hypergeom"), stype=c("fixed", "sequential"), ...)
       {
      # Choose what 'type' to use
        type <- match.arg(type)## Get type, defaut is multinomial
        stype <- match.arg(stype)## Get type, defaut is fixed
      # invoke a new object of that type
        switch(type,
               multinomial={
                          obj <- new("Ocmult.multinomial", rn=rn, type="multinomial", 
                                                            stype=stype, ...)
                          switch(stype,
                                 fixed={
                                        obj@pa <- pmultinom(x=obj@rn-1, size=obj@n, prob=obj@pd)
                                       },
                                 sequential={
                                           obj@pa <- calc.pnmultinom(pd=obj@pd, rn=obj@rn, m=obj@m)
                                           obj@asn <- EWT(pd=obj@pd, rn=obj@rn, m=obj@m)
                                          }
                                 ) 
                           },
               hypergeom={
                         obj <- new("Ocmult.hypergeom", rn=rn, type="hypergeom", 
                                    stype=stype, ...)
                         switch(stype,
                               fixed={
                                      obj@pa <- calc.pmultihyper(obj@pd, rn=obj@rn, 
                                                 n=obj@n, N=obj@N)
                                     },
                               sequential={
                                         obj@pa <- calc.pnmultihyper(pd=obj@pd, 
                                                    rn=obj@rn, m=obj@m, N=obj@N)
                                         obj@asn <- EWTH(pd=obj@pd, rn=obj@rn, m=obj@m, N=obj@N)
                                         }
                                )
                         }
               )
        return(obj)
       }
###############################################################################################
## Initialization of the class
setMethod("initialize", "Ocmult.multinomial",
          function(.Object, rn, n=.Object@n, m, pd, stype=.Object@stype, ...) 
                   {
                      lrn <- length(rn)
                      if(missing(pd))
                      {
                         if(lrn==2){
                           .Object@pd <- as.matrix(expand.grid(seq(0, 0.1, 0.01), seq(0, 0.1, 0.01)))
                           }
                       else    
                         if(lrn > 2)
                           {
                              de <- 1/25/(lrn+1)/lrn*c(1:lrn)
                              x <- matrix(0,nrow=10,ncol=lrn)
                              for (i in c(1:10)) x[i,] <- i*de
                              .Object@pd <- x
                            }
                       }
                    else if(is.vector(pd))
                      {
                        if(lrn > 1)
                          {
                          l.row <- length(pd)/lrn
                          if(l.row < 1)
                          stop("The length of 'rn' must less than or equal to length of 'pd'.")
                          if(l.row > floor(l.row))
                             {
                              l.row <- floor(l.row)
                              l.pd <- lrn*l.row
                              pd <- pd[1:l.pd]                           
                          warning( 
                            paste("The length of the pd vector should be an integer",
                                "multiple of the length of the rn vector.",
                                "\nThe truncated pd in use is: "), t(paste(pd, " ")))          
                             }
                          }
                       if(lrn==2){
                          pdtemp <- matrix(pd, ncol=lrn, byrow=TRUE)
                          .Object@pd <- as.matrix(expand.grid(unique(sort(pdtemp[, 1])), 
                                        unique(sort(pdtemp[, 2]))))
                          }
                       else
                         .Object@pd <- matrix(pd, ncol=lrn, byrow=TRUE)
                      } 
                    else
                      if(lrn==2)
                       .Object@pd <- as.matrix(expand.grid(unique(sort(pd[, 1])), 
                                        unique(sort(pd[, 2]))))
                     else
                      .Object@pd <- pd 
                      if(stype=="fixed")
                        .Object@n <- n 
                      else 
                        .Object@m <- m 
                     .Object@rn <- rn 
                     .Object@stype <- stype          
                      callNextMethod(.Object, ...)## Return to object
                     }
         )

setMethod("initialize", "Ocmult.hypergeom",
          function(.Object, rn, n=.Object@n, m, pd, N=.Object@N, stype=.Object@stype, ...) 
                   { 
                    lrn <- length(rn)
                    if(missing(pd))
                      {
                         if(lrn==2){
                            .Object@pd <- as.matrix(expand.grid(seq(0, 0.1, 0.01), seq(0, 0.1, 0.01)))
                           }
                       else    
                         if(lrn > 2 )
                           {
                              de <- 1/25/(lrn+1)/lrn*c(1:lrn)
                              x <- matrix(0,nrow=10,ncol=lrn)
                              for (i in c(1:10)) x[i,] <- i*de
                              .Object@pd <- x
                             .Object@pd <- round(N*.Object@pd)/N
                            }
                       }
                    else if(is.vector(pd))
                      {
                        if(lrn > 1)
                          {
                          l.row <- length(pd)/lrn
                          l.row <- length(pd)/lrn
                          if(l.row < 1)
                          stop("The length of 'rn' must less than or equal to length of 'pd'.")
                          if(l.row > floor(l.row))
                             {
                              l.row <- floor(l.row)
                              l.pd <- lrn*l.row
                              pd <- pd[1:l.pd]                           
                          warning( 
                            paste("The length of the pd vector should be an integer",
                                "\nmultiple of the length of the rn vector.",
                                "\nThe truncated pd in use is: \n"), t(paste(pd, " ")))           
                             }
                          }
                          if(lrn==2){
                          pdtemp <- matrix(pd, ncol=lrn, byrow=TRUE)
                          .Object@pd <- as.matrix(expand.grid(unique(sort(pdtemp[, 1])), 
                                        unique(sort(pdtemp[, 2]))))
                          }
                       else
                      .Object@pd <- matrix(pd, ncol=lrn, byrow=TRUE)
                      } 
                    else 
                    if(lrn==2)
                       .Object@pd <- as.matrix(expand.grid(unique(sort(pd[, 1])), 
                                        unique(sort(pd[, 2]))))
                    else
                   .Object@pd <- pd 
                    if(stype=="fixed")
                      .Object@n <- n 
                    else 
                      .Object@m <- m 
                   .Object@rn <- rn                  
                   .Object@N <- N
                   .Object@stype <- stype       
                    callNextMethod(.Object, ...)
                   }
        )
###############################################################################


##########################################################################################

## Create show function
setMethod ("show" , "Ocmult",
           function(object)
              {
             if(object@stype=="fixed")
                   {
                    switch(object@type,
                            multinomial=cat(" ", paste(length(object@rn)+1,
                              "-Level Acceptance Sampling Plan Multinomial:", sep=""), "\n"),
                            hypergeom=cat(" ", paste(length(object@rn)+1,
                                "-Level Acceptance Sampling Plan Multivariate Hypergeom: N = ",
                                 object@N, sep=""), "\n" )
                           )
                         cat(" Sample size:   ", paste(object@n), "\n" )

                         if(!is.null(names(object@rn)))
                            {
                            rejnum <- as.matrix(t(object@rn))
                            colnames(rejnum) <- c(names(object@rn))
                            rownames(rejnum) <- c(" Rej. Number(s):")
	                    show(rejnum)
                            }
                         else 
                            {
                            cat(" Rej. Number(s):", paste(object@rn), "\n" )
                            }
                         cat("\n" )
                    }
               else {
                         switch(object@type,
                                multinomial=cat(" ", paste(length(object@rn)+1,
                                "-Level Sequential Acceptance Sampling Plan",
                                " Negative Multinomial:", sep=""), "\n"),
                                 hypergeom=cat(" ", paste(length(object@rn)+1,
                                  "-Level Sequential Acceptance Sampling Plan",
                                  " Negative Multivariate Hypergeom: N = ",
                                 object@N, sep=""), "\n")
                                )
                        cat(" Acc. Number:   ", paste(object@m), "\n" )

                        if(!is.null(names(object@rn))) 
                            {
                            rejnum <- as.matrix(t(object@rn))
                            colnames(rejnum) <- c(names(object@rn))
                            rownames(rejnum) <- c(" Rej. Number(s):")
	                    show(rejnum)
                            }
                         else 
                            {
                            cat(" Rej. Number(s):", paste(object@rn), "\n" )
                            }

                         cat("\n" )
                     }
                  
             }
           )


## Create summary function
setMethod("summary", signature(object="Ocmult"),
          function(object, detail=FALSE)
          {
              if(ncol(object@pd)==1) p.nondef <- 1- object@pd              
              else p.nondef <- 1-rowSums(object@pd)
              pa <- round(object@pa,7)
              pd <- round(object@pd,4)
              p.nondef <- round(p.nondef,4)
              prop <- cbind(pd, p.nondef ,pa)
              l=ncol(object@pd)          

              if(is.null(names(object@rn)))
              defname = as.vector(sapply(" type",  FUN = paste, (1:l)))
              else defname = names(rn)

              rownames(prop) <- rep(" ", length(object@pa))
              if(object@stype=="fixed")
              colnames(prop) <- c(defname, " P.nondef", " P(accept)")
              else
              {
               asn <- round(object@asn,7)
               prop <- cbind(prop, asn)
               colnames(prop) <- c(defname, " P.nondef", " P(accept)", "       ASN") 
              }
             show(object) 
             if(detail){
                      cat(" Detailed acceptance probabilities: \n")
                      show(prop)
                     } 
           if(object@stype=="fixed")
             return(invisible(c(list(n=object@n,  rn=object@rn, P= prop))))
           else
             return(invisible(c(list(m=object@m, rn=object@rn, P= prop))))
           
          }
      )
###########################################################################

## Creating a new generic function for plot

setMethod("plot", c(x = "Ocmult.multinomial", y = "missing"),
          function(x, y, type="o", xlab="Proportion Defective", ylab="P(accept)", 
                   main = main.2dp(x), ...) 
           {
             if(length(x@rn)!=1) 
             stop("The plot for 2-Level acceptance sampling plan only")
             plot(x@pd, x@pa, type=type, xlab=xlab, ylab=ylab, main=main,...)
           }
         )

setMethod("plot", signature(x="numeric", y="Ocmult.multinomial"),
          function(x, y, type="o", ylab="P(accept)", main = main.2dp(y),  ...)
         {
            plot(x, y@pa, type=type, ylab=ylab,  main=main, ...)
          }
         )

setMethod("plot", c(x = "Ocmult.hypergeom", y = "missing"),
          function(x, y, type="p", xlab="Proportion Defective", ylab="P(accept)", 
                   main = main.2dp(x), ...) 
           {
             if(length(x@rn)!=1)
             stop("The plot for 2-Level acceptance sampling plan only")
             plot(x@pd, x@pa, type=type, xlab=xlab, ylab=ylab, main=main,...)
            }
          )

setMethod("plot", signature(x="numeric", y="Ocmult.hypergeom"),
          function(x, y, type="p", ylab="P(accept)", main = main.2dp(y), ...)
         {
            plot(x, y@pa, type=type, ylab=ylab, main=main, ...)
          }
         )

setMethod("persp", c(x = "Ocmult.multinomial"),
          function(x, y,  zlab="P(accept)", xlab="p1", ylab="p2",
                   theta = 20, phi = 30,d=4, expand=0.5, 
                   main = main.3dp(x), ticktype = "detailed",...) 
           {
             if(length(x@rn)!=2)
             stop("The persp plot is only for 3-Level acceptance sampling.")

             p1 <- unique(sort(x@pd[,1]))
             p2 <- unique(sort(x@pd[,2]))
             pa <- matrix(x@pa, nrow=length(p1))
           ## Create three dimension plot   
             persp(p1, p2, pa, theta = theta, phi = phi,d=d, expand=expand,
                   zlab=zlab, xlab=xlab, ylab=ylab, main=main, ticktype = ticktype,...)
           }
         )

setMethod("persp", c(x = "Ocmult.hypergeom"),
          function(x, y,  zlab="P(accept)", xlab="p1", ylab="p2", 
                    theta = 20, phi = 30,d=4, expand=0.5,
                   main = main.3dp(x), ticktype = "detailed",...) 
           {
             if(length(x@rn)!=2)
             stop("The persp plot is only for 3-Level acceptance sampling.")

             p1 <- unique(sort(x@pd[,1]))
             p2 <- unique(sort(x@pd[,2]))
             
             pa <- matrix(x@pa, nrow=length(p1))   
           ## Create three dimension plot   
             persp(p1, p2, pa, theta = theta, phi = phi,d=d, expand=expand,
                   zlab=zlab, xlab=xlab, ylab=ylab, main=main, 
                   ticktype = ticktype,...)
           }
         )

## Creating a new generic function for contour plot 
setMethod("contour", c(x = "Ocmult.multinomial"),
          function(x, y,  nlevel=8, main = main.3dc(x), xlab="p1", ylab="p2",...) 
           {
             if(length(x@rn)!=2)
             stop("The contour plot is only for 3-Level acceptance sampling.")

             p1 <- unique(sort(x@pd[,1]))
             p2 <- unique(sort(x@pd[,2]))

             pa <- matrix(x@pa, nrow=length(p1))               
             contour(p1, p2, pa, nlevel=nlevel, main=main, xlab=xlab, ylab=ylab,...)
           }
         )

setMethod("contour", c(x = "Ocmult.hypergeom"),
          function(x, y,  nlevel=8,  main = main.3dc(x), xlab="p1", ylab="p2",...) 
           {
             if(length(x@rn)!=2)
             stop("The contour plot is only for 3-Level acceptance sampling.")

             p1 <- unique(sort(x@pd[,1]))
             p2 <- unique(sort(x@pd[,2]))
             pa <- matrix(x@pa, nrow=length(p1))        
             contour(p1, p2, pa, nlevel=nlevel, main=main, xlab=xlab, ylab=ylab,...)
           }
         )
#########################################################################
##The function for plot title
main.3dp <- function(obj)
           {
            if(obj@type=="multinomial")
               {
                if(obj@stype=="fixed")
                   return(paste("Multinomial OC Surface with  \nn = ", obj@n, 
                  ", rn = (", obj@rn[1], ",", obj@rn[2], ")", sep=""))
                 else 
                   return(paste("Negative Multinomial OC Surface with \nm = ", 
                          obj@m, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
                }
            else
               {
                if(obj@stype=="fixed")
                   return(paste("Multivariate Hypergeometric OC Surface with \nn = ",
                                             obj@n, ", N = ", obj@N, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
                 else 
                   return(paste("Negative Multivariate Hypergeometric OC Surface with \nm = ", obj@m, 
                   ", N = ", obj@N, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
               }
               
            }

main.3dc <- function(obj)
           {
            if(obj@type=="multinomial")
               {
                if(obj@stype=="fixed")
                   return(paste("Multinomial OC Contour with  \nn = ", 
                          obj@n, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
                 else 
                   return(paste("Negative Multinomial OC Contour with \nm = ", 
                          obj@m, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
                }
            else
               {
                if(obj@stype=="fixed")
                   return(paste("Multivariate Hypergeometric OC Contour with \nn = ", obj@n, 
                    ", N = ", obj@N, ", rn = (", obj@rn[1], ",", obj@rn[2], ")", sep=""))
                 else 
                   return(paste("Negative Multivariate Hypergeometric OC Contour with \nm = ", obj@m,
                              ", N = ", obj@N, ", rn = (", obj@rn[1], ",", obj@rn[2],")", sep=""))
               }
               
            }


 main.2dp <- function(obj)
           {
            if(obj@type=="multinomial")
               {
                if(obj@stype=="fixed")
                   return(paste("Binomial OC Curve with \nn = ", 
                          obj@n, ", rn = ", obj@rn, sep=""))
                 else 
                   return(paste("Negative Binomial OC Curve with \nm = ", 
                          obj@m, ", rn = ", obj@rn[1], sep=""))
                }
            else
               {
                if(obj@stype=="fixed")
                   return(paste("Hypergeometric OC Curve with \nn = ", obj@n, 
                          ", N = ", obj@N, ", rn = ", obj@rn, sep=""))
                 else 
                   return(paste("Negative Hypergeometric OC Curve with \nm = ", obj@m,
                               
              ", N = ", obj@N, ", rn = ", obj@rn[1], sep=""))
               }
           }


###########################################################################################


