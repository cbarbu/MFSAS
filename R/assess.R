###################################################################
# Assessment function is used to assess whether the given plan    #
# can meet the criteria specified in PRP (Producer's Risk Point)  #
# and/or CRP (Consumer's Risk Point).                             #
###################################################################


assess.multi <-	function(rn, n=30, m, N=100, PRP, CRP, type=c("multinomial", "hypergeom"),                          
           stype=c("fixed", "sequential"), print=TRUE)
         {
          if(any(rn < 1))
            stop("The values in 'rn' must not be less than 1.")
          if(missing(PRP) & missing(CRP))
              stop("At least one risk point, PRP or CRP, must be specified")

         type <- match.arg(type)
         stype <- match.arg(stype)
          
         if(stype=="sequential"){
            if(missing(m))
              stop("'m' is missing.")
            if(m < 1)
             stop("'m' must be greater than 0.")
            }
  
         if(missing(PRP))
             {
              l <- length(CRP)
              if(any(CRP < 0)|any(CRP > 1))
                 stop("'CRP' is out of range (0,1).")
              if(l>2)
                 {
                  if(sum(CRP[-l])>1)
                     stop("Sum of risk point must not be greater than 1.")
                  }
              PRP <- rep(NA,l)
              if(stype=="fixed")
                {
                 if(type=="multinomial")
                   pcons <- pmultinom(rn-1, n, CRP[-l])
                 else
                   pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
                 }
               else
                 {
                  if(type=="multinomial")
                    {
                     pcons <- calc.pnmultinom(CRP[-l], rn, m)
                     asncons <- EWT(CRP[-l], rn, m)
                    }
                   else
                     {
                      pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
                      asncons <- EWTH(CRP[-l], rn, m, N)
                     }
                 }
              if(stype=="fixed")
                 result2 <- c(CRP, pcons )
               else
                 result2 <- c(CRP, pcons, asncons)
  
              if(pcons >= CRP[l])
                 plan.meet <- FALSE
              else 
                 plan.meet <- TRUE 

                result <- as.matrix(t(result2))

              if(is.null(names(rn)))
                defname = as.vector(sapply(" type",  FUN = paste, (1:(l-1))))
              else defname = names(rn)

                if(stype=="fixed")
                dimnames(result) <- list("CRP",
                                      c(defname,
                                      " RP P(accept)", " Plan P(accept)"))
               else
                 dimnames(result) <- list("CRP",
                                      c(defname, " RP P(accept)", 
                                        " Plan P(accept)", "       ASN"))

   
             }
          else if(missing(CRP))
             {
              l <- length(PRP)
              if(any(PRP < 0)| any(PRP[l]>1))
                 stop("'PRP' is out of range (0,1).")
              if(l>2)
                {
                 if(sum(PRP[-l])>1)
                    stop("Sum of risk point must not be greater than 1.")
                }
              CRP <- rep(NA,l)
              if(stype=="fixed")
                {
                 if(type=="multinomial")
                    pprod <- pmultinom(rn-1, n, PRP[-l])
                  else
                    pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
                 }
              else
                 {
                  if(type=="multinomial")
                     {
                      pprod <- calc.pnmultinom(PRP[-l], rn, m)
                      asnprod <- EWT(PRP[-l], rn, m)
                      }
                   else
                      {
                       pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
                       asnprod <- EWTH(PRP[-l], rn, m, N)
                      }
                 }
              if(stype=="fixed")
                 result1 <- c(PRP, pprod)
              else
                  result1 <- c(PRP, pprod, asnprod)
 
              if(pprod <= PRP[l])
                 plan.meet <- FALSE
              else 
                 plan.meet <- TRUE 

                result <- as.matrix(t(result1))

              if(is.null(names(rn)))
                defname = as.vector(sapply(" type",  FUN = paste, (1:(l-1))))
              else defname = names(rn)

                if(stype=="fixed")
                dimnames(result) <- list("PRP",
                                      c(defname,
                                      " RP P(accept)", " Plan P(accept)"))
               else
                 dimnames(result) <- list("PRP",
                                      c(defname, " RP P(accept)", 
                                        " Plan P(accept)", "       ASN"))
  
             }
         else
            {
             l <- length(PRP)
             if(any(CRP[-l] <= PRP[-l]))
                stop(
                    "Consumer Risk Point quality must be greater than Producer Risk Point quality.")
             if(any(PRP < 0)| any(PRP > 1))
                stop("'PRP' is out of range (0,1).")
             if(any(CRP < 0)|any( CRP > 1))
                stop("'CRP' is out of range (0,1).")
             if(l>2){if(sum(PRP[-l])>1|sum(CRP[-l])>1)stop(
                     "Sum of risk point must not be greater than 1.")}

             if(stype=="fixed")
                {
                 if(type=="multinomial")
                    {
                     pprod <- pmultinom(rn-1, n, PRP[-l])
                     pcons <- pmultinom(rn-1, n, CRP[-l])
                     }
                  else
                     {
                      pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
                      pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
                      }
                   }
               else
                 {
                  if(type=="multinomial")
                     {
                      pprod <- calc.pnmultinom(PRP[-l], rn, m)
                      pcons <- calc.pnmultinom(CRP[-l], rn, m)
                      asnprod <- EWT(PRP[-l], rn, m)
                      asncons <- EWT(CRP[-l], rn, m)                                
                      }
                  else
                      {
                       pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
                       pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
                       asnprod <- EWTH(PRP[-l], rn, m, N)
                       asncons <- EWTH(CRP[-l], rn, m, N)
                      }
                  }
             if(stype=="fixed")
                {        
                 result1 <- c(PRP, pprod)
                 result2 <- c(CRP, pcons )
                }
             else
                {
                 result1 <- c(PRP, pprod, asnprod)
                 result2 <- c(CRP, pcons, asncons)
                }
             if(pprod <= PRP[l]| pcons >= CRP[l])
                plan.meet <- FALSE
             else 
                plan.meet <- TRUE 
            
                result <- rbind(result1, result2)

              if(is.null(names(rn)))
                defname = as.vector(sapply(" type",  FUN = paste, (1:(l-1))))
              else defname = names(rn)

                if(stype=="fixed")
                dimnames(result) <- list(c("PRP", "CRP" ),
                                      c(defname,
                                      " RP P(accept)", " Plan P(accept)"))
               else
                 dimnames(result) <- list(c("PRP", "CRP" ),
                                      c(defname, " RP P(accept)", 
                                        " Plan P(accept)", "       ASN"))
           }
          if(print)
            {
            if(stype=="fixed")
                {
                 if(type=="multinomial")
                    {
                      cat(" ", paste(l,
                       "-Level Acceptance Sampling Plan Multinomial:", sep=""), "\n")
                        cat(" Sample size:   ", paste(n), "\n" )

                      if(!is.null(names(rn)))
                         {
                         rejnum <- as.matrix(t(rn))
                         colnames(rejnum) <- c(names(rn))
                         rownames(rejnum) <- c(" Rej. Number(s):")
	                 show(rejnum)
                         }
                      else 
                         {
                         cat(" Rej. Number(s):", paste(rn), "\n" )
                         }

                      cat("\n" )
                     }
                  else
                    {
                     cat(" ", paste(l,
                     "-Level Acceptance Sampling Plan Multivariate Hypergeometric: N =",
                      N, sep=""), "\n")
                      cat(" Sample size:   ", paste(n), "\n" )

                      if(!is.null(names(rn)))
                         {
                         rejnum <- as.matrix(t(rn))
                         colnames(rejnum) <- c(names(rn))
                         rownames(rejnum) <- c(" Rej. Number(s):")
	                 show(rejnum)
                         }
                      else 
                         {
                         cat(" Rej. Number(s):", paste(rn), "\n" )
                         }

                      cat("\n" )
                     }
                }
              else
                 {
                   if(type=="multinomial")
                     {
                      cat(" ", paste(l,
                         "-Level Acceptance Sampling Plan Negative Multinomial:", sep=""), "\n")
                      cat(" Acc. Number:   ", paste(m), "\n" )

                      if(!is.null(names(rn)))
                         {
                         rejnum <- as.matrix(t(rn))
                         colnames(rejnum) <- c(names(rn))
                         rownames(rejnum) <- c(" Rej. Number(s):")
	                 show(rejnum)
                         }
                      else 
                         {
                         cat(" Rej. Number(s):", paste(rn), "\n" )
                         }

                      cat("\n" )
                     }
                   else
                     {
                      cat(" ", paste(l,
                         "-Level Acceptance Sampling Plan Negative Multivariate Hypergeometric: N =",
                         N, sep=""), "\n")
                      cat(" Acc. Number:   ", paste(m), "\n" )

                      if(!is.null(names(rn)))
                         {
                         rejnum <- as.matrix(t(rn))
                         colnames(rejnum) <- c(names(rn))
                         rownames(rejnum) <- c(" Rej. Number(s):")
	                 show(rejnum)
                         }
                      else 
                         {
                         cat(" Rej. Number(s):", paste(rn), "\n" )
                         }

                      cat("\n" )
                     }
                 }             

            if(plan.meet)
               {
                cat(" Plan CAN meet desired risk point(s): \n")}
            else 
               {
                cat(" Plan CANNOT meet desired risk point(s): \n")
               }
##            print(formatC(result, digits = 8, format = "f", drop0trailing=TRUE),quote = FALSE)
                show(result)
           }
           if(stype=="fixed")
                {
                 if(type=="multinomial")
                    {
                     return(invisible(c(list(n=n, rn=rn,
                      result=result, plan.meet=plan.meet))))
                    }
                  else
                    {
                     return(invisible(c(list(N=N, n=n, rn=rn,
                     result=result, plan.meet=plan.meet))))
                    }
                  }
                else
                  {
                   if(type=="multinomial")
                    {
                     return(invisible(c(list(rn=rn, m=m,
                     result=result, plan.meet=plan.meet))))
                    }
                   else
                    {
                     return(invisible(c(list(N=N, rn=rn, m=m,
                     result=result, plan.meet=plan.meet))))
                    }
                  }      
         } 

