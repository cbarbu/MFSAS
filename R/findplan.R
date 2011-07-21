find.multi.plan <- function(PRP, CRP, N=100, type= c("multinomial", "hypergeom"),
                       stype= c("fixed", "sequential"))
          {
            type <- match.arg(type)
            stype <- match.arg(stype)
            l=length(PRP)
            k=l-1
            done=0
            if(stype=="fixed")
              {
               if(type=="multinomial")
                 {
                  n <- 0
                   while(done==0)
                    {
                     n=n+1
                     rn=rep(1, k)
                     pcons <- pmultinom(rn-1, n, CRP[-l])
                    if(pcons <= CRP[l])
                     {
                      pprod <- pmultinom(rn-1, n, PRP[-l])
                      if(pprod >= PRP[l])
                            {
                             done=1
                             break
                            }
                          else
                            {
                             func=n.plan(k, rn, n, PRP, CRP, pcons)
                             done=func$done
                             if(done==1)
                              {
                               pprod=func$pprod
                               pcons=func$pcons
                               rn=func$rn
                              }
                            }
                      }
                    }
                  }
                else
                  {
                   n <- 0
                   while(done==0)
                     {
                      n=n+1
                      rn=rep(1, k)
                      pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
                      if(pcons <= CRP[l])
                        {
                         pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
                         if(pprod >= PRP[l])
                            {
                             done=1
                             break
                            }
                          else
                            {
                             func=h.plan(k, rn, n, N, PRP, CRP, pcons)
                             done=func$done
                             if(done==1)
                               {
                                pprod=func$pprod
                                pcons=func$pcons
                                rn=func$rn
                               }
                            }
                         }
                       }
                     }
                 }
               else
                 {
                 if(type=="multinomial")
                 {
                  m <- 0
                  while(done==0)
                     {
                      m=m+1
                      rn=rep(1, k)
                      pprod <- calc.pnmultinom(PRP[-l], rn, m)
                      pcons <- calc.pnmultinom(CRP[-l], rn, m)
                      if(pcons <= CRP[l] & pprod >= PRP[l])
                        {
                         done=1
                         break
                         }
                      else
                         {
                         if(pcons <= CRP[l])
                           {
                            func=nm.plan(k, rn, m, PRP, CRP, pcons)
                            done=func$done
                            if(done==1)
                              {
                               pprod=func$pprod
                               pcons=func$pcons
                               rn=func$rn
                               m=func$m
                              }
                            }
                          }
                        }
                        asnprod <- EWT(PRP[-l], rn, m)
                        asncons <- EWT(CRP[-l], rn, m)
                      }
                  else
                     {
                     m <- 0
                     while(done==0)
                        {
                         m=m+1
                         rn=rep(1, k)
                         pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
                         if(pcons <= CRP[l])
                           {
                             pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
                             if( pprod >= PRP[l])
                               {
                               done=1
                               break
                               }
                             else
                              {
                                func=nh.plan(k, rn, m, PRP, CRP, N, pcons)
                                done=func$done
                               if(done==1)
                                 {
                                  pprod=func$pprod
                                  pcons=func$pcons
                                  rn=func$rn
                                  m=func$m
                                  }
                               }
                            }
                          }
                       asnprod <- EWTH(PRP[-l],  func$rn, func$m, N)
                       asncons <- EWTH(CRP[-l],  func$rn, func$m, N)
                    }
                 }

               cat(" The optimal plan is: \n\n")

               if(stype=="fixed")
	{

                 if(type=="multinomial")
                    {
                      cat(" ", paste(l,
                       "-Level Acceptance Sampling Plan Multinomial:", sep=""), "\n")
                        cat(" Sample size:   ", paste(func$n), "\n" )

                         cat(" Rej. Number(s):", paste(func$rn), "\n" )

                      cat("\n" )
                     }
                  else
                    {
                     cat(" ", paste(l,
                     "-Level Acceptance Sampling Plan Multivariate Hypergeometric: N = ",
                      N, sep=""), "\n")
                      cat(" Sample size:   ", paste(func$n), "\n" )

                         cat(" Rej. Number(s):", paste(func$rn), "\n" )

                      cat("\n" )
                     }


	   result1 <- c(PRP,round(func$pprod,7))
	   result2 <- c(CRP,round(func$pcons,7))
	   result <- rbind(result1,result2)
                    defname = as.vector(sapply(" type",  FUN = paste, (1:(l-1))))
                    dimnames(result) <- list(c("PRP", "CRP" ),
                                      c(defname, " RP P(accept)", 
                                        " Plan P(accept)"))
                    show(result)
                     return(invisible(list(n=func$n, rn=func$rn, p.PRP=func$pprod, p.CRP=func$pcons)))

	}
               else
	{

                   if(type=="multinomial")
                     {
                      cat(" ", paste(l,
                         "-Level Acceptance Sampling Plan Negative Multinomial:", sep=""), "\n")
                      cat(" Acc. Number:   ", paste(m), "\n" )

                         cat(" Rej. Number(s):", paste(rn), "\n" )

                      cat("\n" )
                     }
                   else
                     {
                      cat(" ", paste(l,
                         "-Level Acceptance Sampling Plan Negative Multivariate Hypergeometric: N = ",
                         N, sep=""), "\n")
                      cat(" Acc. Number:   ", paste(m), "\n" )

                         cat(" Rej. Number(s):", paste(rn), "\n" )

                      cat("\n" )
                     }

	   result1 <- c(PRP,round(pprod,7),round(asnprod,7))
	   result2 <- c(CRP,round(pcons,7),round(asncons,7))
	   result <- rbind(result1,result2)
                    defname = as.vector(sapply(" type",  FUN = paste, (1:(l-1))))
                    dimnames(result) <- list(c("PRP", "CRP" ),
                                      c(defname, " RP P(accept)", 
                                        " Plan P(accept)", "       ASN"))
                 show(result)
                 return(invisible(list(m=m, rn=rn, p.PRP=pprod, p.CRP=pcons,
                                ASNp=asnprod, ASNc=asncons)))


	}



          }

## recursive find plan function for multinomial distribution
n.plan <- function(b, rn, n, PRP, CRP, pcons)
        {
         l=length(PRP)
         k=l-1
         done=0
         if(b==1)
           {
            while(pcons <= CRP[l])
              {
               pprod <- pmultinom(rn-1, n, PRP[-l])
               pcons <- pmultinom(rn-1, n, CRP[-l])
               if(pcons <= CRP[l] & pprod >= PRP[l])
                 {
                  done=1
                  break
                  }
                if(sum(rn-1)>= n)
                  break
                else
                  rn[k]=rn[k]+1
                }
                if(done==1)
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                else
                 {
                  if(k!=b)
                     {
                         rn[k]=1
                         rn[k-1]=rn[k-1]+1
                         pprod <- pmultinom(rn-1, n, PRP[-l])
                         pcons <- pmultinom(rn-1, n, CRP[-l])
                    }
                 }
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                }

                else
                while(pcons <= CRP[l])
                 {
                  fuc=n.plan(b-1, rn, n, PRP, CRP, pcons)
                     done = fuc$done
                     pprod=fuc$pprod
                     pcons=fuc$pcons
                   if(done==1)
                     {
                     rn=fuc$rn
                     break
                     }
                    else
                     {
                      rn[k-b+1] <-rn[k-b+1]+1
                      rn[(k-b+2):k] <- 1
                      if(pcons <= CRP[l] & pprod >= PRP[l])
                      {
                      done=1
                      break
                      }
                     }
                   if(sum(rn-1)>= n)
                     break
                    }
                if(done==1)
                return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                else
                  {
                   if(k!=b)
                    {
                    rn[(k-b+1):k] <-1
                    rn[k-b]=rn[k-b]+1
                     pprod <- pmultinom(rn-1, n, PRP[-l])
                     pcons <- pmultinom(rn-1, n, CRP[-l])
                    }
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                   }
        }

## recursive find plan function for multivariate hypergeometric distribution
h.plan <- function(b, rn, n, N, PRP, CRP, pcons)
       {
         l=length(PRP)
         k=l-1
         done=0
          if(b==1)
           {
            while(pcons <= CRP[l])
              {
               pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
               pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
               if(pcons <= CRP[l] & pprod >= PRP[l])
                 {
                  done=1
                  break
                  }
                if(sum(rn-1)>= n)
                  break
                else
                  rn[k]=rn[k]+1
                }
                if(done==1)
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                else
                 {
                  if(k!=b)
                     {
                         rn[k]=1
                         rn[k-1]=rn[k-1]+1
                         pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
                         pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
                    }
                 }
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                }
                else
                while(pcons <= CRP[l])
                 {
                  fuc=fuc=h.plan(b-1, rn, n, N, PRP, CRP, pcons)
                     done = fuc$done
                     pprod=fuc$pprod
                     pcons=fuc$pcons
                   if(done==1)
                     {
                     rn=fuc$rn
                     break
                     }
                    else
                     {
                      rn[k-b+1] <-rn[k-b+1]+1
                      rn[(k-b+2):k] <- 1
                      if(pcons <= CRP[l] & pprod >= PRP[l])
                      {
                      done=1
                      break
                      }
                     }
                   if(sum(rn-1)>= n)
                     break
                    }
                if(done==1)
                return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                else
                  {
                   if(k!=b)
                    {
                    rn[(k-b+1):k] <-1
                    rn[k-b]=rn[k-b]+1
                    pprod <- calc.pmultihyper(PRP[-l], rn, n, N)
                    pcons <- calc.pmultihyper(CRP[-l], rn, n, N)
                    }
                 return(list(n=n, rn=rn, pprod=pprod, pcons=pcons, done=done))
                   }
        }

## recursive find plan function for negative multinomial distribution

nm.plan <- function(b, rn, m, PRP, CRP, pcons)
        {
         l=length(PRP)
         k=l-1
         done=0
         if(b==1)
           {
            while(pcons <= CRP[l]) 
              {
               pprod <- calc.pnmultinom(PRP[-l], rn, m)
               pcons <- calc.pnmultinom(CRP[-l], rn, m)
               if(pcons <= CRP[l] & pprod >= PRP[l])
                  {
                  done=1 
                  break
                  }
                if(rn[k] > m)
                 break
                else
                 rn[k]=rn[k]+1
               }
            if(done==1)
                 return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))
            else
               {
                if(k!=b)
                    {
                     rn[k]=1
                     rn[k-1]=rn[k-1]+1
                     pprod <- calc.pnmultinom(PRP[-l], rn, m)
                     pcons <- calc.pnmultinom(CRP[-l], rn, m)
                    }
                }
                 return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))              
                }
           else
                while(pcons <= CRP[l])
                 {
                  fuc=nm.plan(b-1, rn, m, PRP, CRP, pcons)
                     done = fuc$done
                     pprod=fuc$pprod
                     pcons=fuc$pcons
                   if(done==1)
                     {
                     rn=fuc$rn
                     break
                     }
                    else
                     {
                      rn[k-b+1] <-rn[k-b+1]+1
                      rn[(k-b+2):k] <- 1
                      if(pcons <= CRP[l] & pprod >= PRP[l])
                      {
                      done=1 
                      break
                      }                      
                     }
                   if(rn[k-b+1] > m)
                     break
                    }
                if(done==1)
                return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done)) 
                else
                  {
                   if(k!=b)
                    {
                    rn[(k-b+1):k] <-1
                    rn[k-b]=rn[k-b]+1
                    pprod <- calc.pnmultinom(PRP[-l], rn, m)
                    pcons <- calc.pnmultinom(CRP[-l], rn, m)
                    }
                 return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))
                   }            
        }

## recursive find plan function for negative multivariate hypergeometric distribution
nh.plan <- function(b, rn, m, PRP, CRP, N, pcons)
        {
         l=length(PRP)
         k=l-1
         done=0
         if(b==1)
           {
            while(pcons <= CRP[l]) 
              {
               pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
               pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
               if(pcons <= CRP[l] & pprod >= PRP[l])
                 {
                  done=1 
                  break
                  }
                if(rn[k] > m)
                 break
                else
                 rn[k]=rn[k]+1
                }
                if(done==1)
                 return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))
               else
                 {
                  if(k!=b)
                    {
                     rn[k]=1
                     rn[k-1]=rn[k-1]+1
                     pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
                     pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
                    }
                  }
                 return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))              
                }
           else
                while(pcons <= CRP[l])
                 {
                  fuc=nh.plan(b-1, rn, m, PRP, CRP, N, pcons)
                     pprod=fuc$pprod
                     pcons=fuc$pcons
                     done = fuc$done
                   if(done==1)
                     {
                     rn=fuc$rn
                     break
                     }
                    else
                     {
                      rn[k-b+1] <-rn[k-b+1]+1
                      rn[(k-b+2):k] <- 1
                      if(pcons <= CRP[l] & pprod >= PRP[l])
                      {
                      done=1 
                      break
                      }                      
                     }
                   if(rn[k-b+1] > m)
                     break
                    }
                if(done==1)
                return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done)) 
                else
                  {
                  if(k!=b)
                   {
                    rn[(k-b+1):k] <- 1
                    rn[k-b]=rn[k-b]+1
                    pprod <- calc.pnmultihyper(PRP[-l], rn, m, N)
                    pcons <- calc.pnmultihyper(CRP[-l], rn, m, N)
                   }
                return(list(m=m, rn=rn, pprod=pprod, pcons=pcons, done=done))
                   }            
        }
