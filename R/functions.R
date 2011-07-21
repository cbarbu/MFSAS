###################################################################
## Calculate the pa (the acceptance probabilities) functions.     #
# They are used in the initialzation function when creating       #
# the object.                                                     #
# calc.pmultihyper is used for fixed sampling with multivariate   #
# hypergeometric distribution.                                    #
# calc.pnmultinom is used for sequential sampling with negative   #
# multinomial distribution.                                       #
# calc.pnmultihyper is used for sequential sampling with negative #
# multivariate hypergeometric distribution.                       #
###################################################################        
calc.pmultihyper <- function(pd, rn, n, N)
                  { 
                   M <- round(N*pd)
                   if(((N*pd) < (M-1e-6))||((N*pd) > (M+1e-6)))
                    stop("'N' times 'pd' must be integer numbers.")
                   pa <- pmultihyper(rn-1, n, M, N)
                   return(pa)
                  }
calc.pnmultinom <- function(pd, rn, m)
             {
               R <- rn-1
               pa <- pnmultinom(R, m, prob=pd)
               return(pa)
              }
calc.pnmultihyper <- function(pd, rn, m, N)
             {
              M=round(N*pd)
              if(((N*pd) < (M-1e-6))||((N*pd) > (M+1e-6)))
                    stop("'N' times 'pd' must be integer numbers.")
              R=rn-1                           
              pa <- pnmultihyper(R, m, M, N)                               
              return(pa)
             }

#################################################################
# Calculate the expected sample size (average sampling number)  #
# using the D function for the arguments matching the ones in   #
# the classes.                                                  #
# pd is a vector or matix.                                      #
# rn is a nonnegtive vector.                                    #
#################################################################
EWT <- function(pd, rn, m)
             {
               b=length(rn)
               R=rn
               if(is.matrix(pd))
                  {
                   if(ncol(pd)!=b)
                    return(strwrap("The number of columns in pd must be equal to 
                     the length of rn.", width = 60))

                   ewt <- apply(pd, 1, function(x) EWTD(b, R, m, x))                            
                   
                   }
              else
                   {
                    if(length(pd)!=b)
                     return(strwrap("The length of pd must be equal to the 
                     length of rn.", width = 60))
                                
                    ewt <- EWTD(b, R, m, pd)          
                   }
              return(ewt)
             }

################################################################
# Calculate the expected sample size (average sampling number) # 
# using the HD function for the arguments matching the ones in #
# the classes.                                                 #
# pd is a vector or matix.                                     #
# rn is a nonnegtive vector                                    #
################################################################
EWTH <- function(pd, rn, m, N)
             {
               b=length(rn)
               R=rn
               if(is.matrix(pd))
                  {
                   if(ncol(pd)!=b)
                     return(strwrap("The number of columns in pd must be equal to
                      the length of rn.", width = 60))
                   M <- round(pd*N)
                   ewt <- apply(M, 1, function(x) EWTHD(b, R, m, x, N))             
                   
                   }
              else
                   {
                    if(length(pd)!=b)
                     return(strwrap("The length of pd must be equal to the 
                     length of rn.", width = 60))
                    M <- round(pd*N)            
                    ewt <- EWTHD(b, R, m, M, N)          
                   }
              return(ewt)
             }

################################################################
# Calculate the expected sample size (average sampling number) #
# using D function for the arguments defined as Dirichlet      #
# D function.                                                  #
# b, m are positive integers.                                  #
# R, P are vector arguments.                                   #
#                                                              #
#                                                              #
################################################################

EWTD <- function(b, R, m, P)
        {
         if(length(P)!= b | length(R)!= b)
           stop("The length of 'P1'and 'R' must be equal to 'b'.")
         if(any(P>1)|any(P<0))stop("'p' is out of range (0,1)")
         if(sum(P)>1)stop("Sum of 'P' must not be greater than 1.")         
         if(any( R<=0 ))return(0)
         p0 <- 1-sum(P)
         if(all(P==0))
           return(m)
         else
           {
            RR <- c(m, R)
            PP <- c(p0, P)
            R <- RR[PP!=0]
            P <- PP[PP!=0]
            b <- length(R)      
            EPWT <- 0
            for(i in 1:b)
               {
                ri <- R[i]
                Ri <- R[-i]
                pi <- P[i]
                Pi <- P[-i]
                EPWT <- EPWT+(ri/pi)*RDV(b-1, Ri, ri+1, Pi, pi)
               }
             }
         return(EPWT)         
        }




################################################################
# Calculate the expected sample size (average sampling number) #
# using HD function for the arguments defined as Dirichlet     #
# HD function.                                                 #
# b, m, N are positive integers.                               #
# R, M are vector arguments.                                   #
#                                                              #
#                                                              #
################################################################

EWTHD <- function(b, R, m, M, N)
       {
        if( m <= 0 ) stop("'m' must be greater than 0.")
       if(length(M)!= b | length(R)!= b) 
         stop("The length of 'M'and 'R' must be equal to 'b'.")
        EPWT <- 0
        Mg <- N-sum(M)
        Rm <- c(m, R)
        MM <- c(Mg, M)
        for(i in 1:(b+1))
            {
             ri <- Rm[i]
             Ri <- Rm[-i]
             mi <- MM[i]
             Mi <- MM[-i]
             EPWT <-EPWT+(ri*(N+1)/(mi+1))*RHDV(b, Ri, ri+1, Mi, N+1)
            }
         return(EPWT)
        }


