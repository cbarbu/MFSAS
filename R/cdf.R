#####################################################################
# calculate the cumulative distribution function for multinomial    #
# using recursive algorithms for the Dirichlet J function.          # 
# x is a non-negative integer vector.                               #
# prob is  a vector or matrix.                                      #
# size is a positive integer.                                       #
#                                                                   #                                     
#####################################################################

pmultinom <- function(x, size , prob)
             {
               if(any(x < 0)) 
                  stop("'x' must be non-negative")
               if(size <=0)
                  stop("'size' must be greater than 1.")
               if(any(prob>1)|any(prob<0))
                  stop("'prob' is out of range [0,1]")
               b=length(x)
               if(b==1)
                {
                 if(length(prob)==1)
                   cdf <- JV(x+1, size, prob)
                 else
                  cdf <-sapply(prob, function(P) JV(x+1, size, P))   
                 }
               else
                 {
                   if(is.matrix(prob))  
                      {
                       if(ncol(prob)!=b)
                       return(strwrap("The length of 'x' and the number columns 
                                of 'prob' must be equal."),  width = 60)

                       cdf <-apply(prob, 1, function(P) JV(x+1, size, P))             
                      }
                   else
                      {
                       if(length(prob)!=b)
                       return("'x' and 'prob'  must be equal length vectors.")
             
                       cdf <- JV(x+1, size, prob)
                       }
                 }
            return(cdf)
           }
            
#########################################
#Dirichlet J vector function:           #
# JV function is used to check if the   # 
# input is valid.                       #
# RJV is recursive function for         #
# Dirichlet J function.                 #
#########################################
JV <- function(R, n, P)
      {
       
       if(sum(P) > 1)
       stop("Sum of 'P' must not be greater than 1.")               
       if(any( R<=0 ))return(0)
       if(all(P==0))return(1)
       R <- R[P!=0]
       P <- P[P!=0]
       b <- length(R)  
       jv <- RJV(b, R, n, P)      
       return(jv)
      }

RJV <- function(b, R, n, P)
     {
      if((max(P)== min(P)) & (max(R)== min(R))) ##check if all the values in vector P 
                                                 ##are the same and if all the values
                                                 ##in vector R are the same
          {
           pj = J(b, 0, R[1], n, P[1])
           return(pj)
         }

       R[R>n]=n+1        
       sump=0
       if(R[1]<=0)return(sump)
       else {r1 <- R[1]-1}
       R1 <- R[-1]
       p1 <- P[1]
       P1 <- P[-1]/(1-p1)
       P1[P1>1]=1
       for(i in (0:r1))
          {
           sump <- sump + exp(lgamma(n+1)-lgamma(i+1)-lgamma(n-i+1)+i*log(p1)+(n-i)*log(1-p1)
			+log(RJV(b-1, R1, n-i, P1)))
          }
       return(sump)
     }


#########################################
#J function for the same r and same p.  #
# It is part of Dirichlet J recursive   # 
# function.                             #
#########################################

J <- function(b,j,r,n,p)
      {
       if((j*r) > n)return(0)
                    
        if((j*r)==n)
                   {
	          pj <- exp(lgamma(n+1)-j*lgamma(r+1)+n*log(p))
                     return(pj)
                    }
      
        if(j==b){
                  pj <- 
			exp(lgamma(n+1)-b*lgamma(r+1)-lgamma(n-b*r+1)
			+b*r*log(p)+(n-b*r)*log((1-b*p)))
                  return(pj)
                }

        pj <- (1/(n-j*r))*(n*(1-j*p)*J(b,j,r,n-1,p)-r*(b-j)*J(b,j+1,r,n,p))
        return(pj)
      }

#####################################################################
# calculate the cumulative distribution function for                #
# multivariate hypergeometric using recursive algorithms for the    #
# Dirichlet HJ function.                                            #
# n, N are positive integers.                                       #
# x is a vector of non-negative integers.                           #
# M is a vector or matrix of non-negative integers.                 #
#                                                                   #
#####################################################################

pmultihyper <- function(x, n, M, N)   
             {
              if(any(x < 0)) 
                  stop("'x' must be non-negative")
              if(any(M < 0)) 
                  stop("'M' must be non-negative")
              if(n <=0)
                 stop("'n' must be greater than 0.")
              if(N<=0)
                 stop("'N' must be greater than 0.")
              if(n > N) 
                 stop("'N' must be no less than than 'n'.")
              b=length(x)
              if(b==1)
                {
                 if(length(M)==1)
                   cdf <- HJV(b, x+1, n, M, N)
                 else
                  cdf <-sapply(M, function(y) HJV(b, x+1, n, y, N))    
                 }
               else
                 {
                  if(is.matrix(M))
                     {
                      if(ncol(M)!=b)
                         stop(strwrap("The length of 'x' and the 
                         number columns of 'M' must be equal."),width = 60)
                  
                      cdf <-apply(M, 1, function(y) HJV(b, x+1, n, y, N))                        
                     }
                  else
                     {
                       if(length(M)!=b)
                         stop("'x' and 'M' must be equal length vectors.")
                       cdf <- HJV(b, x+1, n, M, N)           
                     }
                  }
              return(cdf)
             }

#########################################
#Dirichlet HJ vector function:          #
# HJV function is used to check if the  # 
# input is valid.                       #
# RHJV is recursive function for        #
# Dirichlet HJ function.                #
#########################################

HJV <- function(b, R, n, M, N)
      {
       if(sum(M)>N) 
       stop("'N' must be no less than the sum of 'M'.")       
       if(n < 0) return(1) 
       if(any( R<=0 ))return(0)
       hjv <- RHJV(b, R, n, M, N)
       return(hjv)       
      }
RHJV <- function(b, R, n, M, N)
       {
        if(n > N)return(0)    
        if((max(M)== min(M)) & (max(R)== min(R))) ##check if all the values in vector M 
                                                 ##are the same and if all the values 
                                                 ##in vector R are the same
         {
           phj = HJ(b, 0, R[1], n, M[1], N)
           return(phj)
         }
       R[R>n]=n+1        
       sump=0
       if(R[1]<1)return(sump)
       else {r1 <- R[1]-1}
       R1 <- R[-1]
       m1 <- M[1]
       M1 <- M[-1]
       for(i in (0:min(m1,r1)))
          {
           sump <- sump + (((choose(m1,i))*
           (choose(N-m1,n-i)))/(choose(N, n)))*RHJV(b-1, R1, n-i, M1, N-m1)
          }
       return(sump)
       }

##########################################
#HJ function for the same r and same M.  #
# It is part of Dirichlet HJ recursive   # 
# function.                              #
########################################## 

HJ <- function(b, j, r, n, M, N)
      {           
       if(n > N)return(0)
       if(r <= 0)return(0)
       if(r >= M)return(1)
       if((n-b*r) >= (N-b*M))return(0)  #add this              
       if((j*r) > n)return(0)

        if((j*r)==n)
                   {
                     ph <- ((choose(M, r))^j)/(choose(N, n))
                     return(ph)
                    }
      
        if(j==b){
                  ph <- ((choose(M, r))^b)*(choose((N-b*M),
                        (n-b*r)))/(choose(N, n))
                  return(ph)
                }


        ph <- (1/(n-j*r))*(n*(1-((j*(M-r))/(N-n+1)))*HJ(b,j,r,
               n-1,M,N)-r*(b-j)*HJ(b,j+1,r,n,M,N))

        return(ph)
      }

####################################################################
# calculate the cumulative distribution function for negative      #
# multinomial using recursive algorithms for the Dirichlet         #
# D function.                                                      # 
# m is a positive integer.                                         #
# x is a vector of positive integers                               #
# prob is a vector or matrix.                                      #
#                                                                  #        
####################################################################

pnmultinom <- function(x, m, prob)
             {
               if(any(x < 0)) 
                  stop("'x' must be non-negative")
               if(m <=0)
                  stop("'m' must be greater than 0.")
               if(any(prob>1)|any(prob<0))
                  stop("'prob' is out of range [0,1]")
               R=x+1
               b=length(x)
               if(b==1)
                {
                 if(length(prob)==1)
                   {
                   if(prob >= 1) stop("'prob' must be less than 1.")
                   cdf <- DV(b, R, m, prob)
                   }
                 else
                   {
                  if(any(prob >=1 )) stop("'prob' must be less than 1.")
                  cdf <-sapply(prob,function(P) DV(b, R, m, P))   
                   }
                 }
               else
                 {
                  if(is.matrix(prob))
                    {
                     if(ncol(prob)!=b)
                       stop(strwrap("The length of 'x' and the number columns 
                              of 'prob' must be equal."), width = 60)
                    if(any(rowSums(prob) >= 1)) stop("Sum of 'prob' must be less than 1.")  
                    cdf <-apply(prob, 1, function(P) DV(b, R, m, P))                        
                    }
                 else
                    {
                     if(length(prob)!=b)
                        stop("'x' and 'prob'  must be equal length vectors.")             
                     if(sum(prob) >= 1) stop("Sum of 'prob' must be less than 1.")
                     cdf <- DV(b, R, m, prob)
                    }
                 }
            return(cdf)
           }

#########################################
#Dirichlet D vector function:           #
# DV function is used to check if the   # 
# input is valid.                       #
# RDV is recursive function for         #
# Dirichlet D function.                 #
#########################################

DV <- function(b, R, m, P)
       {
       if(sum(P)>1)
       stop("Sum of 'P' must not be greater than 1.")
       if(any( R<=0 ))return(0)
       p0 <- 1-sum(P)
       P1 <- P
       dv <- RDV(b, R, m, P1, p0)
       return(dv)
      }
RDV <- function(b, R, m, P1, p0)
      { 
       if(p0==0) return(0)
       if(p0==1) return(1)     
       if((max(P1)== min(P1)) & (max(R)== min(R))) ##check if all the values in vector P 
                                                 ##are the same and if all the values
                                                 ##in vector R are the same
          
          {
           a=P1[1]/p0
           rdv = D(b, 0, R[1], m, a)
           return(rdv)
          }
       sump=0      
       r1 <- R[1]-1
       R1 <- R[-1]
       p1=P1[1]
       P1 <- P1[-1]
       for(i in 0:r1)
          {
           sump <- 
		sump + exp(lgamma(m+i)-lgamma(i+1)-lgamma(m)+i*log(p1/(p0+p1))
			+m*log(p0/(p0+p1))+log(RDV(b-1, R1,m+i, P1, p0+p1)))
          }
       return(sump)
      }

##########################################
#D function for the same r and same p.   #
# It is part of Dirichlet D recursive    # 
# function.                              #
##########################################

D <- function(b, j, r, m, a)
     {
       if( m <= 0 )return(0)
       if(j==b)
          {
           sump <- 
		exp((lgamma(m+r*b) - b*lgamma(r+1) - lgamma(m)) 
			+m*log(1/(1+a*b))+b*r*log(a/(1+a*b)))
           return(sump)
          }
       if(m > r)
         {
          temp=0
          for(i in 1:r)
              {
               temp <- temp + exp(lgamma(m-i)-lgamma(r-i+1)-lgamma(m-r)-i*log(a)+log(D(b, j+1, r, m-i, a)))
              }
          sump <- (1/choose((m-1), r))*temp
          return(sump)
         }
       sump <- (1/(m+j*r))*(m*(1+j*a)*D(b, j, r, m+1, a)+r*(b-j)*D(b, j+1, r, m, a))
       return(sump)
     }

#############################################################################
# calculate the cumulative distribution function for negative multivariate  #
# hypergeometric using recursive algorithms for the Dirichlet HD function.  #
# b, m, N are positive integers.                                            #
# M, R are vector arguments.                                                #
#                                                                           #
#                                                                           #
#############################################################################

pnmultihyper <- function(x, m, M, N)  
             {
              if(any(x < 0)) 
                  stop("'x' must be non-negative")
              if(any(M < 0)) 
                  stop("'M' must be non-negative")
              if( m <= 0 ) stop("'m' must be greater than 0.")

               R=x+1
               b=length(x)
               if(b==1)
                {
                 if(length(M)==1) 
                    {
                     if(m > N-M) stop("'m' must be less than or equal to N-M.")
                     cdf <- HDV(b, R, m, M, N)
                    } 
                  else 
                    {
                     if(any(m > N-M)) stop("'m' must be less than or equal to N-M.")
                     cdf <-sapply(M,function(y) HDV(b, R, m, y, N))    
                    }
                 }
               else
                 {               
                  if(is.matrix(M))
                     {
                       if(ncol(M)!=b)
                          stop(strwrap("The length of 'x' and the number columns 
                                of 'M' must be equal."),  width = 60)
                       if(any(m > N-rowSums(M))) stop("'m' must be less than or equal to N-sum(M).")  
                       cdf <-apply(M, 1, function(y) HDV(b, R, m, y, N))
                      }
                  else
                      {
                       if(length(M)!=b)
                         stop("'x' and 'M'  must be equal length vectors.")
                       if(m > N-sum(M)) stop("'m' must be less than or equal to N-sum(M).")
                       cdf <- HDV(b, R, m, M, N)
                      }           
                   }
              return(cdf)
             }


##########################################
#Dirichlet HD vector function:           #
# HDV function is used to check if the   # 
# input is valid.                        #
# RHDV is recursive function for         #
# Dirichlet HD function.                 #
##########################################

HDV <- function(b, R, m, M, N)
       {
       if(any( R<=0 ))return(0)
       hdv <- RHDV(b, R, m, M, N)
       return(hdv)
       }

RHDV <- function(b, R, m, M, N)
      {
       if((max(M)== min(M)) & (max(R)== min(R))) ##check if all the values in vector M 
                                                 ##are the same and if all the values 
                                                 ##in vector R are the same
         {
           rhdv = HD(b, 0, R[1], m, M[1], N)
           return(rhdv)
         }
       r1 <- R[1]-1
       R1 <- R[-1]
       m1 <- M[1]
       Mm <- N-sum(M)
       if(Mm < m)return(0)
       M1 <- M[-1]
       sump=0
       for(i in (0:r1))
          {
           sump <- sump + (choose(m1,i)/((choose((N-sum(M1)), (m+i)))*(m+i)))*RHDV(b-1, 
                           R1, (m+i), M1, N)
          }
       rhdv <- (choose(Mm, m))*m*sump
       return(rhdv)
      }


##########################################
#HD function for the same r and same M.  #
# It is part of Dirichlet HD recursive   # 
# function.                              #
##########################################
HD <- function(b, j, r, m, M, N)
       {
       if(r > M) return(1)
       if(N <= m+b*r)return(0)    
       Mm <- N-b*M
       if(Mm < m)return(0)       
       if(j==b)
          {
           phd <- (m/(m+b*r))*(((choose(M,r))^b)*(choose(Mm, m))/choose(N, (m+b*r)))
           return(phd)
          }
       if(m > r)
         {
          temp=0
          for(i in 0:(r-1))
              {
              
               temp <- temp + 
               (m/(m+i-r))*(choose(M,i)/choose(Mm, (m+i-r)))*HD(b, j+1, r, m+i-r, M, N)
             
              }
          phd <- (choose(Mm, m)/choose(M,r))*temp
          return(phd)
         } 
       if(Mm == m)
         {
          temp=0
          for(i in m:(m+b*(r-1)))
              {
          
               temp <- temp + 
               choose(b*M,i-m)/i/choose(N,i)*HJ(b,0,r,i-m,M,b*M)
             
              }
          phd <- m*choose(N-b*M,m)*temp
          return(phd)
         }
         phd <- (1/(m+j*r))*(m*(1+((j*(M-r))/(Mm-m)))
                 *HD(b, j, r, m+1, M, N)+ r*(b-j)*HD(b, j+1, r, m, M, N))
         return(phd)
       }

