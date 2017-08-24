#' @title Probability of presence/absence for multinomial choices
#' @description Given a vector \code{X} of presence (1) or absence (0) or unknown (NA), a number of trials \code{size} and a vector of probabilities \code{prob} returns the probability to have this 
#' @param X The vector of presence (1), absence (0) or unknown (NA) of max size \code{k} the number of variables
#' @param size The number of trials drawing of the k classes
#' @param prob A vector of probabilities corresponding to \code{X}. If not amounting to 1, this function assumes
#'             the existence of an additional class for the differential. 
#' @details Working by recurence considering that p(0,1,1,0,1) = p(0,1,1,0,anything) - p(0)
#'          where NA stands for anything in the calculus. As a consequence, the computation time is 
#'          proportional to 2^n_1 with n_1 the number of choices with presence. It quickly explodes.
#' @examples
#' X <- c(0,1,1,0,1)
#' size <- 5
#' prob <- c(0.1,0.2,0.15,0.1,0.45)
#' ppamultinom(X,size,prob)
#' @export
ppamultinom <- function(X,size,prob,init=TRUE){
    if(init){ # direct call from the user => check things
        if(sum(prob)>1){
            stop("sum(prob)>1, it must be <=1")
        }
        if(length(X)!=length(prob)){
            stop("length(X)!=length(prob)")
        }
        X <- as.integer(X)
        prob <- as.numeric(prob)
    }
    if(sum(X,na.rm=TRUE)>size){
        p <- 0
        return(p)
    }# else if((length(which(X==1|is.na(X)))==0) && size>0){ # delt with ok hereafter


    iNotNull <- which(X!=0)
    if(length(iNotNull)>0){
        Xall <- Xzero <- X
        Xall[iNotNull[1]] <- NA
        Xzero[iNotNull[1]] <- 0
        p <- ppamultinom(Xall,size,prob,init=FALSE)-ppamultinom(Xzero,size,prob,init=FALSE)
    }else{
        p0 <- sum(prob[which(X==0)])
        # pNA <- sum(prob[which(is.na(X))])
        p <- dbinom(size,size,1-p0)
    }

    return(p)
}

