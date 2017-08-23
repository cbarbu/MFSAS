## the k-1 length of p and x is a commodity allowing to compute things
## leaving one parameter free
## it can be specified with k parameters and works good

test_that("pmultinom on k-1 equal pmultinom on k if kth a non informative (>n))",{
library("testthat")
# specifying all k
p <- c(0.1,0.9)
n <- 5
x <- c(2,3)
pAll1 <- pmultinom(x,n,p)

x <- c(2,4)
pAll2 <- pmultinom(x,n,p)

x <- c(2,n)
pAll3 <- pmultinom(x,n,p)


# specifying only k-1 
pFirst <- pmultinom(x[1],n,p[1])

expect_equal(pAll3,pFirst)
})

test_that("pmultinom coherent with dbinom",{
## comparing against the dbinomial
p <- c(0.1,0.7,0.2)
n <- 15
xLow <- c(2,10,3)
pLow <- pmultinom(xLow,n,p)
dLow <- dmultinom(xLow,n,p)
expect_lt(abs(pLow-dLow),10e-15)

xHigh <- c(3,10,3)
pHigh <- pmultinom(xHigh,n,p)

xPos1 <- c(3,10,2)
dPos1 <- dmultinom(xPos1,n,p)

xPos2 <- c(3,9,3)
dPos2 <- dmultinom(xPos2,n,p)

expect_lt(abs(pHigh-pLow-dPos1-dPos2),10e-15)
})

