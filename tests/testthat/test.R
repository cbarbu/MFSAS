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

test_that("ppamultinom is ok with dmultinom",{
X <- c(0,1,1,0,1)
size <- 5
prob <- c(0.1,0.2,0.15,0.1,0.45)
expect_equal(ppamultinom(X,3,prob),dmultinom(X,3,prob))

X <- c(0,0,0,0,1)
expect_equal(ppamultinom(X,size,prob),dmultinom(X*size,size,prob))

X <- c(0,1,0,1,1)
size <- 4
p <- dmultinom(c(0,2,0,1,1),size,prob) + dmultinom(c(0,1,0,2,1),size,prob) + dmultinom(c(0,1,0,1,2),size,prob)
expect_equal(ppamultinom(X,size,prob),p)

size <- 5
p <- dmultinom(c(0,3,0,1,1),size,prob) + dmultinom(c(0,2,0,2,1),size,prob) + dmultinom(c(0,2,0,1,2),size,prob) +
    dmultinom(c(0,1,0,3,1),size,prob) + dmultinom(c(0,1,0,1,3),size,prob) + dmultinom(c(0,1,0,2,2),size,prob)
expect_equal(ppamultinom(X,size,prob),p)
})


