# test of ERE_step_cycle
library("systemicrisk")
set.seed(123787)
gen1stepGibbs <- function(p,lambda,steps=5,k=3){
  n <- dim(p)[1]
  M <- matrix(rbinom(n^2,size=1,prob=p)*rexp(n^2,rate=lambda),nrow=n)
  Mstep <- cloneMatrix(M)
  res <- list(M=M)
  for (i in 1:steps){
    r=sample(0:(n-1),k,replace=FALSE)
    col=sample(0:(n-1),k,replace=FALSE)
    tryCatch(
             ERE_step_cycle(r=r,c=col,L=Mstep, p=p,lambda=lambda),error=function(cond){browser()})
    if (!all(Mstep>=0))
      expect_true(all(Mstep>=0))
  }
  res$Mstep=Mstep
  res
}

test <- function(x,p,lambda){
  n <- dim(p)[1]

  expect_equal(rowMeans(x[1:n^2,]==0),rowMeans(x[(n^2+1):(2*n^2),]==0),tolerance=1e-2,check.attributes = FALSE)
  expect_equal(rowMeans(x[1:n^2,]),rowMeans(x[(n^2+1):(2*n^2),]),tolerance=1e-2,check.attributes = FALSE)

#  expect_equal(table(colSums(x[1:n^2,]==0)),table(colSums(x[(n^2+1):(2*n^2),]==0)),tolerance=1e-2,check.attributes=FALSE)
  for (i in 1:n^2){
    if (p[i]>0){
      expect_equal(mean(x[i,x[i,]>0]),1/lambda[i],tolerance=2e-2,label=paste("Mean conditional on positive value, Initial sample, i=",i,"\n"))
      expect_equal(mean(x[i+n^2,x[i+n^2,]>0]),1/lambda[i],tolerance=2e-2,label=paste("Mean conditional on positive value, after Gibbs steps, i=",i,"\n"))
    }
  }
}
testrun <- function(p,lambda,nrep=1e5,k=2){
  n <- dim(p)[1]
  x <- replicate(nrep,unlist(gen1stepGibbs(p,lambda,k=k)))
  test(x,p,lambda)
}


## test_that("Individual small scale tests",{
##   res <- matrix(0.,nrow=3,ncol=3)
##   p <- matrix(0.2,nrow=3,ncol=3)
##   diag(p)=0;
##   lambda <- matrix(1,nrow=3,ncol=3)
##   gen_3by3_diagonal(r=c(0,0,0),c=c(0,0,0),lambda,p,res,eps=1e-10)
##   expect_equal(res,matrix(0,nrow=3,ncol=3))

##   for (i in 1:3)
##     for (j in 1:3){
##       r <- rep(0,3); col <- rep(0,3);
##       r[i] <- 1;
##       col[j] <- 1;
##       if (i!=j){
##         gen_3by3_diagonal(r=r,c=col,lambda,p,res,eps=1e-10)
##         expect_equal(res,outer(r,col))
##       }else{
##         expect_error(gen_3by3_diagonal(r=r,c=col,lambda,p,res,eps=1e-10)) ## no matrix exists
##       }
##     }
## })

test_that("All probabilities equal, ps equal to 1, lambdas equal to 1",{
  skip_on_cran()

  p <- matrix(1,nrow=3,ncol=3)
  diag(p)=0;
  lambda <- matrix(1,nrow=3,ncol=3)
  testrun(p=p,lambda=lambda,k=3)
  testrun(p=p,lambda=lambda,k=2)
})



test_that("All probabilities equal, lambdas equal to 1",{
  skip_on_cran()

  p <- matrix(0.4,nrow=3,ncol=3)
  lambda <- matrix(1,nrow=3,ncol=3)
  testrun(p=p,lambda=lambda,k=2)
  testrun(p=p,lambda=lambda,k=3)
})

test_that("All probabilities equal, lambdas equal to 1",{
  skip_on_cran()

  p <- matrix(0.9,nrow=3,ncol=3)
  lambda <- matrix(1,nrow=3,ncol=3)
  testrun(p=p,lambda=lambda,k=2)
  testrun(p=p,lambda=lambda,k=3)
})


test_that("All probabilities and lambdas equal",{
  skip_on_cran()

  p <- matrix(0.8,nrow=3,ncol=3)
  lambda <- matrix(2,nrow=3,ncol=3)
  testrun(p=p,lambda=lambda,k=2)
  testrun(p=p,lambda=lambda,k=3)
})

test_that("Example with varying p and lambda",{
  skip_on_cran()

  p <- matrix(seq(0.1,1,length.out=9),nrow=3,ncol=3)
  lambda <- matrix(1/(1:9),nrow=3,ncol=3)
  testrun(p=p,lambda=lambda,k=2)
  testrun(p=p,lambda=lambda,k=3)
})

test_that("Example with varying high p and lambda",{
  skip_on_cran()

  p <- matrix(seq(0.6,1,length.out=9),nrow=3,ncol=3)
  diag(p)=0;
  lambda <- matrix(1/(1:9),nrow=3,ncol=3)
  testrun(p=p,lambda=lambda)
})

