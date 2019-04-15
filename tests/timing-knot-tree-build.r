rm(list=ls())
setwd("/home/marcin")
library(GpGp); library(Matrix); library(RcppParallel)
library(parallel); library(sparseinv); library(fields)
for (nm in list.files('GPvecchia/R',pattern = "\\.[RrSsQq]$")) {
  source(file.path('GPvecchia/R',nm))#; cat("\n")
}
Rcpp::sourceCpp('GPvecchia/src/U_NZentries.cpp')
Rcpp::sourceCpp('GPvecchia/src/MaxMin.cpp')


set.seed(1988)
spatial.dim=2# number of spatial dimensions
n=50**2  # number of observed locs
m=40

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

# n=4^2
# grid.oneside=seq(0,1,length=round(sqrt(n)))
# locs=as.matrix(expand.grid(grid.oneside,grid.oneside))


sig2=1; range=1; smooth=0.5
me.var = 1e-8
covparms =c(sig2,range,smooth)
covfun = function(locs1, locs2=NULL) {
  if(is.null(locs2)){
    sig2*MaternFun(fields::rdist(locs1),covparms)
  } else {
    #c(sig2*exp(-matrix(fields::rdist.vec(locs1, locs2), ncol=1)))
    c(sig2*MaternFun(matrix(fields::rdist.vec(locs1, locs2), ncol=1),covparms))
  }
}
nuggets=rep(me.var,n)

# simulate observations
if(n <= 1e4) {
  Sigma = covfun(locs) + diag(nuggets)
  Sigma.c = chol(Sigma)
  z=as.numeric(t(Sigma.c)%*%rnorm(n))
} else z=rnorm(n)



# mra.options = list(plots=FALSE, r=c(0,16))
# mra.options = list(plots=FALSE, M=4)
# mra.options = list(plots=FALSE, r=c(0,11))
# mra.options = list(plots=FALSE, r=c(1,5), J=10, M=1)


V = vecchia_specify(locs, m, conditioning = 'mra')#, mra.options=mra.options)

vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')



# exact likelihood
const = dim(locs)[1]*log(2*pi)
logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
quad.form2 = as.numeric(t(z) %*% solve(Sigma) %*% z)
neg2loglike = const + logdet + quad.form2
loglik = -neg2loglike/2

print("Likelihood")
print(paste("True: ",loglik,sep=""))
print(paste("Vecchia: ", vecchia_loglik, sep=""))
