"
Emperical Validation of Theoretical Lower Bounds
Author: mbarnes1@cs.cmu.edu
"
rm(list=ls())
setwd('/home/scratch/mbarnes1/ebLink/R/code/')

compute_errors_wrapper <- function(x) {
  result <- tryCatch({
    errors <- compute_errors(x)
    return(errors)
    }, error = function(e) {
    return(e)
    })
  return(result)
}

compute_errors <- function(param) {
  print(param)
  beta = param[[1]]
  print(beta)
  n = param[[2]]
  print(n)
  thetas = param[[3]]
  source('model.R')
  source('conditionals.R')
  source('ebGibbsSampler.R')
  
  # Draw synthetic data
  L <- length(thetas)  # number of fields
  beta <- rep(beta, L)
  Lambda <- model.draw.Lambda(n)
  Y.c <- model.draw.Y.c(thetas, n)
  z <- model.draw.z(beta, n)
  X.c <- model.draw.X.c(Lambda, Y.c, z, thetas)
  Y.s <- array(rep('x',n), dim=c(n,1))  # dummy strings, sampler will not run without this
  X.s <- array(Y.s[Lambda,], dim=c(n,1))
  z <- cbind(rep(0,n), z)  # zero distortion on dummy string
  
  # Compute theoretical lower bound
  gamma <- 0
  minthetas <- sapply(thetas,min)
  for (i in 1:n) {  # gamma is upper bound on largest KL divergence between any two linkages
    for (prime in 1:n)  {
      ind <- Y.c[i,] != Y.c[prime,]
      gamma_ip <- sum(ind*(1-beta)*log(1/(beta*minthetas)))
      gamma <- max(gamma, gamma_ip)
    }
  }
  pr_error <- 1-(gamma + log(2))/log(n+1)
  bound_expected_error <- pr_error
  
  # Run experiment (exact sampling)
  Lambda.exact <- conditionals.draw.lambda(Y.s,Y.c,z, X.s, X.c)
  error_exactsampling <- sum(Lambda != Lambda.exact)/n
  
  # Run experiment (full Gibbs sampler)
  Lambda.gibbs <- rl.gibbs(file.num=rep(1,n),X.s=X.s,X.c=X.c,num.gs=10000,a=a,b=b,c=c,d=d, M=n, num.gs.burnin=10000, Y.s.0=Y.s, Y.c.0=Y.c)
  error_gibbs <- sum(Lambda != Lambda.gibbs[nrow(Lambda.gibbs),])/n
  errors = c(bound_expected_error, error_exactsampling, error_gibbs)
  return(errors)
}

library(parallel)


# Experiment parameters
n_params = 10
n_vec <- seq(from = 10, to = 1000, length = n_params)
betas <- seq(from = 0.6, to = 0.6, length = n_params)
thetas <- rep(list(list(rep(0.1, 10), rep(0.1, 10), rep(0.1, 10))), n_params)

params <- vector('list', n_params)
for (i in 1:n_params) {
  params[i] <- list(list(betas[i], n_vec[i], thetas[[i]]))
}

# Gibbs sampler parameters
a <-1
b <- 999
c <- 1
d <- function(string1,string2){adist(string1,string2)}
  
# Experiment multiprocessing
no_cores <- detectCores()
cl <- makeCluster(min(no_cores, length(betas)))
clusterExport(cl, c("a", "b", "c", "d"))
errors = parLapply(cl, params, compute_errors_wrapper)
print(errors)
stopCluster(cl)

errors = do.call(rbind, errors)
bound_expected_error = errors[,1]
error_exactsampling = errors[,2]
error_gibbs = errors[,3]

# Print the results
print('Parameters were (in order beta, n, theta)')
print(params)
print('Bound expected error')
print(bound_expected_error)
print('Exact sampling error')
print(error_exactsampling)
print('Gibbs sampling')
print(error_gibbs)


# Plot the results
plot(n_vec, error_exactsampling, type="b", col="red", ylim=c(0, max(max(error_gibbs), max(error_exactsampling))), ann=FALSE, pch=21, lwd=4)
lines(n_vec, bound_expected_error, type='b', col="blue", lty=2, pch=22, lwd=4)
lines(n_vec, error_gibbs, type='b', col="green", lty=3, pch=23, lwd=4)
title(main='Experimental Bound Validation', xlab='n', ylab='Error')

legend('bottomright', c("Exact Sampling","Bound", "Gibbs"), cex=1, 
       col=c("red","blue", "green"), lty=1:3, pch=21:23, lwd=2, bty="n")


