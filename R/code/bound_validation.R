"
Emperical Validation of Theoretical Lower Bounds
Author: mbarnes1@cs.cmu.edu
"
rm(list=ls())
setwd('~/Documents/trafficjam/beka/ebLink/R/code/')
source('model.R')
source('conditionals.R')
source('ebGibbsSampler.R')

# Experiment parameters
n = 20
betas <- seq(from = 0.1, to = 0.2, by = 0.1)
thetas <- list(rep(0.1, 10), rep(0.1, 10), rep(0.1, 10))

# Gibbs sampler parameters
a <-1
b <- 999
c <- 1
d <- function(string1,string2){adist(string1,string2)}
  
L <- length(thetas)  # number of fields
# Sweep through experiment parameters
error_exactsampling <- rep(NA, length(betas))
bound_expected_error <- rep(NA, length(betas))
error_gibbs <- rep(NA, length(betas))

for (ibeta in 1:length(betas)) {
  # Draw synthetic data
  beta <- rep(betas[ibeta], L)
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
  bound_expected_error[ibeta] <- pr_error
  
  # Run experiment (exact sampling)
  Lambda.exact <- conditionals.draw.lambda(Y.s,Y.c,z)
  error_exactsampling[ibeta] <- sum(Lambda != Lambda.exact)/n
  
  # Run experiment (full Gibbs sampler)
  Lambda.gibbs <- rl.gibbs(file.num=rep(1,n),X.s=X.s,X.c=X.c,num.gs=100,a=a,b=b,c=c,d=d, M=n, num.gs.burnin=20, Y.s.0=Y.s, Y.c.0=Y.c)
  error_gibbs[ibeta] <- sum(Lambda != Lambda.gibbs[nrow(Lambda.gibbs),])/n
}

# Print the results
print('Bound expected error')
print(bound_expected_error)
print('Exact sampling error')
print(error_exactsampling)
print('Gibbs sampling')
print(error_gibbs)


# Plot the results
#plot(betas, error_exactsampling, type="b", col="red", ylim=c(0, max(max(error_gibbs), max(error_exactsampling))), ann=FALSE)
#lines(betas, bound_expected_error, type='b', col="blue")
#lines(betas, error_gibbs, type='b', col="green")
#title(main='Experimental Bound Validation', xlab='beta', ylab='Error')

#legend('topleft', c("Exact Sampling","Bound", "Gibbs"), cex=1, 
#       col=c("red","blue", "green"), lty=1:3, lwd=2, bty="n")




