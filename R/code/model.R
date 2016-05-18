model.draw.z <- function(beta, n) {
  "
  Inputs:
    beta - vector corruption probabilities
    n - number of records
  Outputs:
    z - array of size n x length(beta)
  "
  L <- length(beta)  # number of fields 
  z <- array(NA, dim=c(n, L))
  for (i in 1:L) {
    z[,i] <- rbinom(n, size=1, prob=beta[i])
  }
  return(z)
}

model.draw.Y.c <- function(thetas, n) {
  "
  Inputs:
    thetas - List of categorical probability vector (not necessarily all the same length)
  Outputs:
    Y - n x length(theta) array. (b/c number of latent individuals can be at most n)
  "
  L <- length(thetas)
  Y.c <- array(NA, dim=c(n,L))
  for (i in 1:L) {
    theta <- thetas[[i]]
    Y.c[,i] <- which(rmultinom(n, size=1, prob=theta)!=0, arr.ind=T)[,1]
  }
  return(Y.c)
}

model.draw.Lambda <- function(n) {
  "
  Inputs:
    n - number of records (and potential latent individuals)
  Outputs:
    Lambda - vector of length n, which latent individual each record belongs to
  "
  Lambda <- sample(1:n, n)
  return(Lambda)
}

model.draw.X.c <- function(Lambda, Y.c, z, thetas) {
  "
  Inputs:
    Lambda - Linkage structure (vector of length n)
    Y.c - Latent entities (n x numberfields array)
    z - Binary corruption array (n x numberfields array)
    thetas - Field probabilities (numberfields length list of probability vectors)
  Outputs:
    X.c - Observed records (n x numberfields array)
  "
  n <- length(Lambda)
  X.c.corruption <- model.draw.Y.c(thetas, n)
  X.c.uncorrupted <- Y.c[Lambda,]
  X.c <- X.c.uncorrupted*(1-z) + X.c.corruption*z
  return(X.c)
}

model.draw.X.s <- function(Lambda, Y.s, z, thetas, S, d, c) {
  "
  Inputs:
  Lambda - Linkage structure (vector of length n)
  Y.s - Latent entities (n x numberfields array)
  z - Binary corruption array (n x numberfields array)
  thetas - Field probabilities (numberfields length list of probability vectors, equivalent to G)
  S - list of strings in dictionary
  d - string edit distance function
  Outputs:
  X.s - Observed records (n x numberfields array)
  "
  n <- length(Lambda)
  m <- length(thetas)
  D <- d(S)
  thetas <- do.call(rbind, thetas)
  X.s.corruption <- matrix(NaN, n, m)
  for (i in 1:n) {
    term1 = thetas*exp(-c*t(D[,Y.s[Lambda[i],]]))
    probs = sweep(term1,1, rowSums(term1), FUN='/')
    probs = split(t(probs), rep(1:nrow(probs), each=ncol(probs)))
    X.s.corruption[i,] <- model.draw.Y.c(probs, 1)
  }
  X.s.uncorrupted <- Y.s[Lambda,]
  X.s <- X.s.uncorrupted*(1-z) + X.s.corruption*z
  return(X.s)
}


# Example
n = 10  # number of records
m_s = 2
m_c = 3
m = m_c + m_s  # number of fields
beta = runif(m)
print(beta)
z <- model.draw.z(beta, n)
print(z)
theta1 = c(0.5, 0.4, 0.1)
theta2 = c(1)
theta3 = rep(0.1, 10)
thetas = list(theta1, theta2, theta3)
Y.c <- model.draw.Y.c(thetas, n)
print(Y.c)
S = list('aa', 'bb', 'cc', 'ddhhe', 'ee', 'ab', 'cd', 'ef', 'gh', 'azaz')
c <- 1
d <- function(S){adist(S)}
theta1_string = rep(0.1, 10)
theta2_string = c(0.5, 0.1, 0.1, 0.05, 0.05, 0.01, 0.02, 0.03, 0.04, 0.1)
thetas_string = list(theta1_string, theta2_string)
Y.s <- model.draw.Y.c(thetas_string, n)
print(Y.s)
Lambda <- model.draw.Lambda(n)
print(Lambda)
X.c <- model.draw.X.c(Lambda, Y.c, z[,(m_s+1):ncol(z)], thetas)
X.s <- model.draw.X.s(Lambda, Y.s, z[,1:m_s], thetas_string, S, d, c)
print(X.c)
print(X.s)



