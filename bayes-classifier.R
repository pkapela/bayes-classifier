
# Script: bayes-classifier.R
# Date: October, 3rd, 2019
# Author: Piotr Kapela

rm(list = ls())
cat("\014")

# Libraries
library(ggplot2) 
library(gridExtra)
library(MASS)

# Reproduce Randomness
set.seed(1)

# Number of Observations
n <- 1000
p <- 2
X <- matrix(ncol=p, nrow=n)

# Population Parameters
############################
# Variance
o0 <- o1 <- o2 <- 0.5

# Mean
mu0 <- c(0,1)
mu1 <- c(1,0)
mu2 <- c(-1,1)

# Correlation Coefficient
rho_fcn <- function(i) {
  return(0.25 + (i - 1)*0.25)
}
rho0 <- rho_fcn(1)
rho1 <- rho_fcn(2)
rho2 <- rho_fcn(3)

# Covariance Matrix 
C0 <- o0^2 * matrix(c(1, rho0, rho0, 1), nrow=2)
C1 <- o1^2 * matrix(c(1, rho1, rho1, 1), nrow=2)
C2 <- o2^2 * matrix(c(1, rho2, rho2, 1), nrow=2)

# Pi
pi0 <- 0.6
pi1 <- 0.3
pi2 <- 0.1

# Custom Functions
get_class_fcn <- function(x) {
  
  max_val <- det_fcn(x, C0, mu0, pi0)
  
  if( max_val < det_fcn(x,C1,mu1,pi1) ) {
    if( max_val < det_fcn(x, C2, mu2, pi2) ){
      return(2)
    } else {
      return(1)
    }
  } else if ( max_val < det_fcn(x, C2, mu2, pi2) ) {
    return(2)
  } else {
    return(0)
  }
}

det_fcn <- function(x, C, mu, pi) {
  return( -0.5 * x %*% solve(C) %*% x + x %*% solve(C) %*% mu -0.5 * mu %*% solve(C) %*% mu - 0.5 * log( det(C) ) + log(pi) )
}

# Data Generation
y <- as.factor(sample(0:2, size=n, replace=TRUE, prob=c(pi0, pi1, pi2)))
                           
for (i in 1:n){
  mu    <- (y[i]==2)*mu2 + (y[i]==1)*mu1 + (y[i]==0)*mu0
  C     <- (y[i]==2)*C2 + (y[i]==1)*C1 + (y[i]==0)*C0
  X[i,] <- mvrnorm(1, mu, C) 
}

# Plotting
plot(X, col=y, pch=16, xlab="X", ylab="Y", main="Scatter Plot")
legend( "topleft", legend=c("Label 0", "Label 1", "Label 2"), col=c("black", "red", "green"), pch=c(16,16))
grid()

# Estimates of Parameters
############################
n0 <- sum(y==0)
n1 <- sum(y==1)
n2 <- sum(y==2)

X0 <- X[y==0,]
X1 <- X[y==1,]
X2 <- X[y==2,]

mu0.hat <- colMeans(X0)
mu1.hat <- colMeans(X1)
mu2.hat <- colMeans(X2)

pi0.hat <- n0/n
pi1.hat <- n1/n
pi2.hat <- n2/n

X0 <- X0 - mu0.hat
X1 <- X1 - mu1.hat
X2 <- X2 - mu2.hat

C0.hat <- cov(X0)
C1.hat <- cov(X1)
C2.hat <- cov(X2)

get_class_est_fcn <- function(x) {
  
  max_val <- det_fcn(x, C0.hat, mu0.hat, pi0.hat)
  
  if(max_val < det_fcn(x,C1.hat, mu1.hat, pi1.hat)) {
    if(max_val < det_fcn(x, C2.hat, mu2.hat, pi2.hat)){
      return(2)
    } else {
      return(1)
    }
  } else if (max_val < det_fcn(x, C2.hat, mu2.hat, pi2.hat)) {
    return(2)
  } else {
    return(0)
  }
}

# Preparing Data for Plots
X.grid <- expand.grid(x=seq(from=min(X[,1]), to=max(X[,1]), length.out=n), y=seq(from=min(X[,2]), to=max(X[,2]), length.out=n))
# Applying Bayes to the Grid
y.hat.grid <- apply(as.matrix(X.grid), 1, get_class_fcn)
y.hat.grid <- apply(as.matrix(X.grid), 1, get_class_est_fcn)

# Creating Data Frames
datafA <- data.frame(X.grid, as.matrix(y.hat.grid))
datafB <- data.frame(y, X)

# Final Plots
p1 <- ggplot(datafB, aes(x=X1, y=X2, colour=y))+geom_point()+ggtitle("Train data")+xlab("X")+ylab("Y")
p2 <- ggplot(datafA, aes(x=x, y=y, colour=y.hat.grid))+geom_point()+ggtitle("Predicted label")+xlab("X")+ylab("Y")

grid.arrange(p1,  p2, nrow =1)
