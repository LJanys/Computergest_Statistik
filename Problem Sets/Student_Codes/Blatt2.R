
library(MASS)
install.packages(MASS)

################################################################################
##data generating process#######################################################
n <- 1000 

beta <- c(1, 0.1, 0.1)
covMat <- matrix(c(5, 2.5, 2.5, 4), 2, 2)
covMat
epsilon <- rnorm(n, 0, 10)
epsilon

#you need to install package MASS to use the function mvrnorm
# ?mvrnorm
#mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = TRUE, EISPACK = FALSE)

data <- mvrnorm(n = n, mu = c(0, 0), Sigma = covMat, empirical = TRUE)
data

colMeans(data) #near zero
var(data) #correct values

X <- cbind(rep(1, n), data)
X
Y <- X%*%beta+epsilon
Y

beta_hat <- solve((t(X)%*%X))%*%t(X)%*%Y
beta_hat

Y_hat <- X%*%beta_hat
Y_hat

################################################################################
##Overall-F-Test################################################################
#Nullhyp: beta=0

alpha <- 0.05
k <- length(beta) #relevant for degrees of freedom
k
explained <- sum((Y_hat-mean(Y))**2)
explained
resi <- sum((Y-Y_hat)**2) #residuals whithin the sample
resi

F_value <- (explained/(k-1))/(resi/(n-k))
F_value

#Quantile F-distribution
#?qf
#pf(q, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE)
F_crit <- qf(1-alpha, k-1, n-k)
F_crit

reject_F <- (F_value>F_crit)
reject_F

  #interpretation: 
  #TRUE means there is a significant effect, FALSE means the model is useless overall

################################################################################
##t-Test for regression, which variable is significant?#########################

##Nullhyp beta[1]=0###############

#alpha <- 0.05
#k <- length(beta) #relevant for degrees of freedom
sigma_hat <- sqrt(resi/(n-k))
t_value <- (beta_hat[1])/(sqrt(sigma_hat**2*solve(t(X)%*%X)[1,1])) 
t_value
  #Quantile t-distribution
  #?qt()
  #qt(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
t_crit <- qt(1-0.5*alpha, n-k) #we search t_0.975-quantile, because it is a two-sided test and t is symmetrical
t_crit
# t_crit <- abs(qt(0.5*alpha, n-k)) would be correct too, because the t-distribution is symmetrical
reject_t <- (abs(t_value)>t_crit)
reject_t #interpretation: accept NULLhyp, if reject_t==FALSE

##Nullhyp beta[2]=0###############

t_value <- (beta_hat[2])/(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2])) 
t_value
reject_t <- (abs(t_value)>t_crit)
reject_t #interpretation: accept NULLhyp, if reject_t==FALSE

##Nullhyp beta[3]=0##############

t_value <- (beta_hat[3])/(sqrt(sigma_hat**2*solve(t(X)%*%X)[3,3]))
t_value
reject_t <- (abs(t_value)>t_crit)
reject_t

################################################################################
##Confidence Intervall and coverage probability#################################

##CI for beta
#here we use the formula we had from the lecture to calculate our CI
t_crit
#for b[2]
lower= (beta_hat[2])-(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2]))*t_crit
upper= (beta_hat[2])+(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2]))*t_crit
lower
beta_hat[2]
upper 
#for b[3]
lower= (beta_hat[3])-(sqrt(sigma_hat**2*solve(t(X)%*%X)[3,3]))*t_crit
upper= (beta_hat[3])+(sqrt(sigma_hat**2*solve(t(X)%*%X)[3,3]))*t_crit
lower 
beta_hat[3]
upper 

################################################################################
##Coverage probability##########################################################

#simulation method for estimating the coverage probability of a confidence interval: 
#The simulation method has three steps:
#first: Simulate many samples of size n from the population.
#second: Compute the confidence interval for each sample.
#third: Compute the proportion of samples for which the (known) population parameter is contained in the confidence interval. 
#That proportion is an estimate for the empirical coverage probability for the CI.


#first: Simulate many samples of size n from the population.
many <- 100
confidence_interval <- matrix(rep(NA, 4*many), many, 4) 
counter_beta2 <- 0 #we need the variable later on
counter_beta3 <- 0
  #later on, evry row should contain lower and upper (of CI) for both, beta[2] and beta[3].
for(i in 1:many){
    epsilon <- rnorm(n, 0, 10)
    data <- mvrnorm(n = n, mu = c(0, 0), Sigma = covMat, empirical = TRUE)
    X <- cbind(rep(1, n), data)
    Y <- X%*%beta+epsilon
    beta_hat <- solve((t(X)%*%X))%*%t(X)%*%Y
    Y_hat <- X%*%beta_hat
    resi <- sum((Y-Y_hat)**2)
    sigma_hat <- sqrt(resi/(n-k))
  #second: Compute the confidence interval for each sample.
    #lower beta[2]
    confidence_interval[i,1] <-(beta_hat[2])-(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2]))*t_crit
    #upper beta[2]
    confidence_interval[i,2] <-(beta_hat[2])+(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2]))*t_crit
    #lower beta[3]
    confidence_interval[i,3] <-(beta_hat[3])-(sqrt(sigma_hat**2*solve(t(X)%*%X)[3,3]))*t_crit
    #upper beta[3]
    confidence_interval[i,4] <-(beta_hat[3])+(sqrt(sigma_hat**2*solve(t(X)%*%X)[3,3]))*t_crit
  #third: Compute the proportion of samples for which the (known) population parameter is contained in the confidence interval. 
        #That proportion is an estimate for the empirical coverage probability for the CI.
    if(confidence_interval[i,1]<=beta[2]&beta[2]<=confidence_interval[i,2]){counter_beta2 <- counter_beta2+1 }
    else{}
    if(confidence_interval[i,3]<=beta[3]&beta[3]<=confidence_interval[i,4]){counter_beta3 <- counter_beta3+1}
    else{}
}
estimate_coverage_probability_beta2 <- (counter_beta2/many)
estimate_coverage_probability_beta3 <- (counter_beta3/many)
#confidence_interval
estimate_coverage_probability_beta2
estimate_coverage_probability_beta3


################################################################################
##Simulation Study##############################################################
##function: Are both NULLhyp (t and F distribution) rejected?###################
#variables and values, which never change and that are arguments in the function
n <- 1000
beta <- c(1, 0.1, 0.1)
covMat <- matrix(c(5, 2.5, 2.5, 4), 2, 2)
alpha <- 0.05

both_reject <- function(n, alpha, beta, covMat){
  #data generating procress
  epsilon <- rnorm(n, 0, 10)
  data <- mvrnorm(n = n, mu = c(0, 0), Sigma = covMat, empirical = TRUE)
  X <- cbind(rep(1, n), data)
  Y <- X%*%beta+epsilon
  beta_hat <- solve((t(X)%*%X))%*%t(X)%*%Y
  Y_hat <- X%*%beta_hat
  #F-statistic
  k <- length(beta) #relevant for degrees of freedom
  F_crit <- qf(1-alpha, k-1, n-k)
  explained <- sum((Y_hat-mean(Y))**2)
  resi <- sum((Y-Y_hat)**2) #residuals whithin the sample
  F_value <- (explained/(k-1))/(resi/(n-k))
  reject_F <- (F_value>F_crit)
  #t-statistics, only for beta[2]
  sigma_hat <- sqrt(resi/(n-k))
  t_value <- (beta_hat[2])/(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2]))
  t_crit <- qt(1-0.5*alpha, n-k) 
  reject_t <- (abs(t_value)>t_crit)
  #
  return(reject_F&reject_t)
}
both_reject(1000, alpha, beta, covMat)

##loop, repeat 1000 times########################################################

reject_both_vector <- rep(NA, 1000) #initializing

for(i in 1:1000){
  reject_both_vector[i] <- both_reject(1000, alpha, beta, covMat)
}
reject_both_vector
#Interpretation: 
#In most cases we get "FALSE", although our true beta was beta=c(1, 0.1, 0.1). 
#Beta is not zero and our NULLhypothesis is wrong. 
#If we get "False", at least one of the testing statistics don't reject the NULLhypothesis.
#Unfortunately we get "TRUE" only rarely.
#

################################################################################
###simulation study:violating model assumptions#################################

#heteroskedastic
#rewriting DGP:

both_reject_hete <- function(n, alpha, beta, covMat){
  #data generating procress
  epsilon_hete=rep(NA, n)
  for(i in 1:n){
   epsilon_hete[i]<-rnorm(1,0,10*(X[i,2])**2)
  }
    #epsilon_hete is heteroskedastic
  data <- mvrnorm(n = n, mu = c(0, 0), Sigma = covMat, empirical = TRUE)
  X <- cbind(rep(1, n), data)
  Y <- X%*%beta+epsilon_hete
  beta_hat <- solve((t(X)%*%X))%*%t(X)%*%Y
  Y_hat <- X%*%beta_hat
  #F-statistic
  k <- length(beta) #relevant for degrees of freedom
  F_crit <- qf(1-alpha, k-1, n-k)
  explained <- sum((Y_hat-mean(Y))**2)
  resi <- sum((Y-Y_hat)**2) #residuals whithin the sample
  F_value <- (explained/(k-1))/(resi/(n-k))
  reject_F <- (F_value>F_crit)
  #t-statistics, only for beta[2]
  sigma_hat <- sqrt(resi/(n-k))
  t_value <- (beta_hat[2])/(sqrt(sigma_hat**2*solve(t(X)%*%X)[2,2])) 
  t_crit <- qt(1-0.5*alpha, n-k)
  reject_t <- (abs(t_value)>t_crit)
  #
  return(reject_F&reject_t)
}
#both_reject_hete(1000, alpha, beta, covMat)

##loop, repeat 1000 times########################################################
reject_both_vector_hete <- rep(NA, 1000)

for(i in 1:1000){
  reject_both_vector_hete[i] <- both_reject_hete(1000, alpha, beta, covMat)
}
reject_both_vector_hete








