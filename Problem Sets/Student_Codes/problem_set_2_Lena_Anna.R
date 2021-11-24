rm(list = ls()) #removes all variables from current environment
#install.packages("plotrix") #more advanced plots                          
library("plotrix")
library("MASS")
set.seed(100) #set seed so sample can be replicated

######## Problem_set_2 #########

###### Our model ###### 
beta = c(1,0.1,0.1)
N = 1000
K = length(beta)
n_s = 2
k_s = 2
data_s = c(5, 2.5,2.5,4)
mu = c(0,0)
S = matrix(data_s, n_s, k_s)
X = matrix(NA,N,K)
X[,1] = 1
X[,2:3] = mvrnorm(N,mu,S) #simulate from Multivariate Normal Distribution
X

epsilon = rnorm(N,0,10)

Y = X%*%beta + epsilon

###### OLS_estimator ######

beta_hat =solve(t(X)%*%X)%*%t(X)%*%Y #solve inverts matrix
beta_hat

############F-statistic#######################

##F_value_data##
#H_0: all beta zero
#H_1: beta_1 or beta_2 unequal to zero

Y_mean = (1/N)*sum(Y) #mean of population
Y_hat = X%*%beta_hat

#[only for homoskedastic OLS standard errors!]
#F= [sum((y_hat-y_mean)^2)/(k-1)]/[sum((y-y_mean)^2)/(n-k)]
#more general:
#F = [(SSR_restricted - SSR_unrestricted)/q] / [SSR_unrestricted/(n - k_unrestricted - 1)]
#with SSR meaning sum of squared residuals from the (un)restricted regression
#k_unrestricted=number of regressors and q=number of restrictions under Null
#alternative with R^2: F = [(R^2_unrestricted - R^2_restricted)/q] / [(1 - R^2_unrestricted)(n - k_unrestricted - 1)]

f_emp_numerator = sum((Y_hat - Y_mean)^2/(K-1))
f_emp_denominator = sum((Y - Y_mean)^2/(N-K))             
f_emp= f_emp_numerator/f_emp_denominator
f_emp


##F_critical_value##

alpha = 0.05
df1 = K-1
df2 = N-K
f_crit = qf(1-alpha, df1, df2, lower.tail = TRUE)
f_crit

##F_test##
#comparing explained variance to unexplained, F value increasing along with high R^2 and sample size

if (f_emp > f_crit) {
    print("reject null hypothesis, at least one beta coefficient has a significant impact")
  } else {
    print("cannot reject null hypothesis, no beta coefficient has a significant impact")
  }
###F-test hypothesis is more likely to be rejected because the restrictions (all betas are equal to zero) are strong

############t-statistic#######################

##t_value_data##

epsilon_hat = Y - Y_hat
beta_hyp = 0


est_epsilon_var = (1/N)*(sum(epsilon_hat^2)) 


se_coef = c()
for(i in 1:length(beta)){
  se_coef[i] = sqrt(est_epsilon_var*(solve(t(X)%*%X)[i,i]))
}
se_coef


t_value_beta = c()
for (i in 1:length(beta)){
  t_value_beta[i] = abs((beta_hat[i]- beta_hyp)/ se_coef[i])
  i = i+1
}
t_value_beta

##t_critical_value##

alpha = 0.05
df = N-K
t_crit = abs(qt(alpha/2, df,lower.tail = FALSE))
t_crit


##t_test##

for(i in 1: length(beta)){
  if (t_value_beta[i] > t_crit) {
    print("reject null hypothesis")
  } else {
    print("cannot reject null hypothesis")
  } 
  #i = i+1 not needed in for-loop
}

#comment on for-loops: if i should increase by 2,3,... loop over vector

############## Calculating coverage probabilites ###########

coverage_probability_beta_2 = function(signflevel, beta_hat_test, se_coef_test){
  
  beta_test = c(1,0.1,0.1)
  alpha = 1 - signflevel
  t_crit = abs(qt(alpha/2, df,lower.tail = FALSE))
  
  confidence_beta_2_lower = beta_hat_test[2]-se_coef_test[2]*t_crit
  confidence_beta_2_upper = beta_hat_test[2]+se_coef_test[2]*t_crit

  
  if ( confidence_beta_2_upper >= beta_test[2] && confidence_beta_2_lower<= beta_test[2]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
  
coverage_probability_beta_3 = function(signflevel, beta_hat_test, se_coef_test){
    
    beta_test = c(1,0.1,0.1)
    alpha = 1 - signflevel
    t_crit = abs(qt(alpha/2, df,lower.tail = FALSE))
    
    confidence_beta_3_lower = beta_hat_test[3]-se_coef_test[3]*t_crit
    confidence_beta_3_upper = beta_hat_test[3]+se_coef_test[3]*t_crit
    
    
    if ( confidence_beta_3_upper >= beta_test[3] && confidence_beta_3_lower<= beta_test[3]) {
      return(TRUE)
    } else {
      return(FALSE)
    }
}

#better: write general function for coverage probability, not separate for each beta coefficient

N = 1000
count_beta_2 = c()
count_beta_3 = c()
for(i in 1:N){
  
  X_test = matrix(NA,N,3)
  X_test[,1] = 1
  X_test[,2:3] = mvrnorm(N,mu,S)
  epsilon_test = rnorm(N,0,10)
  beta_test = c(1,0.1,0.1)
  
  Y_test = X_test%*%beta + epsilon_test
  
  beta_hat_test =solve(t(X_test)%*%X_test)%*%t(X_test)%*%Y_test
  
  Y_hat_test = X_test%*%beta_hat_test
  
  
  epsilon_hat_test = Y_test - Y_hat_test 
  est_epsilon_var_test = (1/N)*(sum(epsilon_hat_test^2)) 
  
  se_coef_test = c()
  for(k in 1:length(beta)){
    se_coef_test[k] = sqrt(est_epsilon_var_test*(solve(t(X_test)%*%X_test)[k,k]))
  }
  se_coef_test
  
  count_beta_2[i]=(coverage_probability_beta_2(0.95,beta_hat_test,se_coef_test))
  count_beta_3[i]=(coverage_probability_beta_3(0.95,beta_hat_test,se_coef_test))
  
}

coverage_prob_beta_2 = sum(count_beta_2)/N
coverage_prob_beta_3 = sum(count_beta_3)/N
coverage_prob_beta_2
coverage_prob_beta_3



#####simulation study


statistics = function(X_func, beta_func,N, alpha, eps_func){
  
  ###OLS_estimator###
  
  Y = X_func%*%beta_func + eps_func
  
  beta_func_hat =solve(t(X_func)%*%X_func)%*%t(X_func)%*%Y
  
  #F-test
  Y_func_mean = 1/N*sum(Y)
  Y_func_hat = X%*%beta_func_hat
  
  f_func_emp = (sum((Y_hat - Y_mean)^2/(ncol(X_func)-1)))/(sum((Y - Y_mean)^2/(nrow(X_func)-ncol(X_func))))
  
  f_func_crit = qf(1-alpha, ncol(X_func)-1, nrow(X_func)-ncol(X_func), lower.tail = TRUE)
  
  #t_test
  eps_func_hat = Y - Y_func_hat
  
  eps_func_var = 1/N*(sum(eps_func_hat^2))
  
  se_func_coef = c()
  for(i in 1:length(beta_func_hat)){
    se_func_coef[i] = sqrt(eps_func_var*(solve(t(X_func)%*%X_func)[i,i]))
  }
  
  t_func_beta = c()
  for (i in 1:length(beta_func_hat)){
    t_func_beta[i] = (beta_func_hat[i])/ se_func_coef[i]
  }
  
  t_func_crit = abs(qt(alpha/2, nrow(X_func)-ncol(X_func),lower.tail = FALSE))
  
  
  
  if (f_func_emp > f_func_crit && t_func_beta[2] > t_func_crit && t_func_beta[3] > t_func_crit) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}

####running simulation with a normally distributed epsilon 

N_simultation = 1000
simultation = c()
for(i in 1:N_simultation){
  
  X_simultation = matrix(NA,N,K)
  X_simultation[,1] = 1
  X_simultation[,2:3] = mvrnorm(N,mu,S)
  
  epsilon_simulation = rnorm(N,0,10)
  
  Y_simulation = X_simultation%*%beta + epsilon_simulation
  
  beta_hat_simulation =solve(t(X_simultation)%*%X_simultation)%*%t(X_simultation)%*%Y_simulation
  
  simultation[i] = statistics(X_simultation, beta_hat_simulation,1000, 0.05, epsilon_simulation)
  
  i = i+1
}
simultation
sum(simultation)

### Covariates are drawn from one distribution, so multicollinearity might occur. Due to that, the standard deviations of the beta coefficients are high and the null hypothesis of the t-test is unlikely to be rejected.  


####modeling heteroskedasticity 

#These standard errors are heteroskedastic, but not dependent on x

eps_hetero = c()
for (i in 1:N){
  eps_hetero[i] = rnorm(1,0,10+i/10-1)
  i = i+1
}

eps_hetero_random = c()
for (i in 1:N){
  eps_hetero_random[i] = rnorm(1,0,runif(1, min=10, max=30))
  i = i+1
}

plot(epsilon,#plot the results of all samples
     main="homoskedasticity",
     xlab="N",ylab="epsilon",pch=16,
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

plot(eps_hetero,#plot the results of all samples
     main="heteroskedasticity",
     xlab="N",ylab="epsilon",pch=16,
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

plot(eps_hetero_random,#plot the results of all samples
     main="heteroskedasticity",
     xlab="N",ylab="epsilon",pch=16,
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)




N_simulation_hetero = 1000

simulation2 = function(N_sim, eps){
  simulation_hetero = c()
  for(i in 1:N_sim){
    
    X_simulation = matrix(NA,N,K)
    X_simulation[,1] = 1
    X_simulation[,2:3] = mvrnorm(N,mu,S)
    
    
    Y_simulation = X_simulation%*%beta + eps
    
    beta_hat_simulation =solve(t(X_simulation)%*%X_simulation)%*%t(X_simulation)%*%Y_simulation
    
    simulation_hetero[i] = statistics(X_simulation, beta_hat_simulation,1000, 0.05, eps)
  }
  simulation_hetero
  sum(simulation_hetero)
}

simulation2(1000, eps_hetero)


