## Set global parameters ##
n=100
sim=100
var=10
#####Example dgp#########
x<-rnorm(n)
eps<-rnorm(n,0,sqrt(var))
beta<-c(1,0.1)
X<-cbind(rep(1,n),x)
y=X%*%beta+eps
###This function makes a partition of the data and then 
MSE_CV<-function(y,X,k,n){
  ######Divide the data into folds.
  ######Train the model
  ######Validate on the remaining fold
  ######Repeat for all k folds
  ######Estimate the average cross-validation score across folds for both the simple linear model and the 
  ######model with squared term
  ######pick the model with lowest average cross-validation criterion
  index=1:n
  x=n/k
  mse<-c()
  mse_sq=c()
  #####Generate the folds
  for(i in 1:k)
  {
    nam <- paste("folds_", i, sep = "")
    value=sample(index,x)
    assign(nam, value)
    namx=paste("X_", i, sep = "")
    assign(namx, X[value,])
    namy=paste("y_", i, sep = "")
    assign(namy, y[value])
    index=index[-value]
  }
  ###fit model###
  X_sq<-cbind(X,X[,2]^2)
  ###for i=1 means fold 1 is the hold out fold##
  for(i in 1:k)
      {
    hold_out_fold=ls(pattern = "^folds_.$")[i]
    value=get(hold_out_fold)
    X_train=X[-value,]
    y_train=y[-value]
    X_sq_train<-X_sq[-value,]
    beta_hatsq<-solve(t(X_sq_train)%*%X_sq_train)%*%t(X_sq_train)%*%y_train
    beta_hat<-solve(t(X_train)%*%X_train)%*%t( X_train)%*%y_train
    y_hat=X[value,]%*%beta_hat
    y_hat_sq=X_sq[value,]%*%beta_hatsq
    mse[i]=mean(y_hat-y[value])**2
    mse_sq[i]=mean(y_hat_sq-y[value])**2
  }
  m.cv=c(mean(mse),mean(mse_sq))
return(m.cv)
}
######################################################
mod_prob=c()
######For k=2##########################################
k=2
MSE_CV(y,X,k,n)
MSE<-matrix(NA,sim,2)
winner=c()### =1 if the correct (linear) model is selected
for(i in 1:sim)
{
  x<-rnorm(n)
  eps<-rnorm(n,0,sqrt(var))
  X<-cbind(rep(1,n),x)
  y=X%*%beta+eps
MSE[i,]=MSE_CV(y,X,k,n)
winner[i]=MSE[i,1]<MSE[i,2]
}
plot(sort(winner),ylab="winner==1")
mod_prob[1]=mean(winner)
######For k=5##########################################
k=5
MSE_CV(y,X,k,n)
MSE<-matrix(NA,sim,2)
winner=c()###gives yes if the correct model is selected
for(i in 1:sim)
{
  x<-rnorm(n)
  eps<-rnorm(n,0,sqrt(var))
  beta<-c(1,0.1)
  X<-cbind(rep(1,n),x)
  y=X%*%beta+eps
  MSE[i,]=MSE_CV(y,X,k,n)
  winner[i]=MSE[i,1]==min(MSE[i,])
}
lines(sort(winner),col="blue",type="p")
mod_prob[2]=mean(winner)
mod_prob






######For k=4##########################################
k=4
MSE_CV(y,X,k,n)
MSE<-matrix(NA,sim,2)
winner=c()###gives yes if the correct model is selected
for(i in 1:sim)
{
  x<-rnorm(n)
  eps<-rnorm(n,0,sqrt(var))
  beta<-c(1,0.1)
  X<-cbind(rep(1,n),x)
  y=X%*%beta+eps
  MSE[i,]=MSE_CV(y,X,k,n)
  winner[i]=MSE[i,1]==min(MSE[i,])
}
lines(sort(winner),col="red",type="p")
mod_prob[2]=mean(winner)











##Data Generating Process## 
n=100
x<-rnorm(n)
eps<-rnorm(n,0,sqrt(10))
beta<-c(1,0.1)
X<-cbind(rep(1,n),x)
y=X%*%beta+eps

######
X_1<-cbind(X,x^2)
beta_hat1<-solve(t(X_1)%*%X_1)%*%t(X_1)%*%y
beta_hat<-solve(t(X)%*%X)%*%t(X)%*%y
######
x.grid=seq(min(x),max(x),length=100)
plot(x,y)
abline(a=beta_hat[1],b=beta_hat[2],col="red")
lines(x.grid,beta_hat1[1]+x.grid*beta_hat1[2]+x.grid^2*beta_hat1[3],type="l",col="blue")
##
