#Please note:
#We made the two parts of the Problem in two different files, where
#we used the same simple variables.
#So please keep in mind which code/function you executed last.

gauss.k<-function(x){
  return(exp(-x^2/4)/sqrt(4*pi))
}
convo.k<-function(x){
  return(exp(-x^2/2)/sqrt(2*pi))
}
LSCV<-function(h.grid, x){
  LSCV.val<-rep(NA, 30)
  for(k in 1:length(h.grid)){
    sum1=0
    sum2=0
    for(i in 1:n){
      for(j in 1:n){
        sum1 = sum1 + convo.k((x[i]-x[j])/h.grid[k])
      }
    }
    for(i in 1:n){
      for(j in 1:n){
        if(j != i){
          sum2 = sum2 + gauss.k((x[i]-x[j])/h.grid[k])
        }else{
          sum2 = sum2
        }
      }
    }
    LSCV.val[k]<-(1/(h.grid[k]*n^2)*sum1)-(2/(n*(n-1)*h.grid[k])*sum2)
  }
  return(0.1 + 0.1*(match(min(LSCV.val),LSCV.val)-1))
}
NP.density.estimation<-function(x, h){
  f.hat<-c(rep(NA,n))
  for(i in 1:n){
    sum=0
    for(j in 1:n){
      sum = sum+gauss.k((x[i]-x[j])/h)
    }
    f.hat[i]=(1/(n*h))*sum
  }
  return(list(f.hat = f.hat))
}
NP.density.estimation.grid<-function(x.grid, x, h){
  f.hat<-c(rep(NA,length(x.grid)))
  for(i in 1:length(x.grid)){
    sum=0
    for(j in 1:n){
      sum = sum+gauss.k((x.grid[i]-x[j])/h)
    }
    f.hat[i]=(1/(n*h))*sum
  }
  return(list(f.hat = f.hat))
}
###############Silverman's rule of thumb###############

n = 1000
h.grid = c(seq(0.1, 3, by=0.1))
my2 = c(20,10,1)

for(i in 1:length(my2)){
  x <- c(rnorm(n = n/2,mean = 0,sd = 1), rnorm(n = n/2,mean = my2[i],sd = 1))
  print(sd(x))
  
  h.LSCV.min = LSCV(h.grid = h.grid,x = x)
  data = NP.density.estimation(x = x,h = h.LSCV.min)$f.hat
  
  
  h.s = ((4*sd(x)^5)/(3*n))^(1/5)
  if(my2[i]==20){
    x.grid = c(seq(-3,23,by = 26/30))
    data.s = NP.density.estimation.grid(x.grid = x.grid,x = x,h = h.s)$f.hat
  }
  else if(my2[i]==10){
    x.grid = c(seq(-3,13,by = 16/30))
    data.s = NP.density.estimation.grid(x.grid = x.grid,x = x,h = h.s)$f.hat
  }
  else{
    x.grid = c(seq(-3,3, by = 6/30))
    data.s = NP.density.estimation.grid(x.grid = x.grid,x = x,h = h.s)$f.hat
  }
  
  plot(x = x,y = data,type="n",ylab="y", xlab="x",ylim=c(0,0.38),
       main=paste("Density Functions","\nn:", n,"   mu2:", my2[i]), cex.main=1)
  #ordering the points to plot a line:
  lines(x[order(x)], data[order(x)],lwd = 1.3)
  #adding the density function with h.s
  lines(x = x.grid,y = data.s, col = "blue",lwd = 1.3)
  legend("topleft",legend = c(paste("h.s:", round(h.s,3)),paste("h.lscv:", h.LSCV.min) ),
         col = c("blue","black"),lwd = 1.3,lty=1:1,cex=0.8)
}

