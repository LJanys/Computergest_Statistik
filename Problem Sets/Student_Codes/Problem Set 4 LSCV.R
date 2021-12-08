convo.k<-function(x){
  return(exp(-x^2/2)/sqrt(2*pi))
} #Convolution Kernel
gauss.k<-function(x){
  return(exp(-x^2/4)/sqrt(4*pi))
} #Gaussian Kernel

n=100

#1
h.grid<-c(seq(0.01,0.5, by=0.01)) #50 equidistant points
LSCV<-function(h.grid){
  x<-rnorm(n = n, mean = 0, sd = 1)
  LSCV.val<-rep(NA, 50)
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
  return(0.01*(match(min(LSCV.val),LSCV.val)))
  #return(LSCV.val) #for the very last plot
}
LSCV.simplified<-function(h.grid){
  x<-rnorm(n = n, mean = 0, sd = 1)
  LSCV.val<-rep(NA, 50)
  for(k in 1:length(h.grid)){
    sum1=0
    #sum2=0
    for(i in 1:n){
      for(j in 1:n){
        sum1 = sum1 + (
                  ((1/(n^2*h.grid[k]))*convo.k((x[i]-x[j])/h.grid[k]))
                - ((2/(n*(n-1)*h.grid[k]))*gauss.k((x[i]-x[j])/h.grid[k]))
                    )
      }
    }
    LSCV.val[k] = sum1 + (((2/pi)^0.5)/((n-1)*h.grid[k]))
  }
  return(0.01*(match(min(LSCV.val),LSCV.val)))
}

#2
LSCV(h.grid = h.grid)
#LSCV.simplified(h.grid = h.grid)

#system.time(LSCV(h.grid = h.grid))
#system.time(LSCV.simplified(h.grid = h.grid))

#3
LSCV.simulation <- function(nsim){
  h.min.all <- c(rep(NA, nsim))
  for (s in 1:nsim){
    h.min.all[s] <- LSCV(h.grid)
  }
  return(h.min.all)  
}
#LSCV.simulation.simplified <- function(nsim){
  h.min.all <- c(rep(NA, nsim))
  for (s in 1:nsim){
    h.min.all[s] <- LSCV.simplified(h.grid)
  }
  return(h.min.all)  
}

h.min.all = LSCV.simulation(50)
#system.time(LSCV.simulation(20))
#system.time(LSCV.simulation.simplified(20))

#4
mean(h.min.all)
sd(h.min.all)
h.min = min(h.min.all)
h.max = max(h.min.all)

x.grid<-c(seq(-2,2,by=0.1))
NP.density.estimation<-function(n, x, h){
  set.seed(1207)
  X = rnorm(n = n,mean = 0,sd = 1)
  
  f.hat<-c(rep(NA,length(x)))
  for(i in 1:length(x)){
    sum=0
    for(j in 1:n){
      sum = sum+gauss.k((x[i]-X[j])/h)
    }
    f.hat[i]=(1/(n*h))*sum
  }
  return(list(f.hat = f.hat))
}

h.min.data <- NP.density.estimation(n = 100,x = x.grid,h = h.min)$f.hat
h.max.data <- NP.density.estimation(n = 100,x = x.grid,h = h.max)$f.hat

plot(x = x.grid,y = h.min.data,type = "l",lwd = 1.3, col = "red",
     xlab = "x",ylab = "y",ylim = c(0,0.6),
     main = paste("Density Functions","\nn =", n), cex.main = 1)
lines(x = x.grid,y = h.max.data,type = "l",lwd = 1.3, col = "blue")
legend("topleft",legend = c(paste("min of h.min.all:", h.min),paste("max of h.min.all:", h.max) ),
       col = c("red","blue"),lwd = 1.3,lty=1:1,cex=0.7)



################
plot(x = h.grid,y = LSCV(h.grid))
