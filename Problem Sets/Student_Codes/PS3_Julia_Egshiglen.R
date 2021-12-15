library(ggplot2)
library(RColorBrewer)
library(ggpubr)


# Normally distributed data -----------------------------------------------

n <- 100
set.seed(100)
X_sample <- rnorm(n, 0, sqrt(1))
X_sample

h <- 0.5
kernel_est <- function(x, h, X_i){
  summe <- 0
  for(i in 1:length(X_i)){
    local_x <- (x - X_i[i])/h
    K_x <- 3/4*(1- local_x ** 2) * (abs(local_x) <= 1)
    summe = summe + K_x
  }
  return(summe/(h* length(X_i)))
}

kernel_est(X_sample[2], h, X_sample)

kernel_vector <- function(daten, bw){
  kern <- c()
  for(i in 1:length(daten)){
    kern[i] <- kernel_est(daten[i], bw, daten)
  }
  return(kern)
}
kernel_vector(X_sample, 0.5)

#Plot for h=0.5:
data_norm <- as.data.frame(cbind(X_sample, kernel_vector(X_sample, h)))

ggplot(data_norm, aes(x = X_sample, y = V2)) + geom_point() + 
  geom_line() + theme_bw() + xlab('Data') + ylab('f_hat')

###Comparison between different bandwidths####

data_comp <- as.data.frame(cbind(X_sample, kernel_vector(X_sample, 0.1), 
                                 kernel_vector(X_sample, 0.3), kernel_vector(X_sample, 0.5),
                                 kernel_vector(X_sample, 0.7), kernel_vector(X_sample, 1)))

colnames(data_comp) <- c("X_sample", "h01", "h03", "h05", "h07", "h1")

#Plot:

true_dnorm100 <- data.frame(true_x = seq(-2,2,le=100), f=dnorm(seq(-2, 2, le=100)))

colors6 <- c("h01" = "#1B9E77", "h03" = "#D95F02", "h05" = "#7570B3", 
             "h07" = "#E7298A", "h1" = "#E6AB02", "reference" = "#666666")
plot_bw <- ggplot() + 
  geom_line(data=data_comp, aes(x = X_sample, y=h01, color = "h01")) + 
  geom_line(data=data_comp, aes(x = X_sample, y=h03, color = "h03")) +
  geom_line(data=data_comp, aes(x = X_sample, y=h05, color = "h05")) + 
  geom_line(data=data_comp, aes(x = X_sample, y=h07, color = "h07")) + 
  geom_line(data=data_comp, aes(x = X_sample, y=h1, color = "h1")) + 
  geom_line(data=true_dnorm100, aes(x = true_x, y=f, color = "reference"), size=1.2, alpha=0.4) +
  theme_bw() + labs(x = "Data", y = "f_hat", color = "Bandwidths") + 
  scale_color_manual(values = colors6)
plot_bw

###Comparison between different sample sizes####

set.seed(124)
k <- 500
X_k <- rnorm(k, 0, sqrt(1))

density_plot <- function(data1, data2, h, axis_active = FALSE){
  data_local1 <- data.frame(X1 = data1, Y1 = kernel_vector(data1, h)) 
  data_local2 <- data.frame(X2 = data2, Y2 = kernel_vector(data2, h)) 
  color2 <- c("n=100" = "#1B9E77", "n=500" = "#D95F02")
  g_plot <- ggplot() + geom_line(data=data_local1, aes(x = X1, y = Y1, color = "n=100")) +
    geom_line(data=data_local2, aes(x = X2, y = Y2, color = "n=500")) +
    theme_bw() + labs(y = "f_hat", color = "Sample sizes", subtitle=sprintf("h = %s", h)) + scale_color_manual(values = color2) +
    scale_y_continuous(limits= c(0,0.8),breaks=seq(0,0.8,0.2)) + theme(axis.title.x = element_blank())
  if(!axis_active){
    g_plot <- g_plot + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title.y = element_blank())
  }
  return(g_plot)
}

plot_h01 <- density_plot(X_sample, X_k, 0.1, TRUE)
plot_h05 <- density_plot(X_sample, X_k, 0.5)
plot_h1 <- density_plot(X_sample, X_k, 1, TRUE)
plot_h2 <- density_plot(X_sample, X_k, 2)

plot_size <- ggarrange(plot_h01, plot_h05, plot_h1, plot_h2, ncol = 2, nrow = 2, 
                       align = "v", common.legend = TRUE, legend = "bottom")
plot_size
