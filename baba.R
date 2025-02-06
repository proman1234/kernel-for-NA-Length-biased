# Load the necessary package 
# install.packages("quantmod")
library(quantmod)
# Download historical data for Alibaba (BABA) from Yahoo
getSymbols("BABA", src = "yahoo", from = "2023-01-03", to = "2024-12-31") 
# the closing price
stock_data <- Cl(BABA)
#Calculate returns
stock_returns <- dailyReturn(stock_data)
# # Plot the stock returns as a time series plot
n <- length(stock_data)
plot(1:n, stock_returns, ylab= "The returns of BABA stock prices", type = "b")

# Plot the autocorrelation function
acf(stock_returns, main="ACF of Stock Returns")


F1 <- function(h,t,e){
  res <- mean(pnorm((t-e)/h, mean = 0, sd = 1))
  return(res)
}


F2 <- function(t,e)
{
  res <- sum(as.numeric(e <= t))/length(e)
  return(res)
}

F3 <- function(t,e){
  sigmax <- sd(e)
  n <- length(e)
  i <- seq(1,n,by=1)
  h <- 2*sigmax*i^(-1/4)
  res <- mean(pnorm((t-e)/h, mean = 0, sd = 1))
}

f1 <- function(t,e){
  n <- length(e)
  i <- seq(1,n,by=1)
  hi <- (log(i)/i)^(1/6)
  res <- mean(qnorm((t-e)/hi,mean=0,sd=1)/hi)
  return(res)
}

f2 <- function(t,e){
  res <- sum(as.numeric(e <= t))/length(e)
}

f3 <- function(h,t,e){
  res = 1/(sqrt(2*pi)*h)*mean(exp(-(t-e)^2/(2*h^2))) 
  return(res)
}

f4 <- function(h, t, e){
  res = length(which( (t-h) < e & e <= (t+h)))/(2*length(e)*h)
}

gbsg <- stock_returns
plot(density(gbsg))

 hn = sd(gbsg)*(log(n)/n)^(1/4) # bandwidths
# hn = (log(n)/n)^(1/4)
# hn = n1^(-1/5)
m <- 50
t = seq(-0.1, 0.1, length = m)

fhat <- fn <- matrix(NA, m, 1)
for (i in 1:m){
  fhat[i] <- f3(h=hn, t=t[i], e = gbsg)
}

for (i in 1:m){
  fn[i] <- f4(h=hn, t=t[i], e = gbsg)
}

Fn <- Fhat <- matrix(NA, m, 1)
for (i in 1:m){
  Fn[i] <- F2(t=t[i], e = gbsg)
}

for (i in 1:m){
  Fhat[i] <- F1(h=hn, t=t[i], e = gbsg)
}


### the distribution function ###
# F.t <- pnorm(t, mean = mean(gbsg), sd = sd(gbsg), log = FALSE)
# lines(t, F.t, lwd=3, col="black",lty=1)
par(mai=c(0.7,.7,.4,.4),cex=0.8)
plot(t, Fhat, xlab="x", ylab= "distribution functions", type = "n", ylim = c(0, 1.01),lwd=2)
lines(t, Fn, col="blue", lty=3, lwd=2)
lines(t, Fhat, col="red",lty=4,lwd=2)
curve(pnorm(x, mean = mean(gbsg), sd = sd(gbsg)), add = T, lwd=3, col= "black", lty=1)
abline(h=1)
legend("bottomright", c( expression(F[n](x)), expression(hat(F)[n](x)), expression(F(x))), 
       col = c("blue","red","black"), lty = c(3,4,1), lwd=2)


### the density function ###
#f.t <- dnorm(t, mean = mean(gbsg), sd = sd(gbsg), log = FALSE)
# lines(t, f.t, lwd=3, col="black",lty=1)
par(mai=c(0.7,.7,.4,.4), cex=0.8)
plot(t, fhat, xlab="x", ylab= "density functions", type = "n", ylim = c(0, 21),lwd=2)
lines(t, fn, col="blue", lty=3, lwd=2)
lines(t, fhat, col="red",lty=4,lwd=2)
abline(v=0)
curve(dnorm(x, mean = mean(gbsg), sd = sd(gbsg)), add = T, lwd=3, col= "black", lty=1)
# title("(d)  Age > 65")
legend("topright", c( expression(f[n](x)), expression(hat(f)[n](x)), expression(f(x))), 
       col = c("blue","red","black"), lty = c(3,4,1), lwd=2)


### hazard function ########
rhat = rn = matrix(NA, m, 1)
rn <- fn/(1 - Fn)
rhat <- fhat/(1 - Fhat)
par(mai = c(0.7,0.7,0.5,0.5),cex=0.8)
plot(t, rn, type = "n", xlab = "x", ylab = "harzard functions")
# plot(rn)
lines(t,rn,col="orange",lty=2,lwd=2)
lines(t,rhat,col="green",lty=5,lwd=2)
curve(dnorm(x, mean = mean(gbsg), sd = sd(gbsg))/(1-pnorm(x, mean = mean(gbsg), sd = sd(gbsg))),add = T, lwd=3, col= "black", lty=1)
legend("topleft", c(expression(r[n](x)), expression(hat(r)[n](x))), 
       col = c("orange","green"), lty = c(2,5), lwd=2)

