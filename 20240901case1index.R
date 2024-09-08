rm(list = ls(all = TRUE))
library(MASS)
N = 250 
M <- 1 # 模拟循环次数
n = 200 # 样本量
beta = 0.1 # 误差协方差

# 生成协方差
Cov = NULL
mu = rep(0,N)
####### 协方差阵生成函数
Covf <- function(rho,n)
{
  Cov = matrix(NA,n,n)
  for (i in 1:n){
    for (j in 1:n){
      Cov[i,j] <- -rho^abs(i-j)
    }
  }
  res <- Cov
  return(res)
}
Cov = 2*diag(N) + Covf(beta,N)
#生成数据

E <- t(mvrnorm(M, mu, Cov))
X <- exp(E)
X <- matrix(X, N, M)
# X

## generate y
Y <- matrix(NA, n, M)
for (j in 1:M) {
  l <- 0; y <- c()
  while (l < n) {
    y.temp <- sample(X[,j], 1)
    u <- runif(1, 0, 10) ## distribution of Ti
    if (u <= y.temp) {y <- c(y, y.temp)} else {y <- y}
    l <- length(y)
  }
  Y[,j] <- y  
}

k.q = function(u)     ## the quartic kernel function 
{
  res = 15*(1 - u^2)^2/16*(abs(u)<=1)
  return(res)
}

k.e = function(u) ## epanechnikov kernel function
{
  res <- 3*(1-u^2)/4*(abs(u)<=1)
  return(res)
}

K = function(x) ## the kernel distribution function 
{
  res <- integrate(k.q, lower = -1.5, upper = x)$value
  return(res)
}


### K(c(-1,-12))
K(0)
integrate(k.q, lower = -Inf, upper = 0.2)$value

## the true distribution function
F.MA <- function(x){
  res <- pnorm(log(x), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  return(res)
}
F.MA(2)

Fns <- function(x,y)
{
  #  y <- Y[1,]
  #   x=3.9
  n <- length(y)
  hn <- IQR(y, na.rm = FALSE, type = 7)*n^{-1/3}
  # hn
  # mu <- 1/mean(1/y)
  K.temp <- rep(NA,n)
  for (i in 1:n) {
    #  i=4
    K.temp[i] <- K((x-y[i])/hn)
  }
  # K((x-y)/hn); K(-0.988)
  K.temp
  res <- mean(K.temp/y)/mean(1/y)
  return(res)
}

# Fns(x=1, y=Y[1,])

Fnc <- function(x,y)
{
  res <- mean(as.numeric(y <= x)/y)/mean(1/y)
  return(res) # return(j) # return() 不能直接使用，只能在编写函数时使用
}

t = seq(0, 8, by = 0.2)
m = length(t)
F.t <- c()
for (i in 1:m) {
  F.t[i] <- F.MA(t[i])
}

Fnshat = Fnchat = matrix(NA, m, M)
for (j in 1:M){
  #  j=3
  # x=4
  Fnshat.temp = c()
  for (i in 1:m){
    # i=30
    Fnshat.temp[i] <- Fns(x = t[i], y = Y[,j])
  }
  Fns(x=3.9, y = Y[,j])
  Fnshat[,j] <- Fnshat.temp
  print(j)
}

# Fnshat[,224]
Fhatm <- apply(Fnshat, 1, mean)
# plot(t, Fhatm, ylim = c(0,1))
# lines(t, F.t, lwd=3, col="black",lty=1)

for (j in 1:M){
  Fnc.temp = c()
  for (i in 1:m){
    Fnc.temp[i] <- Fnc(x=t[i], y=Y[,j])
  }
  Fnchat[,j] <- Fnc.temp
  print(j)
}

# par(mar=c(1,1,1,1))   #设置边距参数
Fcm <- apply(Fnchat, 1, mean)

#### compute index D(F) ######
DF = DFhat = NULL
DFs.temp = DFc.temp= matrix(NA,N,1)
for (i in 1:N) {
  DFs.temp[i] <- max(Fnshat[,i] - F.t)
  DFc.temp[i] <- max(Fnchat[,i] - F.t)
}
DFs <- mean(DFs.temp)
DFc <- mean(DFc.temp)
Dratio = DFc/DFs

#### compute index MISE ######
ISE1 = NULL
for (i in 1:N){
  ISE1[i] <- sum(( Fnshat[,i] - F.t )^2)*0.1 }
MISE.Fns <- mean(ISE1)

ISE2 = NULL
for (i in 1:N){
  ISE2[i] <- sum(( Fnchat[,i]- F.t )^2)*0.1 }
MISE.Fnc <- mean(ISE2)
Mratio <- MISE.Fnc/MISE.Fns

index = cat(DFs, Dratio, MISE.Fns, Mratio)
