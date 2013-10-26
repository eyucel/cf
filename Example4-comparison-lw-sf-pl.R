#####################################################################################
#
# Example 4: Local level model with parameter learning
#
# Lopes and Tsay (2011)
# Particle Filters and Bayesian Inference in Financial Econometrics.
# Journal of Forecasting.
# 
#####################################################################################
#
# Comparison of three particle filters for the local level model 
# when parameter learning is also taken into account.
#
# For t=1,...,n
#
#    y[t]|x[t]   ~ N(x[t],sig2)
#    x[t]|x[t-1] ~ N(alpha+beta*x[t-1],tau2)
#
# and
#
#    x[0]       ~ N(m0,C0)
#    sig2       ~ IG(a0,A0)
#    alpha|tau2 ~ N(b0[1],tau2*B0[1])
#    beta|tau2  ~ N(b0[2],tau2*B0[2])
#    tau2       ~ IG(nu0/2,nu0*tau20/2)
#
# Known quantities: (m0,C0,a0,A0,b0,B0,nu0,tau20)
#
# Particle filters: 
#    Liu and West filter (Liu and West, 2001)
#    Storvik's filter (Storvik, 2002)
#    Particle learning (Carvalho et al., 2010)
#
######################################################################################
#
# Author : Hedibert Freitas Lopes
#          Associate Professor of Econometrics and Statistics
#          The University of Chicago Booth School of Business
#          5807 South Woodlawn Avenue
#          Chicago, Illinois, 60637
#          Email : hlopes@ChicagoBooth.edu
#          URL: http://faculty.chicagobooth.edu/hedibert.lopes
#
#####################################################################################
rm(list=ls())

getindex = function(x){x[1+x[1]]}

LW = function(y,alphas,betas,tau2s,sig2s,xs,a){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,5,3))
  h2 = 1-a^2
  pars = cbind(alphas,betas,log(tau2s),log(sig2s))
  for (t in 1:n){
    # Resampling
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,4,byrow=T)
    weight  = dnorm(y[t],alphas+betas*xs,sqrt(sig2s+tau2s))
    k       = sample(1:N,size=N,replace=T,prob=weight)
    pars    = ms[k,] + matrix(rnorm(4*N),N,4)%*%chol(h2*vpar)
    alphas  = pars[,1]
    betas   = pars[,2]
    tau2s   = exp(pars[,3])
    sig2s   = exp(pars[,4])
    # Propagating
    vars    = 1/(1/sig2s+1/tau2s)
    means   = vars*(y[t]/sig2s+(alphas+betas*xs[k])/tau2s)
    xs      = rnorm(N,means,sqrt(vars))
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(betas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(tau2s,c(0.025,0.5,0.975))
    quants[t,4,] = quantile(sig2s,c(0.025,0.5,0.975))  
    quants[t,5,] = quantile(exp(xs/2),c(0.025,0.5,0.975))  
  }
  return(quants)
}

Storvik = function(y,alphas,betas,tau2s,sig2s,xs){
  n      = length(y)
  N      = length(xs)
  quants = array(0,c(n,5,3))
  s      = matrix(0,N,9)
  s[,1]  = 1.0/B0[1]
  s[,2]  = 0.0
  s[,3]  = 1.0/B0[2]
  s[,4]  = b0[1]/B0[1]
  s[,5]  = b0[2]/B0[2]
  s[,6]  = nu0
  s[,7]  = nu0*tau20 
  s[,8]  = n0
  s[,9]  = n0*sig20
  m      = s[,1]*s[,3]-s[,2]^2
  num1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
  num2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
  for (t in 1:n){
    # Propagating
    xs1    = xs
    mus    = alphas+betas*xs1
    vars   = 1/(1/sig2s+1/tau2s)
    sds    = sqrt(vars)
    means  = vars*(y[t]/sig2s + mus/tau2s)
    xs     = rnorm(N,means,sds)
    # Resampling
    stdevs = sqrt(sig2s+tau2s)
    weight = dnorm(y[t],mus,stdevs)
    k      = sample(1:N,size=N,replace=T,prob=weight)
    alphas = alphas[k]
    betas  = betas[k]
    tau2s  = tau2s[k]
    sig2s  = sig2s[k]
    mus    = mus[k]
    xs1    = xs1[k]
    xs     = xs[k]
    so     = s[k,]
    num1o  = num1[k]
    num2o  = num2[k]	  
    # Updating sufficient statistics
    s[,1]  = so[,1] + 1
    s[,2]  = so[,2] + xs1
    s[,3]  = so[,3] + xs1^2
    s[,4]  = so[,4] + xs
    s[,5]  = so[,5] + xs*xs1
    s[,6]  = so[,6] + 1
    s[,8]  = so[,8] + 1
    s[,9]  = so[,9] + (y[t]-xs)^2
    # Sampling parameters
    sig2s  = 1/rgamma(N,s[8]/2,s[,9]/2)
    m      = s[,1]*s[,3]-s[,2]^2
    num1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
    num2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
    s[,7]  = so[,7] + (xs-num1-num2*xs1)*xs + (num1o-num1)*so[,4]+(num2o-num2)*so[,5]
    tau2s  = 1/rgamma(N,s[,6]/2,s[,7]/2)
    std    = sqrt(tau2s/m)
    norm   = cbind(rnorm(N,0,std),rnorm(N,0,std))
    alphas = num1 + sqrt(s[,3])*norm[,1]
    betas  = num2 - s[,2]/sqrt(s[,3])*norm[,1]+sqrt(s[,1]-s[,2]^2/s[,3])*norm[,2]  
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(betas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(tau2s,c(0.025,0.5,0.975))  
    quants[t,4,] = quantile(sig2s,c(0.025,0.5,0.975))  
    quants[t,5,] = quantile(exp(xs/2),c(0.025,0.5,0.975))  
  }
  return(quants)
}

PL = function(y,alphas,betas,tau2s,sig2s,xs){
  n      = length(y)
  N      = length(xs)
  quants = array(0,c(n,5,3))
  s      = matrix(0,N,9)
  s[,1]  = 1.0/B0[1]
  s[,2]  = 0.0
  s[,3]  = 1.0/B0[2]
  s[,4]  = b0[1]/B0[1]
  s[,5]  = b0[2]/B0[2]
  s[,6]  = nu0
  s[,7]  = nu0*tau20 
  s[,8]  = n0
  s[,9]  = n0*sig20
  m      = s[,1]*s[,3]-s[,2]^2
  num1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
  num2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
  for (t in 1:n){
    # Resampling
    mus    = alphas+betas*xs
    stdevs = sqrt(sig2s+tau2s)
    weight = dnorm(y[t],mus,stdevs)
    k      = sample(1:N,size=N,replace=T,prob=weight)
    alphas = alphas[k]
    betas  = betas[k]
    tau2s  = tau2s[k]
    sig2s  = sig2s[k]
    mus    = mus[k]
    xs1    = xs[k]
    # Propagating
    vars   = 1/(1/sig2s+1/tau2s)
    sds    = sqrt(vars)
    means  = vars*(y[t]/sig2s + mus/tau2s)
    xs     = rnorm(N,means,sds)
    so     = s[k,]
    num1o  = num1[k]
    num2o  = num2[k]	  
    # Updating sufficient stats
    s[,1]  = so[,1] + 1
    s[,2]  = so[,2] + xs1
    s[,3]  = so[,3] + xs1^2
    s[,4]  = so[,4] + xs
    s[,5]  = so[,5] + xs*xs1
    s[,6]  = so[,6] + 1
    s[,8]  = so[,8] + 1
    s[,9]  = so[,9] + (y[t]-xs)^2
    # Sampling parameters
    sig2s  = 1/rgamma(N,s[8]/2,s[,9]/2)
    m      = s[,1]*s[,3]-s[,2]^2
    num1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
    num2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
    s[,7]  = so[,7] + (xs-num1-num2*xs1)*xs + (num1o-num1)*so[,4]+(num2o-num2)*so[,5]
    tau2s  = 1/rgamma(N,s[,6]/2,s[,7]/2)
    std    = sqrt(tau2s/m)
    norm   = cbind(rnorm(N,0,std),rnorm(N,0,std))
    alphas = num1 + sqrt(s[,3])*norm[,1]
    betas  = num2 - s[,2]/sqrt(s[,3])*norm[,1]+sqrt(s[,1]-s[,2]^2/s[,3])*norm[,2]  
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(betas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(tau2s,c(0.025,0.5,0.975))  
    quants[t,4,] = quantile(sig2s,c(0.025,0.5,0.975))  
    quants[t,5,] = quantile(exp(xs/2),c(0.025,0.5,0.975))  
  }
  return(quants)
}

# Simulated data
set.seed(98765)
n     =  200
alpha =  0.0
beta  =  0.9
tau2  =  0.5
sig2  =  1.0
tau   = sqrt(tau2)
sig   = sqrt(sig2)
x     = rep(alpha/(1-beta),n+1)
true  = c(alpha,beta,tau2,sig2)
names = c("alpha","beta","tau2","sigma2")
for (t in 2:(n+1))
  x[t] = alpha+beta*x[t-1]+tau*rnorm(1)
x = x[2:(n+1)]
y = rnorm(n,x,sig)

par(mfrow=c(1,1))
plot(y,ylim=range(x,y),xlab="Time",ylab="",main="",pch=16)
lines(x,col=2,lwd=2)

# Prior hyperparameters
# ---------------------
m0    = 0.0
C0    = 10
n0    = 10
sig20 = sig2
nu0   = 10
tau20 = tau2
b0    = c(alpha,beta)
B0    = c(1,1)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

# ONE LONG PL FILTER
# ------------------
set.seed(246521)
N      = 100000
xs     = rnorm(N,m0,sC0)
tau2s  = 1/rgamma(N,nu0/2,nu0*tau20/2)
taus   = sqrt(tau2s)
alphas = rnorm(N,b0[1],taus*sB0[1])
betas  = rnorm(N,b0[2],taus*sB0[2])
sig2s  = 1/rgamma(N,n0/2,nu0*sig20/2)
pars   = cbind(alphas,betas,tau2s,sig2s)
par(mfrow=c(2,2))
for (i in 1:4){
  hist(pars[,i],prob=TRUE,xlab="",main=names[i])
  abline(v=true[i],col=2)
}
plm    = PL(y,alphas,betas,tau2s,sig2s,xs)

cols = c(grey(0.5),1,grey(0.5))
ind  = 10:n
n1   = length(ind)
par(mfrow=c(2,2))
for (i in 1:4){
  ts.plot(plm[,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
  abline(h=true[i],lty=2)
}

# Particle filters
# ----------------
set.seed(246521)
N    = 1000
nsim = 100
a    = 0.95
qs   = array(0,c(3,nsim,n,5,3))
for (s in 1:nsim){
  print(s)
  xs         = rnorm(N,m0,sC0)
  sig2s      = 1/rgamma(N,n0/2,n0*sig20/2)
  tau2s      = 1/rgamma(N,nu0/2,nu0*tau20/2)
  taus       = sqrt(tau2s)
  alphas     = rnorm(N,b0[1],taus*sB0[1])
  betas      = rnorm(N,b0[2],taus*sB0[2])
  qs[1,s,,,] = LW(y,alphas,betas,tau2s,sig2s,xs,a)
  qs[2,s,,,] = Storvik(y,alphas,betas,tau2s,sig2s,xs)
  qs[3,s,,,] = PL(y,alphas,betas,tau2s,sig2s,xs) 
}

# Mean square error
quants = c("2.5th","50th","97.5th")
filter = c("LW","Storvik","PL")
mse    = array(0,c(3,n1,4,3))
for (l in 1:3) 
  for (k in 1:4) 
    for (i in 1:nsim) 
      mse[l,,k,]=mse[l,,k,]+(qs[l,i,ind,k,]-plm[ind,k,])^2
sq.mse = sqrt(mse/nsim)  

pdf(file=paste("comparison-",round(100*tau2),".pdf",sep=""),width=10,height=15)
cols = c(grey(0.75),grey(0.5),grey(0.75))
#cols = c(3,2,5)
par(mfrow=c(4,3))
for (k in 1:4){
  L = min(qs[,,ind,k,])
  U = max(qs[,,ind,k,])
  for (l in 1:3){
    plot(qs[l,1,,k,1],ylab=names[k],xlab="Time",main=filter[l],ylim=c(L,U),type="l",col=cols[1])
    for (i in c(1,3,2)) for (j in 1:nsim) lines(qs[l,j,,k,i],col=cols[i])
    for (i in c(1,3,2)) lines(plm[,k,i],lwd=1,col=1)  
  }
}  
dev.off()

pdf(file=paste("boxplot-",round(100*tau2),".pdf",sep=""),width=12,height=10)
cols = c(grey(0.3),gray(0.5),grey(0.7))
U    = c(0.275,0.175,0.45,1.25)
par(mfrow=c(2,2))
for (i in 1:4){
  boxplot(sq.mse[1,,i,1],sq.mse[2,,i,1],sq.mse[3,,i,1],
          sq.mse[1,,i,2],sq.mse[2,,i,2],sq.mse[3,,i,2],
          sq.mse[1,,i,3],sq.mse[2,,i,3],sq.mse[3,,i,3],
          outline=FALSE,ylab="Root MSE",xlab="Percentile",axes=FALSE,col=cols,ylim=c(0,U[i]))
  axis(2);axis(1,at=c(2,5,8),lab=quants);box()
  abline(v=3.5)
  abline(v=6.5)
  title(names[i])
}
legend(0.5,1.2,legend=filters,col=cols,cex=1.25,lwd=c(4,4,4),lty=c(1,1,1),bty="n")
dev.off()


