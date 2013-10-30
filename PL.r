
PL = function(y,alphas,betas,tau2s,xs){
  n      = length(y)
  N      = length(xs)
  quants = array(0,c(n,4,3))
  
  
  
  
  s      = matrix(0,N,7)
  
  
  s[,1] = 1/B0[1]
  s[,2] = 0
  s[,3] = 1/B0[2]
  s[,4] = s[,1]*b0[1]# + s[,2]*b0[2]
  s[,5] = s[,3]*b0[2]# + s[,2]*b0[1]
  s[,6] = c0
  s[,7] = d0
  m      = s[,1]*s[,3]-s[,2]^2
  b1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
  b2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
  
  
  for (t in 1:n){
    # Resampling
    mus    = alphas+betas*xs
    stdevs = exp(mus/2)
    weight = dnorm(y[t],0,stdevs)
    k      = sample(1:N,size=N,replace=T,prob=weight)
    alphas = alphas[k]
    betas  = betas[k]
    tau2s  = tau2s[k]
    
    mus    = mus[k]
    xs1    = xs[k]
    sig2s = exp(xs1/2)
    # Propagating
    vars   = 1/(1/sig2s+1/tau2s)
    sds    = sqrt(vars)
    means  = mus + vars*(y[t]/sig2s)
    xs     = rnorm(N,means,sds)
    so     = s[k,]
    b1o  = b1[k]
    b2o  = b2[k]    
    # Updating sufficient stats
    s[,1]  = so[,1] + 1
    s[,2]  = so[,2] + xs1
    s[,3]  = so[,3] + xs1^2
    s[,4]  = so[,4] + xs
    s[,5]  = so[,5] + xs*xs1
    s[,6]  = so[,6] + 1

    # determinant of Binv matrix
    m      = s[,1]*s[,3]-s[,2]^2 
    
    # b_{t} =  inv(Binv_{t}) x (Binv_{t-1} * b_{t-1}) 
    b1   = (s[,3]*s[,4]-s[,2]*s[,5])/m  
    b2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
    
    s[,7]  = so[,7] + (xs-b1-b2*xs1)*xs + (b1o-b1)*so[,4]+(b2o-b2)*so[,5]
    
    tau2s  = 1/rgamma(N,s[,6]/2,s[,7]/2)
    std    = sqrt(tau2s/m)
    norm   = cbind(rnorm(N,0,std),rnorm(N,0,std))
    
    # cholesky decomposition of B, thus inv(Binv)
    alphas = b1 + sqrt(s[,3])*norm[,1]
    betas  = b2 - s[,2]/sqrt(s[,3])*norm[,1]+sqrt(s[,1]-s[,2]^2/s[,3])*norm[,2]  
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(betas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(tau2s,c(0.025,0.5,0.975))  
    quants[t,4,] = quantile(exp(xs/2),c(0.025,0.5,0.975))  
      
  }
  return(quants)
}

# Simulated data
set.seed(98765)
n     =  1000
alpha =  -.01
beta  =  0.99
tau2  =  0.15^2
sig2  =  1.0
tau   = sqrt(tau2)
x     = rep(alpha/(1-beta),n+1)
true  = c(alpha,beta,tau2,sig2)
names = c("alpha","beta","tau2","y")
for (t in 2:(n+1))
  x[t] = alpha+beta*x[t-1]+tau*rnorm(1)
x = x[2:(n+1)]
y = rnorm(n,0,exp(x/2))

par(mfrow=c(1,1))
plot(y,ylim=range(x,y),xlab="Time",ylab="",main="",pch=16)
lines(x,col=2,lwd=2)

# Prior hyperparameters
# ---------------------
m0    = 0.0


C0 = 10
c0    = 1
d0    = 1
tau20 = tau2
b0    = c(alpha,beta)
B0    = c(1,1)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

# ONE LONG PL FILTER
# ------------------
set.seed(246521)
N      = 10000
xs     = rnorm(N,m0,sC0)
tau2s  = 1/rgamma(N,c0/2,d0/2)
taus   = sqrt(tau2s)
alphas = rnorm(N,b0[1],taus*sB0[1])
betas  = rnorm(N,b0[2],taus*sB0[2])

pars   = cbind(alphas,betas,tau2s)
par(mfrow=c(3,1))
for (i in 1:3){
  hist(pars[,i],prob=TRUE,xlab="",main=names[i])
  abline(v=true[i],col=2)
}
plm    = PL(y,alphas,betas,tau2s,xs)

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
# set.seed(246521)
# N    = 1000
# nsim = 100
# a    = 0.95
# qs   = array(0,c(nsim,n,4,3))
# for (s in 1:nsim){
#   print(s)
#   xs         = rnorm(N,m0,sC0)
#   tau2s      = 1/rgamma(N,c0/2,d0/2)
#   taus       = sqrt(tau2s)
#   alphas     = rnorm(N,b0[1],taus*sB0[1])
#   betas      = rnorm(N,b0[2],taus*sB0[2])
# #   qs[1,s,,,] = LW(y,alphas,betas,tau2s,sig2s,xs,a)
# #   qs[2,s,,,] = Storvik(y,alphas,betas,tau2s,sig2s,xs)
#   qs[s,,,] = PL(y,alphas,betas,tau2s,xs) 
# }
# 
# 
# 
# cols = c(grey(0.5),1,grey(0.5))
# ind  = 10:n
# n1   = length(ind)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   ts.plot(qs[50,,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
#   abline(h=true[i],lty=2)
# }
# 
# 
# # Mean square error
# quants = c("2.5th","50th","97.5th")
# filter = "PL"
# mse    = array(0,c(n1,4,3))
#   for (k in 1:3) 
#     for (i in 1:nsim) 
#       mse[,k,]=mse[,k,]+(qs[i,ind,k,]-plm[ind,k,])^2
# sq.mse = sqrt(mse/nsim)  
# 
# 
# cols = c(grey(0.75),grey(0.5),grey(0.75))
# #cols = c(3,2,5)
# par(mfrow=c(3,1))
# for (k in 1:3){
#   L = min(qs[,ind,k,])
#   U = max(qs[,ind,k,])
# 
#     plot(qs[1,,k,1],ylab=names[k],xlab="Time",main=filter,ylim=c(L,U),type="l",col=cols[1])
#     for (i in c(1,3,2)) for (j in 1:nsim) lines(qs[j,,k,i],col=cols[i])
#     for (i in c(1,3,2)) lines(plm[,k,i],lwd=1,col=1)  
#   
# }  
# 
# 
# 
# cols = c(grey(0.3),gray(0.5),grey(0.7))
# U    = c(0.275,0.175,0.45,1.25)
# par(mfrow=c(2,2))
# for (i in 1:3){
#   boxplot(sq.mse[1,,i,1],sq.mse[2,,i,1],sq.mse[3,,i,1],
#           sq.mse[1,,i,2],sq.mse[2,,i,2],sq.mse[3,,i,2],
#           sq.mse[1,,i,3],sq.mse[2,,i,3],sq.mse[3,,i,3],
#           outline=FALSE,ylab="Root MSE",xlab="Percentile",axes=FALSE,col=cols,ylim=c(0,U[i]))
#   axis(2);axis(1,at=c(2,5,8),lab=quants);box()
#   abline(v=3.5)
#   abline(v=6.5)
#   title(names[i])
# }
# #legend(0.5,1.2,legend=filters,col=cols,cex=1.25,lwd=c(4,4,4),lty=c(1,1,1),bty="n")
# 
# 
