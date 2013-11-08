############################################################################################
#
# STOCHASTIC VOLATILITY AR(1) MODEL
#
# For t=1,2,...,n
#
#    Observation equation: y(t) = exp(h(t)/2)*u(t)     
#    Evolution equation:   h(t) = mu + phi*h(t-1) + tau*e(t) 
# 
# with u(t) and e(t) i.i.d N(0,1). 
#
# Prior specification:
#
#    h(0)    ~ N(m0,C0)
#    tau2    ~ IG(nu0/2,nu0*s02)
#    (mu,ph) ~ N(theta0,tau2*V0)
#
# where m0,C0,nu0,s02,theta0,V0 are known hyperparameters.
#
############################################################################################
#
# HEDIBERT FREITAS LOPES
# Associate Professor of Econometrics and Statistics
# The University of Chicago Booth School of Business
# 5807 South Woodlawn Avenue
# Chicago, Illinois, 60637
# Email : hlopes@ChicagoGSB.edu
# URL: http://faculty.chicagobooth.edu/hedibert.lopes/research/
#
############################################################################################
rm(list=ls())
svm.lw = function(y,alphas,betas,tau2s,xs,a){
  n    = length(y)
  N    = length(alphas)
  h2   = 1-a^2
  pars = cbind(alphas,betas,log(tau2s))
  quants = array(0,c(n,4,3))
  for (t in 1:n){
  	 #print(t)
  	 mpar    = apply(pars,2,mean)
  	 vpar    = var(pars)
  	 ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=TRUE)
  	 mus     = pars[,1]+pars[,2]*xs
  	 w       = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
  	 w1      = exp(w-max(w))
  	 k       = sample(1:N,size=N,replace=TRUE,prob=w1)
  	 ms1     = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
  	 w       = w[k]
  	 xs      = xs[k]
  	 xt      = rnorm(N,ms1[,1]+ms1[,2]*xs,exp(ms1[,3]/2))
  	 w1      = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-w
  	 w1      = exp(w1-max(w1))
  	 k       = sample(1:N,size=N,replace=TRUE,prob=w1)
  	 xs      = xt[k]
  	 pars    = ms1[k,]
  	 quants[t,1,] = quantile(pars[,1],c(0.05,0.5,0.95))
  	 quants[t,2,] = quantile(pars[,2],c(0.05,0.5,0.95))
  	 quants[t,3,] = quantile(exp(pars[,3]),c(0.05,0.5,0.95))
  	 quants[t,4,] = quantile(exp(xs/2),c(0.05,0.5,0.95))
  }
  draws = cbind(pars[,1],pars[,2],exp(pars[,3]))
  return(quants)
}

svm.pl = function(y,alphas,betas,tau2s,xs,b0,B0,c0,d0){
  n      = length(y)
  N      = length(xs)
  nmix   = 7  
  mu     = matrix(c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859),
           N,nmix,byrow=TRUE)
  sig2   = matrix(c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610),
           N,nmix,byrow=TRUE)
  q      = matrix(c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500),
           N,nmix,byrow=TRUE) 
  quants = array(0,c(n,4,3))
  z      = log(y^2)
  s      = matrix(0,N,7)
  s[,1]  = 1.0/B0[1]
  s[,2]  = 0.0
  s[,3]  = 1.0/B0[2]
  s[,4]  = b0[1]/B0[1]
  s[,5]  = b0[2]/B0[2]
  s[,6]  = c0
  s[,7]  = d0
  num1   = rep(b0[1],N)
  num2   = rep(b0[2],N)
  for (t in 1:n){
    #print(t)
    mus    = matrix(alphas+betas*xs,N,nmix)
    probs  = q*dnorm(z[t],mus+mu,sqrt(sig2+tau2s))
    weight = apply(probs,1,sum)
  	  k      = sample(1:N,size=N,replace=TRUE,prob=weight)
  	  xs1    = xs[k]
  	  probs  = probs[k,]
  	  mus    = mus[k,]
  	  alphas = alphas[k]
  	  betas  = betas[k]
  	  tau2s  = tau2s[k]
  	  so     = s[k,]
  	  num1o  = num1[k]
  	  num2o  = num2[k]	  
  	  tau2ss = matrix(tau2s,N,nmix)
  	  vars   = 1/(1/sig2+1/tau2ss)
  	  sds    = sqrt(vars)
  	  means  = vars*((z[t]-mu)/sig2 + mus/tau2ss)
  	  for (i in 1:N){
  	  	  comp  = sample(1:nmix,size=1,prob=probs[i,])
  	    xs[i] = rnorm(1,means[i,comp],sds[i,comp])
    }
  	  s[,1]  = so[,1] + 1
  	  s[,2]  = so[,2] + xs1
  	  s[,3]  = so[,3] + xs1^2
  	  s[,4]  = so[,4] + xs
  	  s[,5]  = so[,5] + xs*xs1
  	  s[,6]  = so[,6] + 1/2
  	  m      = s[,1]*s[,3]-s[,2]^2
  	  num1   = (s[,3]*s[,4]-s[,2]*s[,5])/m
  	  num2   = (s[,1]*s[,5]-s[,2]*s[,4])/m
  	  s[,7]  = so[,7] + (xs-num1-num2*xs1)*xs/2 + ((num1o-num1)*so[,4]+(num2o-num2)*so[,5])/2
  	  tau2s  = 1/rgamma(N,s[,6],s[,7])
  	  std    = sqrt(tau2s/m)
    norm   = cbind(rnorm(N,0,std),rnorm(N,0,std))
    alphas = num1 + sqrt(s[,3])*norm[,1]
    betas  = num2 - s[,2]/sqrt(s[,3])*norm[,1]+sqrt(s[,1]-s[,2]^2/s[,3])*norm[,2]   
  	  quants[t,1,] = quantile(alphas,c(0.05,0.5,0.95))
  	  quants[t,2,] = quantile(betas,c(0.05,0.5,0.95))
  	  quants[t,3,] = quantile(tau2s,c(0.05,0.5,0.95))
  	  quants[t,4,] = quantile(exp(xs/2),c(0.05,0.5,0.95))
  }
  return(quants)
}

############################################################################################
data   = read.table("SA-marketindex.txt",header=TRUE)
n      = nrow(data)
price  = data[,2]
return = data[,3]-mean(data[,3])
ind = seq(1,n,length=9)
dat = data[ind,1]

pdf(file="sv-data.pdf",width=15,height=10)
par(mfrow=c(1,1))
plot(return,ylab="",main="Returns",axes=FALSE,xlab="Days",type="h",lwd=1.5)
axis(2);axis(1,at=ind,lab=dat);box()
for (i in seq(-20,20,by=5))
  abline(h=i,lty=2,col=grey(0.5))
abline(h=0,lwd=3)
dev.off()

# Prior hyperparameters
# ---------------------
m0     = 0
C0     = 1
b0     = c(0,0.95)
B0     = c(0.1,0.05)
nu0    = 4.5
s02    = 0.02777778
c0     = nu0/2
d0     = nu0*s02/2

# Particle filter
set.seed(13121654)
N      = 20000
N1     = 20000
nrep   = 20
quants = array(0,c(4,nrep,n,4,3))
for (rep in 1:nrep){
  print(rep)
  tau2s  = 1/rgamma(N,c0,d0)
  alphas = rnorm(N,b0[1],sqrt(tau2s*B0[1]))
  betas  = rnorm(N,b0[2],sqrt(tau2s*B0[2]))
  priors = cbind(alphas,betas,tau2s)
  xs     = rnorm(N,m0,sqrt(C0))
#   print("lw1")
#   quants[1,rep,,,] = svm.lw(return,alphas,betas,tau2s,xs,0.975)
#   print("lw2")
#   quants[2,rep,,,] = svm.lw(return,alphas,betas,tau2s,xs,0.985)
#   print("lw3")
#   quants[3,rep,,,] = svm.lw(return,alphas,betas,tau2s,xs,0.995)
  print("pl")
  quants[4,rep,,,] = svm.pl(return,alphas[1:N1],betas[1:N1],tau2s[1:N1],xs[1:N1],b0,B0,c0,d0)
}



pdf(file="svm-graph1.pdf",width=15,height=10)
names1 = c("LW-0.975","LW-0.985","LW-0.995","PL")
names = c("alpha","phi","tau2")
par(mfrow=c(3,1))
for (i in 1:3){
  for (j in 4:4){
    plot(quants[j,1,,i,2],ylim=range(quants[,,,i,]),
       main=names1[j],ylab=names[i],axes=FALSE,xlab="Days",type="l")
    axis(2);axis(1,at=ind,lab=dat);box()
    for (rep in 1:nrep)
      for (k in 1:3)
        lines(quants[j,rep,,i,k],col=k)
  }
}
dev.off()

pdf(file="svm-graph2.pdf",width=15,height=10)
names1 = c("LW-0.975","LW-0.985","LW-0.995","PL")
par(mfrow=c(4,1))
i = 4
for (j in 1:4){
  plot(quants[j,1,,i,2],ylim=range(quants[,,,i,]),
       main=names1[j],ylab="Standard deviations",axes=FALSE,xlab="Days",type="l")
  axis(2);axis(1,at=ind,lab=dat);box()
  for (rep in 1:nrep)
    for (k in 1:3)
      lines(quants[j,rep,,i,k],col=k)
}
dev.off()
