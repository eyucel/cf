svm.pl = function(y,alphas,betas,tau2s,xs){#,b0,B0,c0,d0){
  n      = length(y)
  N      = length(xs)
  nmix   = 7  
  mu     = matrix(c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859),
                  N,nmix,byrow=TRUE)
  sig2   = matrix(c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610),
                  N,nmix,byrow=TRUE)
  q      = matrix(c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500),
                  N,nmix,byrow=TRUE) 
  #   q = c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500)
  #   mu = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
  #   sig2 = c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610)
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
    #     print(mus[1,])
    #     scan(n=1)
    probs  = q*dnorm(z[t],mus+mu,sqrt(sig2+tau2s))
    weight = apply(probs,1,sum)
    k      = sample(1:N,size=N,replace=TRUE,prob=weight)
    #     print(probs[1,])
    #     scan(n=1)
    
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
  hists = list(hist(alphas),hist(betas),hist(tau2s))
  return(list(quants=quants,hists=hists))
}




par(mfrow=c(1,1))
usePrices = TRUE
estlen = 6


  data = read.csv('weeklyoilprice.csv')
  pricedata = data[(1:length(data[,2])),2]
  
  
  m = length(pricedata)
  #y is return
  r = (pricedata[2:m]/pricedata[1:(m-1)]-1)*100
  n = length(r)-estlen+1
  #h is volatility
  h = rep(0,n-1)
  y=r[estlen:(m-1)]
  for(i in 1:(n-1)){
    if(y[i]==0){y[i]=0.00001}
    h[i] = sd(r[i:(i+estlen)])
  }
  plot(pricedata,type="l")


plot(h,type="l")
plot(y,type="l")

# Prior hyperparameters
# ---------------------
m0    = 0.0


C0 = 20
c0    = 1
d0    = 1
tau20 = C0
b0    = c(m0,m0)
B0    = c(1,1)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

# ONE LONG PL FILTER
# ------------------
# set.seed(246521)
N      = 10000
xs     = rnorm(N,m0,sC0)
tau2s  = 1/rgamma(N,c0/2,d0/2)
taus   = sqrt(tau2s)
alphas = rnorm(N,b0[1],taus*sB0[1])
betas  = rnorm(N,b0[2],taus*sB0[2])

pars   = cbind(alphas,betas,tau2s)

out    = svm.pl(y-mean(y),alphas,betas,tau2s,xs)
plm = out$quants
names = c("mu","phi","sig2")
cols = c(grey(0.5),1,grey(0.5))
ind  = 10:n
n1   = length(ind)
par(mfrow=c(3,2))
i = 1
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])
i = 2
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])
i = 3
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])



# for (i in 1:3){
#   #ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
#   ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=c(0,1))
#   hist(plm[ind,i,2])
#   #abline(h=true[i],lty=2)
# }

mh  =plm[,4,2]
h05 = plm[,4,1]
h95 = plm[,4,3]
mhm   = mh
h05m  = h05
h95m  = h95
parsm = pars

par(mfrow=c(1,1))
plot(1+h,type="l",ylim=c(0.0,max(1+h95m)),xlab="time",ylab="")
lines(1+h,col=1,lwd=1)
lines(1+mhm,col=2,lwd=1)
lines(1+h05m,col=3,lwd=1)
lines(1+h95m,col=3,lwd=1)
legend('topright',legend=c("h(t)","E(h(t)|data)","90%CI"),col=c(1,2,3),lty=rep(1,3),lwd=rep(1,3),bty="n")
#lines(abs(y)/max(abs(y)),type="h")

sqrt(sum((mh[1:1445]-h[1:1445])^2))