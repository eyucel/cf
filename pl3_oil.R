
svm.pl3 = function(y,alphas,betas,tau2s,ps,qs,xs,zs){#,b0,B0,c0,d0){
  n      = length(y)
  
  N      = length(xs)
  nmix   = 7  
  mu     = matrix(c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859),
                  N,nmix,byrow=TRUE)
  sig2   = matrix(c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610),
                  N,nmix,byrow=TRUE)
  q      = matrix(c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500),
                  N,nmix,byrow=TRUE) 
  
  quants = array(0,c(n,8,3))
  z      = log(y^2)
  #   z = y
  s      = matrix(0,N,15)
  zmean=rep(0,n)
  
  s[,1] = 1/B0[1]
  s[,2] = 0
  s[,3] = 0
  s[,4] = 1/B0[2]
  s[,5] = 0
  s[,6] = 1/B0[2]
  s[,7] = s[,1]*b0[1]+s[,2]*b0[2]+s[,3]*b0[3]
  s[,8] = s[,2]*b0[1]+s[,4]*b0[2]+s[,5]*b0[3]
  s[,9] = s[,3]*b0[1]+s[,5]*b0[2]+s[,6]*b0[3]
  s[,10] = c0
  s[,11] = d0
  s[,12] = 1
  s[,13] = 1
  s[,14] = 1
  s[,15] = 1
  blk.A = s[,4]*s[,6] - s[,5]^2
  blk.B = -(s[,2]*s[,6]-s[,3]*s[,5])
  blk.C = s[,2]*s[,5]-s[,3]*s[,5]
  #   blk.D = blk.B
  blk.E = s[,1]*s[,6]-s[,3]^2
  blk.F = -(s[,3]*s[,5]-s[,2]*s[,3])
  #   blk.G = blk.C
  #   blk.H = blk.F
  blk.I  = s[,1]*s[,4]-s[,2]^2
  
  m    = s[,1]*blk.A + s[,2] * blk.B + s[,3] * blk.C
  b1   = (s[,7]*blk.A + s[,8] * blk.B + s[,9] * blk.C)/m
  b2   = (s[,7]*blk.B + s[,8] * blk.E + s[,9] * blk.F)/m
  b3   = (s[,7]*blk.C + s[,8] * blk.F + s[,9] * blk.I)/m
  zs1 = zs

    for (t in 1:n){
      #     print(t)
      #print(t)
      #     mus1    = matrix(alphas+gammas+betas*xs,N,nmix)
      mus    = matrix(alphas+betas*xs,N,nmix)
      #     print(mus[1,])
      #     scan(n=1)
      #     probs  = q*dnorm(z[t],mus+mu,sqrt(sig2+tau2s))
      
      zero.index= zs == 0
      one.index = zs == 1
      #     
      #     
      switch1 = qs*one.index + zero.index*(1-ps)
      switch0 = 1-switch1 #(1-zs)*ps + zs*(1-qs)
      left = q*(dnorm(z[t],mus+gammas+mu,sqrt(sig2+tau2s)))
      right = q*(dnorm(z[t],mus+mu,sqrt(sig2+tau2s)))
      probs = left*switch1+right*switch0
      #     lefts = apply(left,1,sum)
      #     rights = apply(right,1,sum)
      #     weight = lefts*switch1+rights*switch0
      
      #     probs = q*(dnorm(z[t],mus1+mu,sqrt(sig2+tau2s))*(qs*zs+(1-zs)*(1-ps)) + dnorm(z[t],mus+mu,sqrt(sig2+tau2s))*((1-zs)*ps+zs*(1-qs)))
      
      weight = apply(probs,1,sum)
      k      = sample(1:N,size=N,replace=TRUE,prob=weight)
      #     print(probs[1,])
      #     scan(n=1)
      
      xs1    = xs[k]
      zs1    = zs[k]
      probs  = probs[k,]
      mus    = mus[k,]
      alphas = alphas[k]
      gammas = gammas[k]
      betas  = betas[k]
      tau2s  = tau2s[k]
      ps = ps[k]
      qs = qs[k]
      so     = s[k,]
      b1o  = b1[k]
      b2o  = b2[k]    
      b3o  = b3[k]
      
      
      
      
      
      tau2ss = matrix(tau2s,N,nmix)
      vars   = 1/(1/sig2+1/tau2ss)
      sds    = sqrt(vars)
      
      zero.index = zs1 == 0
      one.index  = zs1 == 1
      
      #     mus    = matrix(alphas+gammas+betas*xs1,N,nmix)
      #     mus1    = matrix(alphas+betas*xs1,N,nmix)
      
      
      #     switch1 = qs*one.index + zero.index*(1-ps)
      #     switch0 = 1-switch1 #(1-zs)*ps + zs*(1-qs)
      #     left = q*(dnorm(z[t],mus1+mu,sqrt(sig2+tau2s)))
      #     right = q*(dnorm(z[t],mus+mu,sqrt(sig2+tau2s)))
      #     probs = left*switch1+right*switch0
      #     lefts = apply(left,1,sum)
      #     rights = apply(right,1,sum)
      #     weight = lefts*switch1+rights*switch0
      #     
      #     pp = lefts*switch1/weight
      #     qq = rights*switch0/weight
      
      
      zs     = rbinom(N,1,1-ps)*zero.index
      zs    = zs+ rbinom(N,1,qs)*one.index
      
      mus = matrix(alphas+gammas*zs+betas*xs1,N,nmix)
      #     mus = mus+gammas*zs
      probs = q*(dnorm(z[t],mus+mu,sqrt(sig2+tau2s)))
      means  = vars*((z[t]-mu)/sig2 + mus/tau2ss)
      for (i in 1:N) { 
        comp  = sample(1:nmix,size=1,prob=probs[i,])
        xs[i] = rnorm(1,means[i,comp],sds[i,comp])
      }
      
      zero.index = zs1 == 0
      one.index = zs1 == 1
      zero.index.2 = zs == 0
      one.index.2 = zs == 1
      
      #     zs = zz[t+1,]
      
      s[,1]  = so[,1] + 1
      s[,2]  = so[,2] + zs
      s[,3]  = so[,3] + xs1
      s[,4]  = so[,4] + zs^2
      s[,5]  = so[,5] + zs*xs1
      s[,6]  = so[,6] + xs1^2
      s[,7]  = so[,7] + xs
      s[,8]  = so[,8] + xs*zs
      s[,9]  = so[,9] + xs*xs1
      s[,10]  = so[,10] + 1
      
      #     s[zero.index,12] = so[zero.index,12]+zs1[zero.index]-zs[zero.index]+1
      #     s[zero.index,13] = so[zero.index,13]+zs[zero.index]
      #     s[one.index,14] = so[one.index,14]+zs1[one.index]-zs[one.index]
      #     s[one.index,15] = so[one.index,15]+zs[one.index]-zs1[one.index]+1
      s[,12] = so[,12]+zero.index*zero.index.2
      s[,13] = so[,13]+zero.index*one.index.2
      s[,14] = so[,14]+one.index*zero.index.2
      s[,15] = so[,15]+one.index*one.index.2
      blk.A = s[,4]*s[,6] - s[,5]^2
      blk.B = -(s[,2]*s[,6]-s[,3]*s[,5])
      blk.C = s[,2]*s[,5]-s[,3]*s[,4]
      #   blk.D = blk.B
      blk.E = s[,1]*s[,6]-s[,3]^2
      blk.F = -(s[,1]*s[,5]-s[,2]*s[,3])
      #   blk.G = blk.C
      #   blk.H = blk.F
      blk.I  = s[,1]*s[,4]-s[,2]^2
      
      
      # determinant of Binv matrix
      m    = s[,1]*blk.A + s[,2] * blk.B + s[,3] * blk.C
      
      # b_{t} =  inv(Binv_{t}) x (Binv_{t-1} * b_{t-1})     
      b1   = (s[,7]*blk.A + s[,8] * blk.B + s[,9] * blk.C)/m
      b2   = (s[,7]*blk.B + s[,8] * blk.E + s[,9] * blk.F)/m
      b3   = (s[,7]*blk.C + s[,8] * blk.F + s[,9] * blk.I)/m
      
      s[,11]  = so[,11] + (xs-b1-b2*zs-b3*xs1)*xs + (b1o-b1)*so[,7]+(b2o-b2)*so[,8]+(b3o-b3)*so[,9]
      
      tau2s  = 1/rgamma(N,s[,10]/2,s[,11]/2)
      std    = sqrt(tau2s/m)
      
      norm   = cbind(rnorm(N,0,std),rnorm(N,0,std),rnorm(N,0,std))
      
      # cholesky decomposition of B, thus inv(Binv)    
      L11 = sqrt(blk.A)
      L21 = blk.B/L11
      L31 = blk.C/L11
      
      L22 = sqrt(blk.E - L21^2)
      L32 = (blk.F-L31*L21)/L22
      
      L33 = sqrt(blk.I - L31^2 - L32^2)
      
      alphas = b1 + L11*norm[,1]
      gammas = b2 + L21*norm[,1]+L22*norm[,2]
      #     gamma.mean = b2[b2>0]
      
      #     gammas = rtnorm(N,mean=b2,lower=0)
      betas  = b3 + L31*norm[,1] + L32 * norm[,2] + L33 * norm[,3]
      
      #     gamma.mean = b2 - 1/s[,4]*(s[,2]*(alphas-b1)+s[,5]*(betas-b3))
      #     gammas = rtnorm(N,gamma.mean,sd = sqrt(tau2s/s[,4]),lower=0)
      
      ps = rbeta(N,s[,12]+1,s[,13]+1)
      qs = rbeta(N,s[,15]+1,s[,14]+1)
    #     print(ps)
    #     print(qs)
    #     scan(n=1)
    # Storing quantiles
    quantvec = c(0.05,0.5,0.95)
    quants[t,1,] = quantile(alphas,quantvec)
    quants[t,2,] = quantile(gammas,quantvec)
    quants[t,3,] = quantile(betas,quantvec)  
    quants[t,4,] = quantile(tau2s,quantvec)  
    quants[t,5,] = quantile(ps,quantvec)  
    quants[t,6,] = quantile(qs,quantvec)  
    quants[t,7,] = quantile(zs1,quantvec)
    quants[t,8,] = quantile(exp(xs/2),quantvec)
    zmean[t] =  mean(zs1)
  }
  hists = list(hist(alphas),hist(gammas),hist(betas),hist(tau2s),hist(ps),hist(qs))
  return(list(quants=quants,zmean=zmean,hists=hists))
}






usePrices = TRUE
estlen = 6


data = read.csv('weeklyoilprice.csv')
pricedata = data[(1:length(data[,2])),2]

par(mfrow=c(1,1))
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

pa = 1
pb = 1
qa = 1
qb = 1
m0    = 0.0
C0 = 20
tau2=C0
c0    = 1
d0    = 1
tau20 = tau2
b0    = c(m0,m0,m0)
B0    = c(1,1,1)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

# ONE LONG PL FILTER
# ------------------
# set.seed(246521)
N      = 10000
xs     = rnorm(N,m0,sC0)
zs = rbinom(N,1,.5)
tau2s  = 1/rgamma(N,c0/2,d0/2)
taus   = sqrt(tau2s)
alphas = rnorm(N,b0[1],taus*sB0[1])
gammas = rtnorm(N,b0[2],taus*sB0[2],lower=0)
betas  = rnorm(N,b0[3],taus*sB0[3])
ps  = rbeta(N,pa,pb)
qs = rbeta(N,qa,qb)
pars   = cbind(alphas,gammas,betas,tau2s,ps,qs)

# for (i in 1:4){
#   hist(pars[,i],prob=TRUE,xlab="",main=names[i])
#   abline(v=true[i],col=2)
# }
print(date())
# plm    = PL(y,alphas,betas,tau2s,ps,qs,xs,zs)
out   = svm.pl3(y,alphas,betas,tau2s,ps,qs,xs,zs)
plm = out$quants
cols = c(grey(0.5),1,grey(0.5))
ind  = 10:n
n1   = length(ind)
par(mfrow=c(3,2))
names = c("alpha0","alpha1","phi","sigma2","p","q")
i = 1
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])
i = 2
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])
i = 3
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])
i = 4
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])

i = 5
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])

i = 6
ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
plot(out$hists[[i]],main=names[i])




# for (i in 1:3){
#   #ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
#   ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=c(0,1))
#   hist(plm[ind,i,2])
#   #abline(h=true[i],lty=2)
# }

mh  =plm[,8,2]
h05 = plm[,8,1]
h95 = plm[,8,3]
mhm   = mh
h05m  = h05
h95m  = h95
parsm = pars

par(mfrow=c(1,1))
plot(1+h,type="l",ylim=c(0.0,max(1+h95m)),xlab="time",ylab="")

lines(1+mhm,col=2,lwd=1,type='l')
lines(1+h05m,col=3,lwd=1)
lines(1+h95m,col=3,lwd=1)
legend('topright',legend=c("h(t)","E(h(t)|data)","90%CI"),col=c(1,2,3),lty=rep(1,3),lwd=rep(1,3),bty="n")
#lines(abs(y)/max(abs(y)),type="h")


sqrt(sum((mh[1:1445]-h[1:1445])^2))