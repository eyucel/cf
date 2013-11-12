library(msm)
svm.pl = function(y,alphas,betas,tau2s,ps,qs,xs,zs){#,b0,B0,c0,d0){
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
  s[,12] = 0
  s[,13] = 0
  s[,14] = 0
  s[,15] = 0
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
  
  for (t in 1:n){
    #print(t)
    mus1    = matrix(alphas+gammas+betas*xs,N,nmix)
    mus    = matrix(alphas+betas*xs,N,nmix)
    #     print(mus[1,])
    #     scan(n=1)
#     probs  = q*dnorm(z[t],mus+mu,sqrt(sig2+tau2s))
    switch1 = qs*zs + (1-zs)*(1-ps)
    switch0 = (1-zs)*ps + zs*(1-qs)
    left = q*(dnorm(z[t],mus1+mu,sqrt(sig2+tau2s)))
    right = q*(dnorm(z[t],mus+mu,sqrt(sig2+tau2s)))
    probs = left*switch1+right*switch0
    lefts = apply(left,1,sum)
    rights = apply(right,1,sum)
#     weight = apply(probs,1,sum)
    weight = lefts*switch1+rights*switch0
#     scan(n=1)
    k      = sample(1:N,size=N,replace=TRUE,prob=weight)
    #     print(probs[1,])
    #     scan(n=1)
    
    xs1    = xs[k]
    zs1    = zs[k]
    probs  = probs[k,]
    mus    = mus[k,]
    mus1   = mus1[k,]
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
    
    
    zero.index = zs1 == 0
    one.index = zs1 == 1
    ll = sum(zero.index)
    
    switch1 = qs*zs + (1-zs)*(1-ps)
    switch0 = (1-zs)*ps + zs*(1-qs)
    
    #     zs = zz[t+1,]
    pp = lefts[k]*switch1/weight[k]
    qq = rights[k]*switch0/weight[k]
    
    
#     print(sum(pp<1))
#     print(sum(qq<1))
    
    zs[zero.index]     = rbinom(ll,1,1-pp[zero.index])
    zs[one.index]      = rbinom(N-ll,1,qq[one.index])

#     scan(n=1)
#     print(zs)
    tau2ss = matrix(tau2s,N,nmix)
    vars   = 1/(1/sig2+1/tau2ss)
    sds    = sqrt(vars)
    mus = matrix(alphas+gammas*zs+betas*xs,N,nmix)
    means  = vars*((z[t]-mu)/sig2 + mus/tau2ss)
    for (i in 1:N){
      comp  = sample(1:nmix,size=1,prob=probs[i,])
      xs[i] = rnorm(1,means[i,comp],sds[i,comp])
    }
    

    
    s[,1]  = so[,1] + 1
    s[,2]  = so[,2] + zs1
    s[,3]  = so[,3] + xs1
    s[,4]  = so[,4] + zs1^2
    s[,5]  = so[,5] + zs1*xs1
    s[,6]  = so[,6] + xs1^2
    s[,7]  = so[,7] + xs
    s[,8]  = so[,8] + xs*zs1
    s[,9]  = so[,9] + xs*xs1
    s[,10]  = so[,10] + 1
    
    s[zero.index,12] = so[zero.index,12]+zs1[zero.index]-zs[zero.index]+1
    s[zero.index,13] = so[zero.index,13]+zs[zero.index]
    s[one.index,14] = so[one.index,14]+zs1[one.index]-zs[one.index]
    s[one.index,15] = so[one.index,15]+zs[one.index]-zs1[one.index]+1
    
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
    
    s[,11]  = so[,11] + (xs-b1-b2*zs1-b3*xs1)*xs + (b1o-b1)*so[,7]+(b2o-b2)*so[,8]+(b3o-b3)*so[,9]

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
    
        gamma.mean = b2 - 1/s[,4]*(s[,2]*(alphas-b1)+s[,5]*(betas-b3))
        gammas = rtnorm(N,gamma.mean,sd = sqrt(tau2s/s[,4]),lower=0)
    
    ps = rbeta(N,s[,12]+1,s[,13]+1)
    qs = rbeta(N,s[,15]+1,s[,14]+1)
    
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(gammas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(betas,c(0.025,0.5,0.975))  
    quants[t,4,] = quantile(tau2s,c(0.025,0.5,0.975))  
    quants[t,5,] = quantile(ps,c(0.025,0.5,0.975))  
    quants[t,6,] = quantile(qs,c(0.025,0.5,0.975))  
    quants[t,7,] = quantile(zs,c(0.025,0.5,0.975))
    quants[t,8,] = quantile(exp(xs/2),c(0.025,0.5,0.975))
    zmean[t] =  mean(zs)
  }
  return(list(quants=quants,zmean=zmean))
}
PL = function(y,alphas,betas,tau2s,ps,qs,xs,zs){
  n      = length(y)
  N      = length(xs)
  quants = array(0,c(n,7,3))
  
  
  
#   zz = matrix(zs,n,N)
  s  = matrix(0,N,15)
  # matrix Binv looks like
  # | s1   s2    s3 |
  # | s2   s4    s5 |
  # | s3   s5    s6 |
  
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
  

  for (t in 1:n){
#     zs = zz[t,]
    # Resampling
    mus    = alphas+gammas+betas*xs
    stdevs = exp(mus/2)
    stdevs1 = exp((alphas+betas*xs)/2)
    weight = dnorm(y[t],0,stdevs)*(qs+(1-zs)*(1-ps)) + dnorm(y[t],0,stdevs1)*((1-zs)*ps+zs*(1-qs))
#     weight = dnorm(y[t],0,stdevs)
    k      = sample(1:N,size=N,replace=T,prob=weight)
    alphas = alphas[k]
    betas  = betas[k]
    tau2s  = tau2s[k]
    gammas = gammas[k]
    ps = ps[k]
    qs = qs[k]
    
    mus    = mus[k]
    xs1    = xs[k]
    zs1    = zs[k]
    sig2s = exp(xs1)
    # Propagating
    vars   = 1/(1/sig2s+1/tau2s)
    sds    = sqrt(vars)
    means  =  vars*(mus/tau2s +y[t]/sig2s)
    zero.index = zs1 == 0
    one.index = zs1 == 1
    ll = sum(zero.index)
    #     zs = zz[t+1,]
    zs[zero.index]     = rbinom(ll,1,1-ps[zero.index])
    zs[one.index]      = rbinom(N-ll,1,qs[one.index])
#     mus    = alphas+gammas*zs+betas*xs1
    xs     = rnorm(N,means,sds)
    xs = alphas+betas*xs + gammas*zs + rnorm(N,0,sqrt(tau2s))

#     print(zs)
#     print(qs)
#     print(ps)
    so     = s[k,]
#     print(s[,12:15])
#     print(so[,12:15])
#     print(k)
#     scan(n=1)
    b1o  = b1[k]
    b2o  = b2[k] 
    b3o  = b3[k]
    # Updating sufficient stats
    s[,1]  = so[,1] + 1
    s[,2]  = so[,2] + zs1
    s[,3]  = so[,3] + xs1
    s[,4]  = so[,4] + zs1^2
    s[,5]  = so[,5] + zs1*xs1
    s[,6]  = so[,6] + xs1^2
    s[,7]  = so[,7] + xs
    s[,8]  = so[,8] + xs*zs1
    s[,9]  = so[,9] + xs*xs1
    s[,10]  = so[,10] + 1
    
    s[zero.index,12] = so[zero.index,12]+zs1[zero.index]-zs[zero.index]+1
    s[zero.index,13] = so[zero.index,13]+zs[zero.index]
    s[one.index,14] = so[one.index,14]+zs1[one.index]-zs[one.index]
    s[one.index,15] = so[one.index,15]+zs[one.index]-zs1[one.index]+1
#     
#     print(s[,12:15])
#     print(zs1)
#     print(zs)
#     print(zs1[zero.index]-zs[zero.index]+1)
#     print(zs[zero.index])
#     print(zs1[one.index]-zs[one.index])
#     print(zs[one.index]-zs1[one.index]+1)
#     scan(n=1)
    blk.A = s[,4]*s[,6] - s[,5]^2
    blk.B = -(s[,2]*s[,6]-s[,3]*s[,5])
    blk.C = s[,2]*s[,5]-s[,3]*s[,4]
    #   blk.D = blk.B
    blk.E = s[,1]*s[,6]-s[,3]^2
    blk.F = -(s[,1]*s[,5]-s[,2]*s[,3])
    #   blk.G = blk.C
    #   blk.H = blk.F
    blk.I  = s[,1]*s[,4]-s[,2]^2
    
    
#     ti = scan(n=1)
#     print(s[ti,])
    
    
    # determinant of Binv matrix
    m    = s[,1]*blk.A + s[,2] * blk.B + s[,3] * blk.C
    
#     AA = matrix(c(s[ti,1],s[ti,2],s[ti,3],s[ti,2],s[ti,4],s[ti,5],s[ti,3],s[ti,5],s[ti,6]),3,3)
#     
#     print(solve(AA))
#     print(AA)
#     print(matrix(c(blk.A[ti],blk.B[ti],blk.C[ti], blk.B[ti],blk.E[ti], blk.F[ti], blk.C[ti], blk.F[ti], blk.I[ti]),3,3)/m[ti])
#      print(c(blk.A[1]/m[1],blk.B[1]/m[1],blk.C[1]/m[1]))
#      print(solve(AA))
    # b_{t} =  inv(Binv_{t}) x (Binv_{t-1} * b_{t-1})     
    b1   = (s[,7]*blk.A + s[,8] * blk.B + s[,9] * blk.C)/m
    b2   = (s[,7]*blk.B + s[,8] * blk.E + s[,9] * blk.F)/m
    b3   = (s[,7]*blk.C + s[,8] * blk.F + s[,9] * blk.I)/m
#

#     print('b2')
#     print((s[,7]*blk.B)[ti])
# print((blk.F)[ti])
# print((s[,9])[ti])
#     print(m[ti])
#     print(b2[ti])
#     F.new = rbind(1,zs1[ti],xs1[ti])
#     C.old = diag(B0,3)
#     Q = 1
# #     print(zs1[ti])
#     D1 = t(F.new) %*% C.old%*% F.new + Q 
#       
#       
#       #       C0 = diag(B0,2)
#       #       print( C0 %*%F.new[,j]  %*% solve(D1) %*% (xs[1]-t(F.new[,j])%*%b0))
#       CF = C.old %*% F.new
# #       
#       D.new = t(F.new)  %*% CF +Q
#       D.inv = solve(D.new)
# #       print(D.inv)
#     
#       m.old = b0
# #       print(D.new)
# #       print(CF %*% D.inv %*% (xs[j]-t(F.new[,j]) %*% m.old[j,]))
# #       m.old = m.old + CF %*% D.inv %*% (xs[1]-t(F.new) %*% m.old)
# #     print(D.inv %*% (xs[3]-t(F.new) %*% m.old))
#     
#       C.old = C.old - CF %*% D.inv  %*% t(F.new) %*% C.old
#       d0 = d0 + t(xs[ti]-t(F.new)%*% m.old)%*% D.inv %*% (xs[ti]-t(F.new) %*% m.old)
#     m.old = m.old + CF %*% D.inv %*% (xs[ti]-t(F.new) %*% m.old)
#     print(d0)
#     print(C.old)
#     print(chol(C.old))
#     print(m.old)
#     print(c(b1[ti],b2[ti],b3[ti]))
#     print(F.new)
#     print(m.old)
#     print(c(b1o[3],b2o[3],b3o[3]))
#     print(c(zs1[3],xs1[3]))
#     print(((xs-b1-b2*zs1-b3*xs1)*xs+(b1o-b1)*so[,7]+(b2o-b2)*so[,8]+(b3o-b3)*so[,9])[3])
#       print(m.old)
#     print(c(b1[1],b2[1],b3[1]))
#     scan(n=1)
#     
    #print(xs)
    s[,11]  = so[,11] + (xs-b1-b2*zs1-b3*xs1)*xs + (b1o-b1)*so[,7]+(b2o-b2)*so[,8]+(b3o-b3)*so[,9]
#     print(s[ti,11])
    tau2s  = 1/rgamma(N,s[,10]/2,s[,11]/2)
    std    = sqrt(tau2s/m)
    
    norm   = cbind(rnorm(N,0,std),rnorm(N,0,std),rnorm(N,0,std))
    
    # cholesky decomposition of B, thus inv(Binv)
    
    L11 = sqrt(blk.A)
#     print(L11[ti])
    L21 = blk.B/L11
#     print(L21[ti])
    L31 = blk.C/L11
#     print(L31[ti])
    
    L22 = sqrt(blk.E - L21^2)
#     print(L22[ti])
    L32 = (blk.F-L31*L21)/L22
#     print(L32[ti])
    L33 = sqrt(blk.I - L31^2 - L32^2)
#     print(L33[ti]))
    
    alphas = b1 + L11*norm[,1]
    gammas = b2 + L21*norm[,1]+L22*norm[,2]
#     gamma.mean = b2[b2>0]
    
#     gammas = rtnorm(N,mean=b2,lower=0)
    betas  = b3 + L31*norm[,1] + L32 * norm[,2] + L33 * norm[,3]
    
#     gamma.mean = b2 - 1/s[,4]*(s[,2]*(alphas-b1)+s[,5]*(betas-b3))
#     gammas = rtnorm(N,gamma.mean,sd = sqrt(tau2s/s[,4]),lower=0)
    
    ps = rbeta(N,s[,12]+1,s[,13]+1)
    qs = rbeta(N,s[,15]+1,s[,14]+1)
    # Storing quantiles
    quants[t,1,] = quantile(alphas,c(0.025,0.5,0.975))
    quants[t,2,] = quantile(gammas,c(0.025,0.5,0.975))
    quants[t,3,] = quantile(betas,c(0.025,0.5,0.975))  
    quants[t,4,] = quantile(tau2s,c(0.025,0.5,0.975))  
    quants[t,5,] = quantile(ps,c(0.025,0.5,0.975))  
    quants[t,6,] = quantile(qs,c(0.025,0.5,0.975))  
    quants[t,7,] = quantile(zs,c(0.025,0.5,0.975))

  }
  
  return(quants)
}

# Simulated data
# set.seed(98765)
n     =  1000
alpha =  -1
gamma = 1.5
beta  =  0.7
tau2  =  .4
sig2  =  1.0
tau   = sqrt(tau2)
p0 = .99
q0 = .985

x     = rep(alpha/(1-beta),n+1)
S = rep(0,n+1)
S[1] = rbinom(1,1,p0/(1-q0)/100)
true  = c(alpha,gamma,beta,tau2,p0,q0)
names = c("alpha","gamma","beta","tau2","p","q")


for (t in 2:(n+1))
  {
  if (S[t-1] == 0)
    S[t] = rbinom(1,1,1-p0)
  if (S[t-1] == 1)
    S[t] = rbinom(1,1,q0)

  x[t] = alpha+gamma*S[t]+beta*x[t-1]+tau*rnorm(1)}
x = x[2:(n+1)]
y = rnorm(n,0,exp(x/2))
# print(S)
par(mfrow=c(1,1))
plot(y,ylim=range(x,y),xlab="Time",ylab="",main="",pch=16)
lines(x,col=2,lwd=2)

# data = read.csv('OILPRICE.csv')
# prices = data[,2]
# n = length(prices)
# returns = diff(prices)/prices[1:(n-1)]
# logreturns = log(prices[2:n]/prices[1:(n-1)])[300:(n-1)]
# n = length(logreturns)
# y = returns[300:n]+rnorm(n-299,0,.00001)
# Prior hyperparameters
# ---------------------
m0    = 0.0


C0 = 10
c0    = 1
d0    = 1
tau20 = tau2
b0    = c(alpha,gamma,beta)
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
ps  = rbeta(N,1,1)
qs = rbeta(N,1,1)
pars   = cbind(alphas,gammas,betas,tau2s,ps,qs)
par(mfrow=c(2,2))
for (i in 1:4){
  hist(pars[,i],prob=TRUE,xlab="",main=names[i])
  abline(v=true[i],col=2)
}
print(date())
# plm    = PL(y,alphas,betas,tau2s,ps,qs,xs,zs)
out   = svm.pl(y,alphas,betas,tau2s,ps,qs,xs,zs)
plm = out$quants

print(date())
cols = c(grey(0.5),1,grey(0.5))
ind  = 20:length(y)
n1   = length(ind)
par(mfrow=c(3,2))
for (i in 1:6){
  ts.plot(plm[ind,i,],xlab="",ylab="",main=names[i],col=cols,ylim=range(plm[ind,i,]))
  abline(h=true[i],lty=2)
}


# par(mfrow=c(1,1))
# ts.plot(cbind(2*plm[,8,2],-2*plm[,8,2]))
# lines(returns*100,col='red')

# ts.plot(out$zmean)
# par(mfrow=c(1,1))
# ts.plot(plm[,7,2],xlab="",ylab="",main="states",col=1,ylim=c(-1,2),lwd=2)
# points(S[2:(n+1)],col='red',ylim=c(-1,2),lwd=1)


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





################ ETC
# 
# for (j in 1:N){
#   D1 = F.new[,j] %*% C.old[j,,]%*% F.new[,j] + Q 
#   
#   print(D1)
#   #       C0 = diag(B0,2)
#   #       print( C0 %*%F.new[,j]  %*% solve(D1) %*% (xs[1]-t(F.new[,j])%*%b0))
#   CF = C.old[j,,] %*% F.new[,j]
#   
#   D.new = t(F.new[,j])  %*% CF +Q
#   D.inv = solve(D.new)
#   print(D.new)
#   print(CF %*% D.inv %*% (xs[j]-t(F.new[,j]) %*% m.old[j,]))
#   m.old[j,] = m.old[j,] + CF %*% D.inv %*% (xs[j]-t(F.new[,j]) %*% m.old[j,])
#   C.old[j,,] = C.old[j,,] - CF %*% D.inv  %*% t(F.new[,j]) %*% C.old[j,,]
#   print(m.old[1,])
#   print(b1[1])
#   print(b2[1])
#   scan(n=1)
# }