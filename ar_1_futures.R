data = read.csv('WTI3.csv')
y = data[,2:13] 
dates = as.Date(data[,1], format = '%m/%d/%y')
data.cleaned = data[which(dates>=as.Date("1/31/2007",format = '%m/%d/%Y')),]

prompt_month = rev(data.cleaned$P1)


prompt_month.lag = diff(prompt_month)
n = length(prompt_month.lag)
pct_chg = prompt_month.lag/prompt_month[1:n]

p = 2
X = cbind(rep(1,n-1),pct_chg[1:(n-1)])
Y = pct_chg[2:n]

tX = t(X)%*%X
tXY = t(X)%*%Y

# Priors
#####################
m0 = rep(0,p)
C0 = diag(100,p)
iC0 = solve(C0)


nu0 = 3
nusig0 = nu0*2

######################


M = 5000
BETAs = matrix(0,M,p)
SIGs  = rep(0,M)

## Initial value
#################

sigma2 = 1/rgamma(1,nu0/2,nusig0/2)

for(i in 1:M){
  
  
  ## Posterior Full Conditional for Beta
  ######################################
  
  C1 = solve(tX/sigma2 + iC0)
  m1 = C1%*%(tXY/sigma2 + iC0%*%m0)
  
  beta = m1 + t(chol(C1))%*%rnorm(p)
  
  
  ## Posterior Full Conditional for sigma
  #######################################
  
  e = Y - X%*%(beta)
  
  
  nu1 = n + nu0
  nusig1 = nusig0 + sum(e^2)
  
  sigma2 = 1/rgamma(1,nu1/2,nusig1/2)
  
  
  ## Save
  BETAs[i,] = beta
  SIGs[i] = sigma2
  
  print(i)
  
  
}

beta_hat = apply(BETAs[200:M,],2,mean)
sigma2_hat = mean(SIGs[200:M])


# k = 10
# r_k  = rep(0,M)
# for(i in 1:M){
#     sum_beta = 0
#     sum_beta_e = 0
#     for (j in 0:k){
#         sum_beta = sum_beta + beta_hat[2]^j
#         sum_beta_e = sum_beta_e + beta_hat[2]^(k-j)*rnorm(1,0,sqrt(sigma2_hat))
#     }
#     r_k[i] = beta_hat[1] * sum_beta + beta_hat[2]^(k+1)*r[n]+sum_beta_e
# }


par(mfrow=c(3,3))
for(i in 1:p) {
  hist(BETAs[,i],main=paste("Beta",i))
  #     abline(v=MLE[i],col=2,lwd=2)           
}  

par(mfrow=c(1,1))
hist(sqrt(SIGs),main="sigma")
# abline(v=sqrt(sig2MLE),col=2,lwd=2)	



par(mfrow=c(2,2))
for(i in 1:p) {
  plot(BETAs[,i],type="l")
  #     abline(h=MLE[i],col=2,lwd=2)
}

plot(SIGs,type="l")
# abline(h=sig2MLE,col=2,lwd=2)
par(mfrow=c(1,1))