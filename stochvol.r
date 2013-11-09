#-------------------------------------------------------
# Univariate FFBS: 
# y(t)     ~ N(alpha(t)+F(t)*theta(t);V(t))        
# theta(t) ~ N(gamma+G*theta(t-1);W)               
#-------------------------------------------------------
ffbsu = function(y,F,alpha,V,G,gamma,W,a1,R1,nd=1){
  n = length(y)
  if (length(F)==1){F = rep(F,n)}
  if (length(alpha)==1){alpha = rep(alpha,n)}
  if (length(V)==1){V = rep(V,n)}
  a = rep(0,n)
  R = rep(0,n)
  m = rep(0,n)
  C = rep(0,n)
  B = rep(0,n-1)
  H = rep(0,n-1)
  # time t=1
  a[1] = a1
  R[1] = R1
#   f    = alpha[1]+F[1]*a[1]
#   Q    = R[1]*F[1]**2+V[1]
#   A    = R[1]*F[1]/Q
#   m[1] = a[1]+A*(y[1]-f)
#   C[1] = R[1]-Q*A**2
#   # forward filtering
#   for (t in 2:n){
#     a[t] = gamma + G*m[t-1]
#     R[t] = C[t-1]*G**2 + W
#     f    = alpha[t]+F[t]*a[t]
#     Q    = R[t]*F[t]**2+V[t]
#     A    = R[t]*F[t]/Q
#     m[t] = a[t]+A*(y[t]-f)
#     C[t] = R[t]-Q*A**2
#     B[t-1] = C[t-1]*G/R[t]
#     H[t-1] = sqrt(C[t-1]-R[t]*B[t-1]**2)
#   }
  # backward sampling
  theta = matrix(0,nd,n)
  
  theta = ffbs.cpp(y,F,a,R,m,C,B,H,alpha,V,G,gamma,W,theta)
#   theta[,n] = rnorm(nd,m[n],sqrt(C[n]))
#   for (t in (n-1):1)
#     theta[,t] = rnorm(nd,m[t]+B[t]*(theta[,t+1]-a[t+1]),H[t])
  if (nd==1){
    theta[1,]
  }
  else{
    theta
  }
}

ffbs.cpp = cxxfunction(
  signature(y = "numeric",
            F = "numeric",
            a = "numeric",
            R = "numeric",
            m = "numeric",
            C = "numeric",
            B = "numeric",
            H = "numeric",
            alpha = "numeric",
            V = "numeric",
            G = "double",
            gamma = "numeric",
            W = "double",            
            theta = "numeric"
  ), 
  body = '
  #include <math.h>
  Rcpp::RNGScope scope;
  Rcpp::NumericVector ycpp(y);
  Rcpp::NumericVector Fcpp(F);
  Rcpp::NumericVector acpp(a);
  Rcpp::NumericVector Rcpp(R);
  Rcpp::NumericVector mcpp(m);
  Rcpp::NumericVector Ccpp(C);
  Rcpp::NumericVector Bcpp(B);
  Rcpp::NumericVector Hcpp(H);
  Rcpp::NumericVector alphacpp(alpha);
  Rcpp::NumericVector Vcpp(V);
  Rcpp::NumericVector thetacpp(theta);
  
  
  int n = ycpp.size();
  double Q;
  double A;
  double f;
  double gammacpp = as<double>(gamma);
  
  double GG = as<double>(G);
  double WW = as<double>(W);
  int nd=1;
  
  
  f    = alphacpp[0]+Fcpp[0]*acpp[0];
  Q    = Rcpp[0]*pow(Fcpp[0],2)+Vcpp[0];
  A    = Rcpp[0]*Fcpp[0]/Q;
  mcpp[0] = acpp[0]+A*(ycpp[0]-f);
  Ccpp[0] = Rcpp[0]-Q*pow(A,2);
  
  
  for (int t = 1; t < n; t++){
  acpp[t] = gammacpp + GG*mcpp[t-1];
  Rcpp[t] = Ccpp[t-1]*pow(GG,2) + WW;
  f    = alphacpp[t]+Fcpp[t]*acpp[t];
  Q    = Rcpp[t]*pow(Fcpp[t],2)+Vcpp[t];
  A    = Rcpp[t]*Fcpp[t]/Q;
  mcpp[t] = acpp[t]+A*(ycpp[t]-f);
  Ccpp[t] = Rcpp[t]-Q*pow(A,2);
  Bcpp[t-1] = Ccpp[t-1]*GG/Rcpp[t];
  Hcpp[t-1] = sqrt(Ccpp[t-1]-Rcpp[t]*pow(Bcpp[t-1],2));
  
  
  }
  
  
  
  thetacpp[n-1] = as<double>(rnorm(nd,mcpp[n-1],sqrt(Ccpp[n-1])));
  for (int t = n-2; t >= 0; t--){
  thetacpp[t] = as<double>(rnorm(nd,mcpp[t]+Bcpp[t]*(thetacpp[t+1]-acpp[t+1]),Hcpp[t]));
  

  
  }
  //Rcpp::Rcout << counts[0]<<" "<< counts[1] << " "<< counts[2] << " "<< counts[3] <<std::endl;
  /*return Rcpp::List::create(Rcpp::Named( "theta" ) = thetacpp,
  Rcpp::Named( "S" ) = Scpp,
  Rcpp::Named( "counts" ) = counts);*/
  return thetacpp;
  ', plugin="Rcpp")




#-------------------------------------------------------
# y    = X*beta + u    u ~ N(0,sig2*I_n)
# beta|sig2 ~ N(b,sig2*inv(A))
#      sig2 ~ IG(v/2,v*lam/2)
#-------------------------------------------------------
fixedpar = function(y,X,A,b,v,lam){
  n     = length(y)
  k     = ncol(X)
  par1  = (v+n)/2
  var   = solve(crossprod(X,X)+A)
  mean  = matrix(var%*%(crossprod(X,y)+crossprod(A,b)),k,1)
  par2  = v*lam + sum((y-crossprod(t(X),mean))^2)
  par2  = (par2 + crossprod(t(crossprod(mean-b,A)),mean-b))/2
  sig2  = 1/rgamma(1,par1,par2)
  var   = var*sig2
  mean  = mean + crossprod(t(chol(var)),rnorm(2))
  return(c(mean,sig2))
}

#----------------------------------------------------------------------------
# Sample Z from 1,2,...,k, with P(Z=i) proportional to q[i]N(mu[i],sig2[i])
#----------------------------------------------------------------------------
ncind = function(y,mu,sig,q){
  w = dnorm(y,mu,sig)*q
  return(sample(1:length(q),size=1,prob=w/sum(w)))
}

#------------------------------------------------
# Quantile statistics
#------------------------------------------------
quant05 = function(x){quantile(x,0.05)}
quant95 = function(x){quantile(x,0.95)}

#------------------------------------------------
# Sampling the log-volatilities in a 
# standard univariate stochastic volatility model
#------------------------------------------------
svu = function(y,lambda,gamma,G,W,a1,R1){
  mu   = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
  sig2 = c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610)
  q    = c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500)
  y    = log(y^2)
  sig  = sqrt(sig2)
  z    = sapply(y-lambda,ncind,mu,sig,q)
  ffbsu(y,1.0,mu[z],sig2[z],G,gamma,W,a1,R1)
}



##############################################################################
#  STOCHASTIC VOLATILITY + FORWARD FILTERING BACKWARD SAMPLING               #
##############################################################################
set.seed(1243)
gamma     = -0.00645
G         = 0.99
W         = 0.15^2
lambda0   = 0.0
n         = 1000
lambda    = rep(0,n)
error     = rnorm(n,0.0,sqrt(W))
lambda[1] = gamma+G*lambda0+error[1]
for (t in 2:n)
  lambda[t] = gamma+G*lambda[t-1]+error[t]
h = exp(lambda/2)
y = rnorm(n,0,h)
plot(y,type="l")

# Prior hyperparameters
b         = c(0,0)
A         = diag(0.00000001,2)
v         = 0.00000002
lam       = 1.00000000
a1        = 0
R1        = 1000

# Initial values
lambdas   = lambda
pars      = c(gamma,G,sqrt(W))
M0        = 2000
M         = 2000

date()
for (iter in 1:(M0+M)){
  lambda  = svu(y,lambda,gamma,G,W,a1,R1)
  X       = cbind(1,lambda[2:n])
  param   = fixedpar(lambda[1:(n-1)],X,A,b,v,lam)
  gamma   = param[1]
  G       = param[2]
  W       = param[3]
  pars    = rbind(pars,c(param[1:2],sqrt(param[3])))
  lambdas = rbind(lambdas,lambda)
}
date()
hs  = exp(lambdas[(M0+1):(M0+M),])
mh  = apply(hs,2,mean)
h05 = apply(hs,2,quant05)
h95 = apply(hs,2,quant95)
hsm   = hs
mhm   = mh
h05m  = h05
h95m  = h95
parsm = pars

names = c("mu","phi","sigma2")
par(mfrow=c(3,2))
for (i in 1:3){
  plot(parsm[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  abline(h=parsm[1,i],col=2)
  hist(parsm[(M0+1):(M0+M),i],xlab="",ylab="",main="")
  abline(v=parsm[1,i],col=2)
}

par(mfrow=c(1,1))
plot(1+h,type="l",ylim=c(0.0,max(1+h95m)),xlab="time",ylab="")
lines(1+h,col=1,lwd=1)
lines(1+mhm,col=2,lwd=1)
lines(1+h05m,col=3,lwd=1)
lines(1+h95m,col=3,lwd=1)
legend(400,4,legend=c("h(t)","E(h(t)|data)","90%CI"),col=c(1,2,3),lty=rep(1,3),lwd=rep(1,3),bty="n")
lines(abs(y)/max(abs(y)),type="h")

# Comparison
par(mfrow=c(1,1))
plot(h,type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
lines(h,col=1,lwd=1.5)
lines(mhn,col=2,lwd=1.5)
lines(mhm,col=3,lwd=1.5)
legend(400,2.5,legend=c("h(t)","normal approximation","mixture normals"),
       col=c(1,2,3),lty=rep(1,3),lwd=rep(1.5,3),bty="n")

# Comparison
par(mfrow=c(2,2))
plot(1:250,h[1:250],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
lines(1:250,h[1:250],col=1,lwd=1.5)
lines(1:250,mhn[1:250],col=2,lwd=1.5)
lines(1:250,mhm[1:250],col=3,lwd=1.5)
plot(251:500,h[251:500],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
lines(251:500,h[251:500],col=1,lwd=1.5)
lines(251:500,mhn[251:500],col=2,lwd=1.5)
lines(251:500,mhm[251:500],col=3,lwd=1.5)
legend(300,2.5,legend=c("h(t)","normal approximation","mixture normals"),
       col=c(1,2,3),lty=rep(1,3),lwd=rep(1.5,3),bty="n")
plot(501:750,h[501:750],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
lines(501:750,h[501:750],col=1,lwd=1.5)
lines(501:750,mhn[501:750],col=2,lwd=1.5)
lines(501:750,mhm[501:750],col=3,lwd=1.5)
plot(501:750,h[751:n],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
lines(501:750,h[751:n],col=1,lwd=1.5)
lines(501:750,mhn[751:n],col=2,lwd=1.5)
lines(501:750,mhm[751:n],col=3,lwd=1.5)





