library(inline)

#-------------------------------------------------------
# Univariate FFBS: 
# y(t)     ~ N(alpha(t)+F(t)*theta(t);V(t))        
# theta(t) ~ N(gamma+G*theta(t-1);W)               
#-------------------------------------------------------
ffbsu = function(y,F,alpha,V,G,gamma,W,a1,R1,S,lambda,p,q,nd=1){
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
#     a[t] = gamma[S[t]+1] + G*m[t-1]
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
  theta = rep(0,n)
#   theta[,n] = rnorm(nd,m[n],sqrt(C[n]))
#   for (t in (n-1):1)
#     theta[,t] = rnorm(nd,m[t]+B[t]*(theta[,t+1]-a[t+1]),H[t])

  ZZZZ = ffbs.cpp(y,F,a,R,m,C,B,H,alpha,V,G,gamma,W,S,theta,lambda,p,q)

#   theta = ZZZZ$theta
#   if (nd==1){
#     theta[1,]
#   }
#   else{
#     theta
#   }
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
            S = "numeric",
            theta = "numeric",
            lambda = "numeric",
            p = "double",
            q = "double"
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
//Rcpp::NumericVector Gcpp(G);
Rcpp::NumericVector gammacpp(gamma);
Rcpp::IntegerVector Scpp(S);
Rcpp::NumericVector lambdacpp(lambda);
Rcpp::NumericVector pred(y);
Rcpp::NumericVector post(y);
Rcpp::NumericVector Pstate(y);
Rcpp::NumericVector thetacpp(theta);


Rcpp::IntegerVector counts(4);
int n = ycpp.size();
double Q;
double A;
double f;
double likelihood_S0;
double likelihood_S1;
double prior = 0.5;
double GG = as<double>(G);
double WW = as<double>(W);
double pp = as<double>(p);
double qq = as<double>(q);
double m;
double sd;
double num;
double den;
int nd=1;


  f    = alphacpp[0]+Fcpp[0]*acpp[0];
  Q    = Rcpp[0]*pow(Fcpp[0],2)+Vcpp[0];
  A    = Rcpp[0]*Fcpp[0]/Q;
  mcpp[0] = acpp[0]+A*(ycpp[0]-f);
  Ccpp[0] = Rcpp[0]-Q*pow(A,2);

post[0] = pp*prior + (1-pp)*(1-prior);

for (int t = 1; t < n; t++){
    acpp[t] = gammacpp[Scpp[t]] + GG*mcpp[t-1];
    Rcpp[t] = Ccpp[t-1]*pow(GG,2) + WW;
    f    = alphacpp[t]+Fcpp[t]*acpp[t];
    Q    = Rcpp[t]*pow(Fcpp[t],2)+Vcpp[t];
    A    = Rcpp[t]*Fcpp[t]/Q;
    mcpp[t] = acpp[t]+A*(ycpp[t]-f);
    Ccpp[t] = Rcpp[t]-Q*pow(A,2);
    Bcpp[t-1] = Ccpp[t-1]*GG/Rcpp[t];
    Hcpp[t-1] = sqrt(Ccpp[t-1]-Rcpp[t]*pow(Bcpp[t-1],2));

    pred[t] = pp*post[t-1] + (1-pp)*(1-post[t-1]);
    m = GG*lambdacpp[t-1];
    sd = sqrt(WW);
    likelihood_S0 = R::dnorm(ycpp[t],0,sqrt(exp(gammacpp[0]+m)),0);
    likelihood_S1 = R::dnorm(ycpp[t],0,sqrt(exp(gammacpp[1]+m)),0);
    post[t] = likelihood_S0 * pred[t] /(likelihood_S0*pred[t]+likelihood_S1*(1-pred[t]));
    //Rcpp::Rcout << likelihood_S0 << " " << likelihood_S1 << " " <<post[t] << std::endl;
  }


  Scpp[n-1] = as<int>(rbinom(1,1,post[n-1]));
  thetacpp[n-1] = as<double>(rnorm(nd,mcpp[n-1],sqrt(Ccpp[n-1])));
  for (int t = n-2; t >= 0; t--){
    thetacpp[t] = as<double>(rnorm(nd,mcpp[t]+Bcpp[t]*(thetacpp[t+1]-acpp[t+1]),Hcpp[t]));
    
      // really hacky code. either 0-p or 1-p, in the 0-p case you get cancellations. 
      num = (Scpp[t+1]-pp)*post[t];
      den = num + (Scpp[t+1]-qq)*(1-post[t]);
      //Rcpp::Rcout << num/den <<std::endl;
      Scpp[t] = as<int>(rbinom(1,1,1-num/den));

      if ( Scpp[t] < Scpp[t+1] ){
          counts[1]++;
      }
      else if ( Scpp[t] > Scpp[t+1] ){
          counts[2]++;
      }
      else if (Scpp[t] == 0){
          counts[0]++;}
      else {
          counts[3]++;
      }
    
  }
//Rcpp::Rcout << counts[0]<<" "<< counts[1] << " "<< counts[2] << " "<< counts[3] <<std::endl;
return Rcpp::List::create(Rcpp::Named( "theta" ) = thetacpp,
                              Rcpp::Named( "S" ) = Scpp,
                              Rcpp::Named( "counts" ) = counts);
//return thetacpp;
', plugin="Rcpp")


#-------------------------------------------------------
# y    = X*beta + u    u ~ N(0,sig2*I_n)
# beta|sig2 ~ N(b,sig2*inv(A))
#      sig2 ~ IG(v/2,v*lam/2)
#-------------------------------------------------------
fixedpar = function(y,X,A,b,v,lam,counts){
  n     = length(y)
  k     = ncol(X)
  par1  = (v+n)/2
  var   = solve(crossprod(X,X)+A)
  mean  = matrix(var%*%(crossprod(X,y)+crossprod(A,b)),k,1)
  par2  = v*lam + sum((y-crossprod(t(X),mean))^2)
  par2  = (par2 + crossprod(t(crossprod(mean-b,A)),mean-b))/2
  sig2  = 1/rgamma(1,par1,par2)
  var   = var*sig2
  mean  = mean + crossprod(t(chol(var)),rnorm(k))
  p = rbeta(1,counts[2]+1,counts[1]+1)
  q = rbeta(1,counts[4]+1,counts[3]+1)
  return(c(sort(mean[1:2]),mean[3],sig2,sort(c(p,q))))
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
svu = function(y,lambda,gamma,G,W,a1,R1,S,p,q){
  mu   = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
  sig2 = c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610)
  pi    = c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500)
  y    = log(y^2)
  sig  = sqrt(sig2)
  z    = sapply(y-lambda,ncind,mu,sig,pi)
  
  return(ffbsu(y,1.0,mu[z],sig2[z],G,gamma,W,a1,R1,S,lambda,p,q))
}

sampleS = function(){
  prior = 0.5
  pred = rep(0,n)
  post = rep(0,n)
  
  p = rbeta(1,1,1)
  q = rbeta(1,1,1)
  
  Pstate[1] = p*prior + q*(1-prior)
  
  # Filter forward
  for (t in 2:n){
    pred[t] = p*Pstate[t-1] + q*(1-Pstate[t-1])
    likelihood_S0 = dnorm(lambda[t],alpha[1]+phi*lambda[t-1],w)
    likelihood_S1 = dnorm(lambda[t],alpha[2]+phi*lambda[t-1],w)
    post[t] = likelihood_S0 * pred[t] /(likelihood_S0 * pred[t]+likelihood_S1*(1-pred[t]))
  }
  
  # backwards sample
  S[n] = rbinom(1,1,post[n])
  for (t in (n-1):2){
      num = trans[1,S[t+1]+1]*post[t]
      den = num + trans[2,S[t+1]+1]*(1-post[t])
      S[t] = rbinom(1,1,num/den)
  }
  
  
  
  trans =  matrix(c(0.3,0.7,0.5,0.5),2,2)
}

##############################################################################
# Defining the mixture of 7 normals that approximate                         #
# the log-chi-square with one degree of freedom                              #
##############################################################################
# set.seed(1576)
mu     = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
sig2   = c(  5.79596, 2.61369, 5.17950,0.16735, 0.64009,0.34023, 1.26261)
q      = c(  0.00730, 0.10556, 0.00002,0.04395, 0.34001,0.24566, 0.25750)
mm     = sum(q*mu)
vv     = sum(q*(sig2+mu^2))-(mm)^2
M      = 10000
x      = seq(-20,10,length=M)
den    = exp(-(x)/2)*exp(-0.5*exp(x))*exp(x)/sqrt(2*pi)
mix    = rep(0,M)
# for (i in 1:M) 
#   mix[i] = sum(q*dnorm(x[i],mu,sqrt(sig2)))
# norm   = dnorm(x,mm,sqrt(vv))
# par(mfrow=c(1,1))
# plot(x,den,ylab="density",main="",type='n')
# lines(x,den,col=1,lty=1,lwd=2)
# lines(x,mix,col=2,lty=2,lwd=2)
# legend(-15,0.2,legend=c("log chi^2_1","mixture of normals"),lty=1:2,col=1:2,lwd=c(2,2),bty="n")


##############################################################################
#  STOCHASTIC VOLATILITY + FORWARD FILTERING BACKWARD SAMPLING               #
##############################################################################
# set.seed(12435)
gamma     = cbind(-.1, .1)
G         = 0.99
W         = .15^2
lambda0   = 0.0
n         = 1000
lambda    = rep(0,n)
error     = rnorm(n,0.0,sqrt(W))
S = rep(0,n)
S[1] = rbinom(1,1,.5)
S[1] = 0
p0 = .3
q0 = .25
lambda[1] = gamma[S[1]+1]+G*lambda0+error[1]
for (t in 2:n){
  if (S[t-1] == 0)
    S[t] = rbinom(1,1,1-p0)
  if (S[t-1] == 1)
    S[t] = rbinom(1,1,1-q0)

  lambda[t] = gamma[S[t]+1]+G*lambda[t-1]+error[t]
}
#print(S)
h = exp(lambda/2)
y = rnorm(n,0,h)
# plot(y,type="l")

# Prior hyperparameters
b         = c(0,0,0)
A         = diag(0.00000001,3)
v         = 0.00000002
lam       = 1.00000000
a1        = 0
R1        = 1000



lambdas   = lambda
pars      = c(gamma,G,sqrt(W),p0,q0)
M0        = 2000
M         = 2000

# Initial values
p = rbeta(1,1,1)
q = rbeta(1,1,1)

print(date())
for (iter in 1:(M0+M)){
  
  output  = svu(y,lambda,gamma,G,W,a1,R1,S,p,q)
  lambda = output$theta
  S = output$S
  counts = output$counts
  X       = cbind((1-S)[2:n],S[2:n],lambda[2:n])
  param   = fixedpar(lambda[1:(n-1)],X,A,b,v,lam,counts)
  alpha0  = param[1]
  alpha1  = param[2]
  G       = param[3]
  W       = param[4]
  p       = param[5]
  q       = param[6]
  pars    = rbind(pars,c(param[1:3],sqrt(param[4]),param[5:6]))
  lambdas = rbind(lambdas,lambda)
}
print(date())
hs  = exp(lambdas[(M0+1):(M0+M),])
mh  = apply(hs,2,mean)
h05 = apply(hs,2,quant05)
h95 = apply(hs,2,quant95)
hsm   = hs
mhm   = mh
h05m  = h05
h95m  = h95
parsm = pars

names = c("alpha0","alpha1","phi","sigma2","p","q")
par(mfrow=c(3,2))
for (i in 1:3){
  plot(parsm[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  abline(h=parsm[1,i],col=2)
  hist(parsm[(M0+1):(M0+M),i],xlab="",ylab="",main="")
  abline(v=parsm[1,i],col=2)
  print(names[i])
  print(mean(parsm[(M0+1):(M0+M),i]))
}
par(mfrow=c(3,2))
for (i in 4:6){
  plot(parsm[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  abline(h=parsm[1,i],col=2)
  hist(parsm[(M0+1):(M0+M),i],xlab="",ylab="",main="")
  abline(v=parsm[1,i],col=2)
  print(names[i])
  print(mean(parsm[(M0+1):(M0+M),i]))
}
# 
# par(mfrow=c(1,1))
# plot(1+h,type="l",ylim=c(0.0,max(1+h95m)),xlab="time",ylab="")
# lines(1+h,col=1,lwd=1)
# lines(1+mhm,col=2,lwd=1)
# lines(1+h05m,col=3,lwd=1)
# lines(1+h95m,col=3,lwd=1)
# legend(400,4,legend=c("h(t)","E(h(t)|data)","90%CI"),col=c(1,2,3),lty=rep(1,3),lwd=rep(1,3),bty="n")
# lines(abs(y)/max(abs(y)),type="h")
# 
# # Comparison
# par(mfrow=c(1,1))
# plot(h,type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
# lines(h,col=1,lwd=1.5)
# lines(mhn,col=2,lwd=1.5)
# lines(mhm,col=3,lwd=1.5)
# legend(400,2.5,legend=c("h(t)","normal approximation","mixture normals"),
# col=c(1,2,3),lty=rep(1,3),lwd=rep(1.5,3),bty="n")
# 
# # Comparison
# par(mfrow=c(2,2))
# plot(1:250,h[1:250],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
# lines(1:250,h[1:250],col=1,lwd=1.5)
# lines(1:250,mhn[1:250],col=2,lwd=1.5)
# lines(1:250,mhm[1:250],col=3,lwd=1.5)
# plot(251:500,h[251:500],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
# lines(251:500,h[251:500],col=1,lwd=1.5)
# lines(251:500,mhn[251:500],col=2,lwd=1.5)
# lines(251:500,mhm[251:500],col=3,lwd=1.5)
# legend(300,2.5,legend=c("h(t)","normal approximation","mixture normals"),
# col=c(1,2,3),lty=rep(1,3),lwd=rep(1.5,3),bty="n")
# plot(501:750,h[501:750],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
# lines(501:750,h[501:750],col=1,lwd=1.5)
# lines(501:750,mhn[501:750],col=2,lwd=1.5)
# lines(501:750,mhm[501:750],col=3,lwd=1.5)
# plot(501:750,h[751:n],type="l",ylim=c(0.0,max(h,mhn,mhm)),xlab="time",ylab="")
# lines(501:750,h[751:n],col=1,lwd=1.5)
# lines(501:750,mhn[751:n],col=2,lwd=1.5)
# lines(501:750,mhm[751:n],col=3,lwd=1.5)
# 
# 
# 
# 
# 

