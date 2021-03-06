### Final Model - WTI Crude Estimates


### Clear Data and Start Timer
rm(list = ls())
set.seed(17)
start.time = proc.time()

### Load Libraries
library(Matrix)
library(MASS)
library(rgl)
library(pscl)
library(inline)
library(Rcpp)

### Load Data
normalize <- function(x) {
  x = x - mean(x)
  vec.range = max(x) - min(x)
  x = x / vec.range
} 

normalizeUnit <- function(x) {
  x = x - min(x)
  vec.max = max(x)
  x = x/vec.max
}

data = read.csv('WTI3.csv')
data = data[9:nrow(data),]
y = data[,2:13] 

dates = as.Date(data[,1], format = '%m/%d/%Y')

days = weekdays(dates)
days.wednesday = rep(0, length(days))
for(i in 1:length(days)){
  if(days[i] == 'Wednesday'){
    days.wednesday[i] = 1
  }
}
days.wednesday = which(days.wednesday == 1)

y.spread = matrix(0, nrow(y), (ncol(y)-1)) 
for(i in 1:nrow(y)){
  for(j in 1:(ncol(y)-1)){
    y.spread[i, j] = y[i, (j+1)] - y[i, (j)] - 0.27
  }
}
y.cumspread = t(apply(y.spread, 1, cumsum))

### Create Y and Z matrices

# Y52
nx.test = 11
nt.test = 52
x.test = 1:nx.test
t.test = 1:nt.test
  
nx.train = 11
nt.train = 52
x.train = 1:nx.train
t.train = 1:nt.train + nt.test
t.train.plot = 1:nt.train

y.train = y.cumspread[days.wednesday[t.train], x.train]
y.test = y.cumspread[days.wednesday[t.test], x.test]

surface3d(t.test, x.test, y.test, col='green', alpha = .5)
surface3d(t.test, x.test, y.test, col='green', alpha = .5, front = "line", back = "line")
surface3d(t.train, x.train, y.train, col='green', alpha = .5)
surface3d(t.train, x.train, y.train, col='green', alpha = .5, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='green', alpha = .2, front = "line", back = "line")

t.vec.train = rep(1:nt.train, nx.train)
x.vec.train = NULL
for(i in 1:nx.train){
  x.vec.train = c(x.vec.train, rep(i,nt.train))
}

y.train.vec = as.vector(y.train)

t.vec.test = rep(1:nt.test, nx.test)
x.vec.test = NULL
for(i in 1:nx.test){
  x.vec.test = c(x.vec.test, rep(i,nt.test))
}

y.test.vec = as.vector(y.test)

# Economic Data
X1 = rev(read.csv('Import-Export.csv')[,2])
X1.train = X1[t.train]
X1.test = X1[t.test]
X2 = rev(read.csv('Demand.csv')[,2])
X2.train = X2[t.train]
X2.test = X2[t.test]
X3 = rev(read.csv('Stocks.csv')[,3])
X3.train = X3[t.train]
X3.test = X3[t.test]
X4 = rev(read.csv('refining.csv')[,3])
X4.train = X4[t.train]
X4.test = X4[t.test]

Z.train.class = cbind(normalize(X1.train), normalize(X2.train), normalize(X3.train), normalize(X4.train))

S.train = rep(0, nt.train)
for(i in 1:nrow(y.train)){
  if(y.train[i,1] > max(y.train[i,6:nx.train])){
    S.train[i] = 1
  }
}
St.train = rep(S.train, nx.train)
St.cont = which(St.train == 0)
St.back = which(St.train == 1)

p = 1 + 4*nx.train

Z.train = matrix(0,nt.train*nx.train, p)

Z.train[,1] = rep(1, nt.train*nx.train)

for(i in 1:nx.train){
  ind = which(x.vec.train == i)
  Z.train[ind,(i+1)] = normalize(X1.train)
}

for(i in 1:nx.train){
  ind = which(x.vec.train == i)
  Z.train[ind,(i+ 1 + nx.train)] = normalize(X2.train)
}

for(i in 1:nx.train){
  ind = which(x.vec.train == i)
  Z.train[ind,(i+1 + 2*nx.train)] = normalize(X3.train)
}

for(i in 1:nx.train){
  ind = which(x.vec.train == i)
  Z.train[ind,(i+1 +3* nx.train)] = normalize(X4.train)
}

Z.test = matrix(0,nt.test*nx.test, p)

Z.test[,1] = rep(1, nt.test*nx.test)

for(i in 1:nx.test){
  ind = which(x.vec.test == i)
  Z.test[ind,(i+1)] = normalize(X1.test)
}

for(i in 1:nx.test){
  ind = which(x.vec.test == i)
  Z.test[ind,(i+ 1 + nx.test)] = normalize(X2.test)
}

for(i in 1:nx.test){
  ind = which(x.vec.test == i)
  Z.test[ind,(i+1 + 2*nx.test)] = normalize(X3.test)
}

for(i in 1:nx.test){
  ind = which(x.vec.test == i)
  Z.test[ind,(i+1 +3* nx.test)] = normalize(X4.test)
}

pdf('Plot_Data.pdf', width = 12, height = 9)
par(mfrow = c(2,1))
plot(normalizeUnit(y.train[,1]), typ = 'l', lwd = 3, main = 'Front Month')
lines(normalizeUnit(X1.train), col = 2, lty = 2)
lines(normalizeUnit(X2.train), col = 4, lty = 3)
lines(normalizeUnit(X3.train), col = 6, lty = 4)
lines(normalizeUnit(X4.train), col = 8, lty = 5)

legend(45, .5, c("Actual", "Import", "Demand", "Stocks", "Refining"), col = c(1,2,4,6,8), lwd = c(3,0,0,0,0),lty = c(1,2,3,4,5), box.col = "white",bg = "white")

plot(normalizeUnit(y.train[,11]), typ = 'l', lwd = 3, main = 'Back Month')
lines(normalizeUnit(X1.train), col = 2, lty = 2)
lines(normalizeUnit(X2.train), col = 4, lty = 3)
lines(normalizeUnit(X3.train), col = 6, lty = 4)
lines(normalizeUnit(X4.train), col = 8, lty = 5)
dev.off()



# if(tcpp[i] == tcpp[j]){
#   Kcpp[nrows * j + i] = Kcpp[nrows * j + i] + parcpp[2] * exp(-0.5 * ((xcpp[i] - xcpp[j])/parcpp[3]) * ((xcpp[i] - xcpp[j])/parcpp[3]));
# }

### Define Functions

K_Calc.cpp = cxxfunction(signature(par = "numeric", x = "numeric", t= "numeric", K = "numeric"), body = '
Rcpp::NumericMatrix Kcpp(K);
Rcpp::NumericVector xcpp(x);
Rcpp::NumericVector tcpp(t);
Rcpp::NumericVector parcpp(par);

int nrows = Kcpp.nrow();
int ncolumns = Kcpp.ncol();

for (int i = 0; i < nrows; i++){
    for (int j = 0; j < ncolumns; j++){
        if(xcpp[i] == xcpp[j]){
            Kcpp[nrows * j + i] = Kcpp[nrows * j + i] + parcpp[0] * exp(-0.5 * ((tcpp[i] - tcpp[j])/parcpp[1]) * ((tcpp[i] - tcpp[j])/parcpp[1]));
        }

        if(tcpp[i] == tcpp[j]){
          Kcpp[nrows * j + i] = Kcpp[nrows * j + i] + parcpp[2] * exp(-0.5 * ((xcpp[i] - xcpp[j])/parcpp[3]) * ((xcpp[i] - xcpp[j])/parcpp[3]));
        }
    }
}
return Kcpp;
', plugin="Rcpp")

log.like = function(par, y, nx, nt, x.vec, t.vec, sig2, K){
  theta.t = par[1:2]
  theta.x = par[3:4]
  sig2.diag = diag(1/sig2, nx*nt, nx*nt)
  K = matrix(0, nrow(K), ncol(K))
  K = K_Calc.cpp(par, x.vec, t.vec, K) + sig2.diag
#   K.inv = ginv(K, 1e-20)
  K.inv = solve(K)
  K.det = det(K)
  ll = .5 * t(y) %*% K.inv %*% y + 0.5 * log(K.det) + nx*nt/2*log(2*pi)
  return(ll)
}

GPparam.opt = function(par, y.GP, nx, nt, x.vec, t.vec, sig2, K){
  hyperparam = par
  par.opt = nlminb(start = hyperparam, log.like, y = y.GP, nx = nx, nt = nt, x.vec = x.vec, t.vec = t.vec, sig2 = sig2, K = K, lower = rep(0,6), control = list(iter.max = 10, trace = 0, rel.tol = .01))
  print('Optimizing Done')
  theta.t = par.opt$par[1:2]
  theta.x = par.opt$par[3:4]
  par = c(theta.t, theta.x)
  return(par)
}

RSSB.lasso = function(par, y, z, lam){
  B = par
  B.mat = as.matrix(B, 1, p)
  SSE = sum((y - z %*% B.mat)**2)
  pen = sqrt(sum(B**2))
  obj = SSE + lam*pen
}

### Setup Model Parameters

### Non-split regression
B.full = ginv(t(Z.train) %*% Z.train, 1e-20) %*% t(Z.train) %*% y.train.vec
y.full = Z.train %*% B.full

MSE.full = sum((y.train.vec - y.full)**2)

surface3d(t.train, x.train, y.train, col='green', alpha = .4)
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.full), nt.train, nx.train), col='orange', alpha = .4)


### Least Squares
B.cont = ginv(t(Z.train[St.cont,]) %*% Z.train[St.cont,], 1e-20) %*% t(Z.train[St.cont,]) %*% y.train.vec[St.cont]
B.back = ginv(t(Z.train[St.back,]) %*% Z.train[St.back,], 1e-20) %*% t(Z.train[St.back,]) %*% y.train.vec[St.back]

y.fromB = matrix(0,nx.train*nt.train,1)
y.fromB[St.cont] = Z.train[St.cont,] %*% B.cont
y.fromB[St.back] = Z.train[St.back,] %*% B.back

MSE.fromB = sum((y.train.vec - y.fromB)**2)


### Gaussian Process
K.train = matrix(0,nx.train*nt.train,nx.train*nt.train)
sig2 = 1
sig2.diag = diag(1/sig2, nx.train*nt.train, nx.train*nt.train)
par.GP = c(10, 5, 1, 5)
y.forGP = y.train.vec - y.fromB

par.opt = GPparam.opt(par.GP, as.vector(y.forGP), nx.train, nt.train, x.vec.train, t.vec.train, sig2, K.train)
par.GP = par.opt

K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GP = matrix(0, nrow(K), ncol(K))
K.GP = K_Calc.cpp(par.GP, x.vec.train, t.vec.train, K.GP)    
E = ginv(sig2.diag  + ginv(K.GP, 1e-20), 1e-20)
y.fromGP = E %*% sig2.diag %*% y.forGP

y.all = y.fromB + y.fromGP
MSE.all = sum((y.train.vec - y.all)**2)

### Just Gaussian Process
sig2 = 1
sig2.diag = diag(1/sig2, nx.train*nt.train, nx.train*nt.train)
K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GP = matrix(0, nrow(K), ncol(K))
K.GP = K_Calc.cpp(par.GP, x.vec.train, t.vec.train, K.GP) + sig2.diag
K.GP.inv = ginv(K.GP, 1e-20)
E = ginv(sig2.diag  + ginv(K.GP, 1e-20), 1e-20)
B.GP = solve(t(Z.train) %*% K.GP.inv %*% Z.train) %*% t(Z.train) %*% K.GP.inv %*% y.train.vec

y.fromGPOnly = Z.train %*% B.GP + E %*% sig2.diag %*% (y.train.vec - Z.train %*% B.GP)
  
MSE.GPOnly = sum((y.train.vec - y.fromGPOnly)**2)
# 
# surface3d(t.train, x.train, y.train, col='green', alpha = .4)
# surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.fromGPOnly), nt.train, nx.train), col='orange', alpha = .4)
#   
# surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.fromGPOnly), nt.train, nx.train), col='orange', alpha = .6)
# surface3d(t.train, x.train, matrix(as.vector(y.full), nt.train, nx.train), col='blue', alpha = .4)
# 
# pdf('Plot_GP_Full_BackMonth.pdf', width = 12, height = 9)
# par(mfrow = c(2,1))
# plot(y.train[,11], typ = 'l', lty = 2, main = 'Back Month - GP vs Full', xlab = '', ylab = 'Spread')
# lines(matrix(as.vector(y.full), nt.train, nx.train)[,11], col = 'red')
# lines(matrix(as.vector(y.fromGPOnly), nt.train, nx.train)[,11], col = 'blue')
# 
# plot(matrix(as.vector(y.full), nt.train, nx.train)[,11] - y.train[,11], typ = 'l', , col = 'red', main = 'Back Month - GP vs Full Error', xlab = 'Historical Time', ylab = 'Error')
# lines(matrix(as.vector(y.fromGPOnly), nt.train, nx.train)[,11] - y.train[,11], col = 'blue')
# dev.off()
# 
# ### MSE Summary
MSE.full
MSE.fromB
# MSE.fromB.gibbs
# MSE.fromB.lasso
MSE.all
MSE.GPOnly
# 


# pdf('Plot_Coeffs_GP_Full.pdf', width = 12, height = 9)
# par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
# plot(B.GP[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta', ylim = c(-1.5, 1.4))
# points(B.full[2:12], col = 'red', pch = 3)
# abline(h = 0, col = 'grey')
# plot(B.GP[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
# points(B.full[13:23], col = 'red', pch = 3)
# abline(h = 0, col = 'grey')
# plot(B.GP[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
# points(B.full[24:34], col = 'red', pch = 3)
# abline(h = 0, col = 'grey')
# plot(B.GP[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
# points(B.full[35:45], col = 'red', pch = 3)
# abline(h = 0, col = 'grey')
# mtext("Coeffcients - GP Only vs Full", outer = TRUE, cex = 1.5)
# dev.off()
# 
# surface3d(t.train, x.train, y.train, col='green', alpha = .4)
# surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.fromB), nt.train, nx.train), col='orange', alpha = .4)
# 
# surface3d(t.train, x.train, y.train, col='green', alpha = .4)
# surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.fromB.lasso), nt.train, nx.train), col='blue', alpha = .4)
# 
# surface3d(t.train, x.train, y.train, col='green', alpha = .4)
# surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.fromB.gibbs), nt.train, nx.train), col='blue', alpha = .4)
# 
# surface3d(t.train, x.train, y.train, col='green', alpha = .5)
# surface3d(t.train, x.train, y.train, col='green', alpha = .5, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
# surface3d(t.train, x.train, matrix(as.vector(y.all), nt.train, nx.train), col='orange', alpha = .5)
# 
# ### Next Step is in contango
# 
<<<<<<< HEAD
# ### Predict Bayes
=======
### Predict Bayes
>>>>>>> 44e11bc49c8dc69ff008fa58d04599305c376aad

y.pred = y.test[52,]
z.pred = Z.test[which(t.vec.test == 52),]


# Full
y.pred.full = z.pred %*% B.full
y.pred.B = z.pred %*% B.cont
# y.pred.gibbs = z.pred %*% B.cont.gibbs
# y.pred.lasso = z.pred %*% B.cont.lasso
y.pred.GP = z.pred %*% B.GP
<<<<<<< HEAD
# 
# ### Prediction linear plus GP
=======

### Prediction linear plus GP
>>>>>>> 44e11bc49c8dc69ff008fa58d04599305c376aad
y.forGP = y.train.vec - y.fromB
K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GPFull = matrix(0, nrow(K), ncol(K))
par.GPFull = c(10, 5, 1, 5)
sig2 = 1
sig2.diag = diag(1/sig2, nx.train*nt.train, nx.train*nt.train)
par.opt = GPparam.opt(par.GPFull, y.forGP, nx.train, nt.train, x.vec.train, t.vec.train, sig2, K.GPFull)
par.GPFull = par.opt

K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GP = matrix(0, nrow(K), ncol(K))
K.GP = K_Calc.cpp(par.GPFull, x.vec.train, t.vec.train, K.GP) + sig2.diag
K.GP.inv = ginv(K.GP, 1e-20)
B.GP = solve(t(Z.train) %*% K.GP.inv %*% Z.train) %*% t(Z.train) %*% K.GP.inv %*% y.forGP

g.star.GPadd = rep(0,11)

for(i in 1:11){
  x.vec.pred = c(x.vec.train, (1:11)[-i])
  t.vec.pred = c(t.vec.train, rep(53,10))
  sig2 = 1
  sig2.diag = diag(1/sig2, 582, 582)
  K = matrix(0,582,582)
  K.GP = matrix(0, nrow(K), ncol(K))
  K.GP = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.GP) + sig2.diag
  K.GP.inv = ginv(K.GP, 1e-20)
  
  x.vec.pred = c(x.vec.train, (1:11)[-i], (1:11)[i])
  t.vec.pred = c(t.vec.train, rep(53,11))
  
  K.pre = matrix(0, (583),(583))
  K.pre = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.pre)
  
  K.star = as.matrix(K.pre[,583][1:(582)])
  
  y.pred.vec = c(y.forGP, matrix(y.forGP, nt.train, nx.train)[1,][-i])
  f.star = t(K.star) %*% (K.GP.inv) %*% y.pred.vec
  H.star = Z.test[which(t.vec.test == 52)[i],]
  Z.pred = rbind(Z.train, Z.test[which(t.vec.test == 52)[-i],])
  
  R.pred = H.star - t(Z.pred) %*% K.GP.inv %*% K.star
  g.star.GPadd[i] = f.star + t(R.pred) %*% B.GP
}

y.pred.BGP = y.pred.B + g.star.GPadd
points(y.pred.B + g.star.GPadd)
<<<<<<< HEAD
# 
# ### Prediction just GP
# K = matrix(0,nx.train*nt.train,nx.train*nt.train)
# K.GPFull = matrix(0, nrow(K), ncol(K))
# par.GPFull = c(10, 5, 1, 5)
# sig2 = 1
# sig2.diag = diag(1/sig2, nx.train*nt.train, nx.train*nt.train)
# par.opt = GPparam.opt(par.GPFull, y.train.vec, nx.train, nt.train, x.vec.train, t.vec.train, sig2, K.GPFull)
# par.GPFull = par.opt
# 
# K = matrix(0,nx.train*nt.train,nx.train*nt.train)
# K.GP = matrix(0, nrow(K), ncol(K))
# K.GP = K_Calc.cpp(par.GPFull, x.vec.train, t.vec.train, K.GP) + sig2.diag
# K.GP.inv = ginv(K.GP, 1e-20)
# B.GP = solve(t(Z.train) %*% K.GP.inv %*% Z.train) %*% t(Z.train) %*% K.GP.inv %*% y.train.vec
# 
# 
# 
# 
# 
# g.star = rep(0,11)
# for(i in 1:11){
#   x.vec.pred = c(x.vec.train, (1:11)[-i])
#   t.vec.pred = c(t.vec.train, rep(53,10))
#   sig2 = 1
#   sig2.diag = diag(1/sig2, 582, 582)
#   K = matrix(0,582,582)
#   K.GP = matrix(0, nrow(K), ncol(K))
#   K.GP = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.GP) + sig2.diag
#   K.GP.inv = ginv(K.GP, 1e-20)
#   
#   x.vec.pred = c(x.vec.train, (1:11)[-i], (1:11)[i])
#   t.vec.pred = c(t.vec.train, rep(53,11))
#   
#   K.pre = matrix(0, (583),(583))
#   K.pre = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.pre)
#   
#   K.star = as.matrix(K.pre[,583][1:(582)])
#   
#   y.pred.vec = c(y.train.vec, y.train[1,][-i])
#   f.star = t(K.star) %*% (K.GP.inv) %*% y.pred.vec
#   H.star = Z.test[which(t.vec.test == 52)[i],]
#   Z.pred = rbind(Z.train, Z.test[which(t.vec.test == 52)[-i],])
#   
#   R.pred = H.star - t(Z.pred) %*% K.GP.inv %*% K.star
#   g.star[i] = f.star + t(R.pred) %*% B.GP
# }
# pdf('Plot_Predict.pdf', width = 12, height = 9)
# plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction', xlab = 'Month', ylab = 'Spread' )
# points(y.test[52,])
# lines(y.pred.full, col = 2)
# points(y.pred.full, pch = 2, col = 2)
# lines(y.pred.B, col = 3)
# points(y.pred.B, pch = 3, col = 3)
# lines(y.pred.gibbs, col = 4)
# points(y.pred.gibbs, pch = 4, col = 4)
# lines(y.pred.lasso, col = 5)
# points(y.pred.lasso, pch = 5, col = 5)
# lines(y.pred.BGP, col = 6)
# points(y.pred.BGP, pch = 6, col = 6)
# lines(g.star, col = 7)
# points(g.star, pch = 7, col = 7)
# legend(1, -1, c("Actual", "Full", "LS", "Gibbs", "Lasso", "LS-GP", "GP"), col = c(1,2,3,4,5,6,7), pch = c(1,2,3,4,5,6,7), lwd = c(3,0,0,0,0,0,0))
# dev.off()
# 
# pdf('Plot_Predict_GPs.pdf', width = 12, height = 9)
# plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction - LS/GP vs GP', xlab = 'Month', ylab = 'Spread' )
# points(y.test[52,])
# lines(y.pred.BGP, col = 6)
# points(y.pred.BGP, pch = 6, col = 6)
# lines(g.star, col = 7)
# points(g.star, pch = 7, col = 7)
# legend(1, -1, c("Actual", "LS-GP", "GP"), col = c(1,6,7), pch = c(1,6,7), lwd = c(3,0,0))
# dev.off()
# 
# prederror.full = sum((y.pred.full - y.test[52,])**2)
# prederror.LS = sum((y.pred.B - y.test[52,])**2)
=======

### Prediction just GP
K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GPFull = matrix(0, nrow(K), ncol(K))
par.GPFull = c(10, 5, 1, 5)
sig2 = 1
sig2.diag = diag(1/sig2, nx.train*nt.train, nx.train*nt.train)
par.opt = GPparam.opt(par.GPFull, y.train.vec, nx.train, nt.train, x.vec.train, t.vec.train, sig2, K.GPFull)
par.GPFull = par.opt

K = matrix(0,nx.train*nt.train,nx.train*nt.train)
K.GP = matrix(0, nrow(K), ncol(K))
K.GP = K_Calc.cpp(par.GPFull, x.vec.train, t.vec.train, K.GP) + sig2.diag
K.GP.inv = ginv(K.GP, 1e-20)
B.GP = solve(t(Z.train) %*% K.GP.inv %*% Z.train) %*% t(Z.train) %*% K.GP.inv %*% y.train.vec


g.star = rep(0,11)
for(i in 1:11){
  x.vec.pred = c(x.vec.train, (1:11)[-i])
  t.vec.pred = c(t.vec.train, rep(53,10))
  sig2 = 1
  sig2.diag = diag(1/sig2, 582, 582)
  K = matrix(0,582,582)
  K.GP = matrix(0, nrow(K), ncol(K))
  K.GP = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.GP) + sig2.diag
  K.GP.inv = ginv(K.GP, 1e-20)
  
  x.vec.pred = c(x.vec.train, (1:11)[-i], (1:11)[i])
  t.vec.pred = c(t.vec.train, rep(53,11))
  
  K.pre = matrix(0, (583),(583))
  K.pre = K_Calc.cpp(par.GPFull, x.vec.pred, t.vec.pred, K.pre)
  
  K.star = as.matrix(K.pre[,583][1:(582)])
  
  y.pred.vec = c(y.train.vec, y.train[1,][-i])
  f.star = t(K.star) %*% (K.GP.inv) %*% y.pred.vec
  H.star = Z.test[which(t.vec.test == 52)[i],]
  Z.pred = rbind(Z.train, Z.test[which(t.vec.test == 52)[-i],])
  
  R.pred = H.star - t(Z.pred) %*% K.GP.inv %*% K.star
  g.star[i] = f.star + t(R.pred) %*% B.GP
}

plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction', xlab = 'Month', ylab = 'Spread' )
points(y.test[52,])
lines(y.pred.full, col = 2)
points(y.pred.full, pch = 2, col = 2)
lines(y.pred.B, col = 3)
points(y.pred.B, pch = 3, col = 3)
lines(y.pred.BGP, col = 6)
points(y.pred.BGP, pch = 6, col = 6)
lines(g.star, col = 7)
points(g.star, pch = 7, col = 7)
legend(1, -1, c("Actual", "Full", "LS", "LS-GP", "GP"), col = c(1,2,3,4,5,6,7), pch = c(1,2,3,4,5,6,7), lwd = c(3,0,0,0,0,0,0))



plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction - LS/GP vs GP', xlab = 'Month', ylab = 'Spread' )
points(y.test[52,])
lines(y.pred.BGP, col = 6)
points(y.pred.BGP, pch = 6, col = 6)
lines(g.star, col = 7)
points(g.star, pch = 7, col = 7)
legend(1, -1, c("Actual", "LS-GP", "GP"), col = c(1,6,7), pch = c(1,6,7), lwd = c(3,0,0))


prederror.full = sum((y.pred.full - y.test[52,])**2)
prederror.LS = sum((y.pred.B - y.test[52,])**2)
>>>>>>> 44e11bc49c8dc69ff008fa58d04599305c376aad
# prederror.gibbs = sum((y.pred.gibbs - y.test[52,])**2)
# prederror.lasso = sum((y.pred.lasso - y.test[52,])**2)
prederror.BGP = sum((y.pred.BGP - y.test[52,])**2)
prederror.GP= sum((g.star - y.test[52,])**2)

prederror = c(prederror.full,prederror.LS, prederror.BGP, prederror.GP)

par(mfrow = c(1,1))
lab = c('Full', 'LS', 'LS - GP', 'GP')
barplot(prederror, axisnames = TRUE, names.arg = lab, main = 'Prediction Error')




