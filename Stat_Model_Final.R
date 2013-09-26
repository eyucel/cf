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
  par.opt = nlminb(start = hyperparam, log.like, y = y.GP, nx = nx, nt = nt, x.vec = x.vec, t.vec = t.vec, sig2 = sig2, K = K, lower = rep(0,6), control = list(iter.max = 10, trace = 1, rel.tol = .01))
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


### Gibbs Regression Estimation
M = 5000
at2 = 1
bt2 = 1
t2.cont = rigamma(1, alpha = at2/2, beta = bt2/2)
t2.back = rigamma(1, alpha = at2/2, beta = bt2/2)
B.cont.gibbs.track = matrix(0,p,M)
t2.cont.gibbs.track = rep(0,M)
B.back.gibbs.track = matrix(0,p,M)
t2.back.gibbs.track = rep(0,M)

for(i in 1:M){
  # B
  varB.cont = solve(t(Z.train[St.cont,]) %*% Z.train[St.cont,] + diag(1/t2.cont, p))
  meanB.cont = varB.cont %*% t(Z.train[St.cont,]) %*% y.train.vec[St.cont]
  B.cont.gibbs = chol(varB.cont) %*% rnorm(p) + meanB.cont
  
  # t2
  t2.cont.gibbs = rigamma(1, alpha = (p + at2)/2, beta = 0.5 * (t(B.cont.gibbs) %*% B.cont.gibbs + bt2))
  
  B.cont.gibbs.track[,i] = B.cont.gibbs
  t2.cont.gibbs.track[i] = t2.cont.gibbs
  
  # B
  varB.back = solve(t(Z.train[St.back,]) %*% Z.train[St.back,] + diag(1/t2.back, p))
  meanB.back = varB.back %*% t(Z.train[St.back,]) %*% y.train.vec[St.back]
  B.back.gibbs = chol(varB.back) %*% rnorm(p) + meanB.back
  
  # t2
  t2.back.gibbs = rigamma(1, alpha = (p + at2)/2, beta = 0.5 * (t(B.back.gibbs) %*% B.back.gibbs + bt2))
  
  B.back.gibbs.track[,i] = B.back.gibbs
  t2.back.gibbs.track[i] = t2.back.gibbs
  
  if(i %% 100 == 0){print(i)}
}

B.cont.gibbs.mean  = apply(B.cont.gibbs.track[,2500:M], 1, mean)
t2.cont.gibbs.mean  = mean(t2.cont.gibbs.track[2500:M])
B.back.gibbs.mean  = apply(B.back.gibbs.track[,2500:M], 1, mean)
t2.back.gibbs.mean  = mean(t2.back.gibbs.track[2500:M])

y.fromB.gibbs = matrix(0,nx.train*nt.train,1)
y.fromB.gibbs[St.cont] = Z.train[St.cont,] %*% B.cont.gibbs.mean
y.fromB.gibbs[St.back] = Z.train[St.back,] %*% B.back.gibbs.mean

MSE.fromB.gibbs = sum((y.train.vec - y.fromB.gibbs)**2)

### Least Squares
B.cont = ginv(t(Z.train[St.cont,]) %*% Z.train[St.cont,], 1e-20) %*% t(Z.train[St.cont,]) %*% y.train.vec[St.cont]
B.back = ginv(t(Z.train[St.back,]) %*% Z.train[St.back,], 1e-20) %*% t(Z.train[St.back,]) %*% y.train.vec[St.back]

y.fromB = matrix(0,nx.train*nt.train,1)
y.fromB[St.cont] = Z.train[St.cont,] %*% B.cont
y.fromB[St.back] = Z.train[St.back,] %*% B.back

MSE.fromB = sum((y.train.vec - y.fromB)**2)


### Lasso Regression
p.lasso = 45
par.lasso = rep(0, p.lasso)
lam.lasso = seq(.1,5.1,.2)
K.fold = 6

# Contango
y.lasso.cont = y.train.vec[St.cont]
z.lasso.cont = Z.train[St.cont,]
n.train = length(y.train.vec[St.cont])
ind.mat = matrix(0, K.fold, n.train/K.fold)
sample.lasso.full = 1:n.train
sample.lasso = sample.lasso.full

for(i in 1:K.fold){
  ind.mat[i,] = sample(x = sample.lasso, size = n.train/K.fold)
  sample.lasso = sample.lasso.full[-ind.mat[1:i,]]
}

MSE.lasso = matrix(0, K.fold, length(lam.lasso))

for(i in 1:K.fold){
  for(j in 1:length(lam.lasso)){
    par.opt.lasso = nlminb(start = par.lasso, RSSB.lasso, y = y.lasso.cont[ind.mat[i,]], z = z.lasso.cont[ind.mat[i,],], lam = lam.lasso[j], control = list(iter.max = 100, trace = 1))
    B.lasso = par.opt.lasso$par
    
    MSE.lasso[i,j] = mean((y.lasso.cont[-ind.mat[i,]] - z.lasso.cont[-ind.mat[i,],] %*% as.matrix(B.lasso))**2)
#     print(i)
  }
}

MSE.lasso = apply(MSE.lasso, 2, sum)

par(mfrow = c(1,1))
plot(lam.lasso, MSE.lasso)
lam.min = lam.lasso[which.min(MSE.lasso)]
par.opt.lasso = nlminb(start = par.lasso, RSSB.lasso, y = y.lasso.cont, z = z.lasso.cont, lam = lam.min, control = list(iter.max = 100, trace = 1))
B.cont.lasso = par.opt.lasso$par

# Backwardation
y.lasso.back = y.train.vec[St.back]
z.lasso.back = Z.train[St.back,]
n.train = length(y.train.vec[St.back])
ind.mat = matrix(0, K.fold, n.train/K.fold)
sample.lasso.full = 1:n.train
sample.lasso = sample.lasso.full

for(i in 1:K.fold){
  ind.mat[i,] = sample(x = sample.lasso, size = n.train/K.fold)
  sample.lasso = sample.lasso.full[-ind.mat[1:i,]]
}

MSE.lasso = matrix(0, K.fold, length(lam.lasso))

for(i in 1:K.fold){
  for(j in 1:length(lam.lasso)){
    par.opt.lasso = nlminb(start = par.lasso, RSSB.lasso, y = y.lasso.back[ind.mat[i,]], z = z.lasso.back[ind.mat[i,],], lam = lam.lasso[j], control = list(iter.max = 100, trace = 1))
    B.lasso = par.opt.lasso$par
    
    MSE.lasso[i,j] = mean((y.lasso.back[-ind.mat[i,]] - z.lasso.back[-ind.mat[i,],] %*% as.matrix(B.lasso))**2)
    print(i)
  }
}

MSE.lasso = apply(MSE.lasso, 2, sum)

par(mfrow = c(1,1))
plot(lam.lasso, MSE.lasso)
lam.min = lam.lasso[which.min(MSE.lasso)]
par.opt.lasso = nlminb(start = par.lasso, RSSB.lasso, y = y.lasso.back, z = z.lasso.back, lam = lam.min, control = list(iter.max = 100, trace = 1))
B.back.lasso = par.opt.lasso$par

y.fromB.lasso = matrix(0,nx.train*nt.train,1)
y.fromB.lasso[St.cont] = Z.train[St.cont,] %*% B.cont.lasso
y.fromB.lasso[St.back] = Z.train[St.back,] %*% B.back.lasso

MSE.fromB.lasso = sum((y.train.vec - y.fromB.lasso)**2)

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
  
MSE.GP0nly = sum((y.train.vec - y.fromGPOnly)**2)

surface3d(t.train, x.train, y.train, col='green', alpha = .4)
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.fromGPOnly), nt.train, nx.train), col='orange', alpha = .4)
  
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.fromGPOnly), nt.train, nx.train), col='orange', alpha = .6)
surface3d(t.train, x.train, matrix(as.vector(y.full), nt.train, nx.train), col='blue', alpha = .4)

pdf('Plot_GP_Full_BackMonth.pdf', width = 12, height = 9)
par(mfrow = c(2,1))
plot(y.train[,11], typ = 'l', lty = 2, main = 'Back Month - GP vs Full', xlab = '', ylab = 'Spread')
lines(matrix(as.vector(y.full), nt.train, nx.train)[,11], col = 'red')
lines(matrix(as.vector(y.fromGPOnly), nt.train, nx.train)[,11], col = 'blue')

plot(matrix(as.vector(y.full), nt.train, nx.train)[,11] - y.train[,11], typ = 'l', , col = 'red', main = 'Back Month - GP vs Full Error', xlab = 'Historical Time', ylab = 'Error')
lines(matrix(as.vector(y.fromGPOnly), nt.train, nx.train)[,11] - y.train[,11], col = 'blue')
dev.off()

### MSE Summary
MSE.full
MSE.fromB
MSE.fromB.gibbs
MSE.fromB.lasso
MSE.all
MSE.GP0nly

MSE = c(MSE.full, MSE.fromB, MSE.fromB.gibbs, MSE.fromB.lasso, MSE.all, MSE.GP0nly)
pdf('Plot_ModelError.pdf', width = 12, height = 9)
par(mfrow = c(1,1))
lab = c('Full', 'LS', 'Gibbs', 'Lasso', 'LS - GP', 'GP')
barplot(MSE, , axisnames = TRUE, names.arg = lab, main = 'Model Error')
dev.off()

pdf('Plot_ContCoeffs_LS_Lasso.pdf', width = 12, height = 9)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
plot(B.cont[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta')
points(B.cont.lasso[2:12], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.cont[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
points(B.cont.lasso[13:23], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.cont[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
points(B.cont.lasso[24:34], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.cont[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
points(B.cont.lasso[35:45], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
mtext("Contango Coeffcients - LS and Lasso", outer = TRUE, cex = 1.5)
dev.off()

pdf('Plot_BackCoeffs_LS_Lasso.pdf', width = 12, height = 9)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
plot(B.back[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta', ylim = c(-1.5, 1.4))
points(B.back.lasso[2:12], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.back[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
points(B.back.lasso[13:23], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.back[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
points(B.back.lasso[24:34], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.back[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
points(B.back.lasso[35:45], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
mtext("Backwardation Coeffcients - LS and Lasso", outer = TRUE, cex = 1.5)
dev.off()

pdf('Plot_ContCoeffs_Gibbs.pdf', width = 12, height = 9)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
plot(B.cont.gibbs[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta', ylim = c(-1.5, 1.4))
abline(h = 0, col = 'grey')
plot(B.cont.gibbs[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
plot(B.cont.gibbs[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
plot(B.cont.gibbs[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
mtext("Contango Coeffcients - Gibbs", outer = TRUE, cex = 1.5)
dev.off()

pdf('Plot_BackCoeffs_Gibbs.pdf', width = 12, height = 9)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
plot(B.back.gibbs[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta', ylim = c(-1.5, 1.4))
abline(h = 0, col = 'grey')
plot(B.back.gibbs[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
plot(B.back.gibbs[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
plot(B.back.gibbs[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
abline(h = 0, col = 'grey')
mtext("Contango Coeffcients - Gibbs", outer = TRUE, cex = 1.5)
dev.off()

pdf('Plot_Coeffs_GP_Full.pdf', width = 12, height = 9)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
plot(B.GP[2:12], main = 'Imports', xlab = 'Month', ylab = 'Beta', ylim = c(-1.5, 1.4))
points(B.full[2:12], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.GP[13:23], main = 'Demand', xlab = 'Month', ylab = 'Beta')
points(B.full[13:23], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.GP[24:34], main = 'Stocks', xlab = 'Month', ylab = 'Beta')
points(B.full[24:34], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
plot(B.GP[35:45], main = 'Refining', xlab = 'Month', ylab = 'Beta')
points(B.full[35:45], col = 'red', pch = 3)
abline(h = 0, col = 'grey')
mtext("Coeffcients - GP Only vs Full", outer = TRUE, cex = 1.5)
dev.off()

surface3d(t.train, x.train, y.train, col='green', alpha = .4)
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.fromB), nt.train, nx.train), col='orange', alpha = .4)

surface3d(t.train, x.train, y.train, col='green', alpha = .4)
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.fromB.lasso), nt.train, nx.train), col='blue', alpha = .4)

surface3d(t.train, x.train, y.train, col='green', alpha = .4)
surface3d(t.train, x.train, y.train, col='green', alpha = .4, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.fromB.gibbs), nt.train, nx.train), col='blue', alpha = .4)

surface3d(t.train, x.train, y.train, col='green', alpha = .5)
surface3d(t.train, x.train, y.train, col='green', alpha = .5, front = "line", back = "line")
surface3d(t.train, x.train, matrix(0, nt.train, nx.train), col='blue', alpha = .2, front = "line", back = "line")
surface3d(t.train, x.train, matrix(as.vector(y.all), nt.train, nx.train), col='orange', alpha = .5)

### Next Step is in contango

### Predict Bayes

y.pred = y.test[52,]
z.pred = Z.test[which(t.vec.test == 52),]


# Full
y.pred.full = z.pred %*% B.full
y.pred.B = z.pred %*% B.cont
y.pred.gibbs = z.pred %*% B.cont.gibbs
y.pred.lasso = z.pred %*% B.cont.lasso
y.pred.GP = z.pred %*% B.GP

### Prediction linear plus GP
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
pdf('Plot_Predict.pdf', width = 12, height = 9)
plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction', xlab = 'Month', ylab = 'Spread' )
points(y.test[52,])
lines(y.pred.full, col = 2)
points(y.pred.full, pch = 2, col = 2)
lines(y.pred.B, col = 3)
points(y.pred.B, pch = 3, col = 3)
lines(y.pred.gibbs, col = 4)
points(y.pred.gibbs, pch = 4, col = 4)
lines(y.pred.lasso, col = 5)
points(y.pred.lasso, pch = 5, col = 5)
lines(y.pred.BGP, col = 6)
points(y.pred.BGP, pch = 6, col = 6)
lines(g.star, col = 7)
points(g.star, pch = 7, col = 7)
legend(1, -1, c("Actual", "Full", "LS", "Gibbs", "Lasso", "LS-GP", "GP"), col = c(1,2,3,4,5,6,7), pch = c(1,2,3,4,5,6,7), lwd = c(3,0,0,0,0,0,0))
dev.off()

pdf('Plot_Predict_GPs.pdf', width = 12, height = 9)
plot(y.test[52,], typ = 'l', lwd = 3, ylim = c(-2.2, 0.6), main = 'Prediction - LS/GP vs GP', xlab = 'Month', ylab = 'Spread' )
points(y.test[52,])
lines(y.pred.BGP, col = 6)
points(y.pred.BGP, pch = 6, col = 6)
lines(g.star, col = 7)
points(g.star, pch = 7, col = 7)
legend(1, -1, c("Actual", "LS-GP", "GP"), col = c(1,6,7), pch = c(1,6,7), lwd = c(3,0,0))
dev.off()

prederror.full = sum((y.pred.full - y.test[52,])**2)
prederror.LS = sum((y.pred.B - y.test[52,])**2)
prederror.gibbs = sum((y.pred.gibbs - y.test[52,])**2)
prederror.lasso = sum((y.pred.lasso - y.test[52,])**2)
prederror.BGP = sum((y.pred.BGP - y.test[52,])**2)
prederror.GP= sum((g.star - y.test[52,])**2)

prederror = c(prederror.full,prederror.LS, prederror.gibbs, prederror.lasso, prederror.BGP, prederror.GP)
pdf('Plot_PredictError.pdf', width = 12, height = 9)
par(mfrow = c(1,1))
lab = c('Full', 'LS', 'Gibbs', 'Lasso', 'LS - GP', 'GP')
barplot(prederror, axisnames = TRUE, names.arg = lab, main = 'Prediction Error')
dev.off()

# Diagnostics
start.time = proc.time()
Rprof("profile1.out")
Rprof(NULL)
end.time = proc.time()
total.time = end.time - start.time
total.time
summaryRprof("profile1.out")



