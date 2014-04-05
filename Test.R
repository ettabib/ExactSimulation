## Defining Data
T=1
x=100
source("exp.R")
source("sin.R")
source("black.R")
source("Asian.R")
source("Myfunction.R")
source("Euler.R")
Z <- Algorithme.1(PAS=2^-10)
plot(Z$phi.Z-min.phi,type="l")
plot(Z$Z,type="l")
points(y=Z$V,x=Z$U.index)
N=2^10
X <- euler(N=N,alpha=a,X0=100,T=1)
plot(X,type="l")
mean(Z$Z)
mean(X)
