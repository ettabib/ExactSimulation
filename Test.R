## Defining Data
source("exp.R")
source("sin.R")
source("black.R")
source("Asian.R")
source("Myfunction.R")
source("Euler.R")

### test Sin
test.loi.euler.sin <- function(Ndisc=1000,Nsimu=100,f=a){
  Z <- c()
  for(i in 1:Nsimu){
    path <- euler(N=Ndisc,alpha=f,X0=X0,T=1)
    Z <- c(Z,path[length(path)])
  }
  return (Z)
}

test.loi.exat.sin <- function(Nsimu=100,Ndisc=1000){
  Z<-c()
  for(i in 1:Nsimu){
    path <- Algorithme.2.bis(N=Nsimu,NT=Ndisc,T=1,x0=X0)$Z
    Z <- c(Z,path[length(path)])
  }
  return (Z)
}

path.euler <- test.loi.euler.sin()
path.exact <- test.loi.exat.sin(Nsimu=10000)
