##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title  
##' @param u 
##' @param T  
##' @return 
##' @author Ettabib Mohammad
h <- function(u,T=T){
    return (1/sqrt(2*pi)*exp(-(u-x)^2/(2*T)))
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param n the number of realisations 
##' @return a realisation of the function h
##' @author Ettabib Mohammad
sim.h <- function(n=1){
    return (rnorm(n=n,mean=x,sd=sqrt(T)))
}

##' .. content for \description{} (no empty lines) ..
##' build a 
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @param y 
##' @param N 
##' @param t0 
##' @param t1 
##' @return path as vector
##' @author Ettabib Mohammad
sim.path.brow.t0.t1 <- function(x,y,t0,t1,pas){
    if (t1-t0>pas){
        sigma = sqrt(t1 - t0) / 2
        Y <- rnorm(n = 1,mean = (x + y) / 2,sd = sigma)
        pathl <- sim.path.brow.t0.t1(x = x,y = Y,t0 = t0,t1 = (t0 + t1) / 2,pas = pas)
        pathr <- sim.path.brow.t0.t1(x = Y,y = y,t0=(t0 + t1) / 2,t1 = t1,pas = pas)
        return (c(pathl,Y,pathr))
    }else{
        return (c())
    }
}


## Simulation of a brownian path knowing that x and y
                                        #plot(sim.path.brow.t0.t1(x=0,y=1,t0=0,t1=1,pas=2^-10),type="l")



##' this is the implementation of the first algorithme of Beskos
##' @title 
##' @param X 
##' @param T0 
##' @param T1 
##' @param PAS 
##' @return 
##' @author Ettabib Mohammad
Algorithme.1 <- function(X=x,T0=0,T1=T,PAS=2^-4){
    res <- list()
    res$V <- c()
    res$Z <- c()
    res$U.index <- c()
    res$phi.Z <- c()
    while(TRUE){     
                                        #browser()
        N <- c()
        M <- c()
        U <- c()
        V <- c()
        Z <- c()
        Y <- c()
        Y <- sim.h()
        Z <- sim.path.brow.t0.t1(x = X,y = Y,t0 = T0,t1 = T1,pas = PAS)
        index.min <- which.min(Z)
        t.min <- T0 + index.min * PAS
        M <- max(phi(Z[index.min:length(Z)]) - min.phi)
        N <- rpois(n = 1, lambda = T * M)
        if(N>0){
            U <- runif(n = N, min = 0,max = T)
            V <- runif(n = N, min = 0,max = M)
            i=floor((U-T0)/PAS)+1
            if(length(i)>0){
                if(phi(Z[i])-min.phi<V){
                    res$V <- V
                    res$Z <- Z
                    res$phi.Z <- phi(Z)
                    res$U.index <- i
                    return (res)
                }
            }
        }
    }
}

library(statmod)
dens.inv.invgauss <- function(mu,lambda,u){
    return (sqrt(lambda/(2*pi*u))*exp(lambda/(2*mu))*exp(-lambda/(2*mu^2)*(1/u+u*mu^2)))
}

g <- function(u,lambda){
    return (-lambda/2*exp(-lambda/2*u))
}

rinv.invgauss <- function(n=1,mu,lambda){
                                        # on simule l'inverse de invgauss en prenant g ~ exp(-lambda/2)
    c <- dens.inv.invgauss(mu=mu,lambda=lambda,u=lambda/(mu^2))/g(u=lambda/(mu^2),lambda=lambda)
    while(TRUE){
        y <- rexp(n=1,rate=lambda/2)
        u <- runif(n=1,min=0,max=1)
        if(-c*u*g(y,lambda)<dens.inv.invgauss(mu=mu,lambda=lambda,u=y)){
            return(y)
        }
    }
}

Algorithme.min.2 <- function(a,b,T){
    #browser()  
  c1=(a-b)^2/T
    c2=b^2/(2*T)
    while(TRUE){
        U <- runif(n=1,min=0,max=1)
        Z1 <- (a-sqrt(2*rexp(n=1,rate=1)+a^2))/2
        I1 <- rinvgauss(n=1,mu=sqrt(c1/c2),lambda=2*c1)
        I2 <- rinv.invgauss(n=1,mu=sqrt(c2/c1),2*c2)
        if(U < (1+sqrt(c1/c2))^-1){
            V <- I1
        }else{
            V <- I2
        }
        Z2=T/(1+V)
        res <- c()
        res$min <- Z1
        res$t.min <- Z2
        return (res)
    }
}
## Z <- Algorithme.1()
## plot(Z,type="l")





