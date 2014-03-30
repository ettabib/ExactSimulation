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
plot(sim.path.brow.t0.t1(x=0,y=1,t0=0,t1=1,pas=2^-10),type="l")



##' this is the implementation of the first algorithme of Beskos
##' @title 
##' @param X 
##' @param T0 
##' @param T1 
##' @param PAS 
##' @return 
##' @author Ettabib Mohammad
Algorithme.1 <- function(X=x,T0=0,T1=T,PAS=2^-4){
    while(TRUE){
        ##browser()
        Z <- c()
        Y <- c()
        Y <- sim.h()
        Z <- sim.path.brow.t0.t1(x = X,y = Y,t0 = T0,t1 = T1,pas = PAS)
        index.min <- which.min(Z)
        m <- T0+index.min*PAS
        M <- max(phi(Z[index.min:length(Z)])-min.phi)
        N <- rpois(n=1,lambda=T*M)
        U <- runif(n=N,min=0,max=T)
        V <- runif(n=N,min=0,max=M)
        i=floor((U-T0)/PAS)+1
        if(N>0){
            if(phi(Z[i])-min.phi<V){
                return (Z)
            }
        }
    }
}

## Z <- Algorithme.1()
## plot(Z,type="l")





