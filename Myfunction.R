

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
sim.path.brow.t0.t1 <- function(X,Y,t0,t1,pas){
    if (t1-t0>pas){
        sigma = sqrt(t1 - t0) / 2
        Y <- rnorm(n = 1,mean = (X + Y) / 2,sd = sigma)
        pathl <- sim.path.brow.t0.t1(X = X,Y = Y,t0 = t0,t1 = (t0 + t1) / 2,pas = pas)
        pathr <- sim.path.brow.t0.t1(X = Y,Y = Y,t0=(t0 + t1) / 2,t1 = t1,pas = pas)
        return (c(pathl,Y,pathr))
    }else{
        return (c())
    }
}

sim.brow.x.y <- function(t,x,y,t0,t1){
  m <- x * (t1-t) / (t1-t0) + y * (t-t0) / (t1-t0)
  v <- (t-t0) * (t1-t) / (t1-t0)
  return ( rnorm(n=1,mean=m,sd=sqrt(v)))
}

sim.brow.C <- function(x,y,Time){
  Z <- c()
  time <- sort(unique(Time))
  Z[length(time)] = y; Z[1] = x;
  for(i in 2:(length(time)-1)){
    Z[i] <- sim.brow.x.y(t=time[i],x=Z[i-1],y=y,t0=time[i-1],t1=time[length(time)])
  }
  return (Z)
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

g <- function(u,lambda) {
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

Algorithme.min.2 <- function(a,T,X0=X0){
    #browser()  
        a <- a - X0
        U <- runif(n=1,min=0,max=1)
        Z1 <- (a-sqrt(2*T*rexp(n=1,rate=1)+a^2))/2
        b <- Z1
        c1=(a-b)^2/T
        c2=b^2/(2*T)
        I1 <- rinvgauss(n=1,mu=sqrt(c1/c2),lambda=2*c1)
        I2 <- rinv.invgauss(n=1,mu=sqrt(c2/c1),lambda=2*c2)
        if(U < (1+sqrt(c1/c2))^-1){
            V <- I1
        }else{
            V <- I2
        }
        Z2=T/(1+V)
        res <- c()
        res$min <- Z1 + X0
        res$t.min <- Z2
        return (res)
}

R <- function(t,delta){
  W.br.1 <- rnorm(n=1,mean=0,sd=sqrt(t))-t*rnorm(n=1,mean=0,sd=1)
  W.br.2 <- rnorm(n=1,mean=0,sd=sqrt(t))-t*rnorm(n=1,mean=0,sd=1)
  W.br.3 <- rnorm(n=1,mean=0,sd=sqrt(t))-t*rnorm(n=1,mean=0,sd=1)
  return (sqrt((delta*t*W.br.1)^2+(W.br.2)^2+(W.br.3)^2))
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param mT the minimum
##' @param theta the time at minimum
##' @param WT the final point of the brownian motion
##' @param s the times of discretization
##' @return 
##' @author Mohammad
decompose <- function(mT,theta,WT,s,W0){
  S <- sort(s)
  Z <- c()
  for(i in 1:length(S)){
    delta.1 <- -mT/sqrt(theta)
    delta.2 <- (WT-mT)/sqrt(T-theta)
    if(S[i]<= theta){
      Z[i] <- sqrt(theta)*R(delta=delta.1,t=(theta-S[i])/theta)+mT+W0
    }else{
      Z[i] <- sqrt(T-theta)*R(delta=delta.2,t=(S[i]-theta)/(T-theta))+mT+W0
    }
  }
  return (Z)
}

Algorithme.2 <-function(T0=0,T=T,PAS=2^-4){
  while(TRUE){
  Z <- c()
  # draw the ending point ZT of the process Z with respect to the density h
  Z.T <- sim.h()
  # simulate th minimum m of the process Z giving Z.T
  m <- Algorithme.min.2(a=Z.T,T=T,X0=X0)
  #Fix an upper bound M(m)
  M.m <- Max.phi(m$min)-min.phi
  #Draw N with poisson distribution T*M
  N <- rpois(n=1,lambda=T*M.m)
  U <- runif(n=N,min=0,max=T)
  V <- runif(n=N,min=0,max=M.m)
  # Fill in the path of Z at the remaining times U
  i=floor((U-T0)/PAS)+1
  Z <- decompose(mT=m$min,theta=m$t.min,WT=Z.T,s=U,W0=X0)
  # evaluate the number of points under M.m
  if(prod(phi(Z)-min.phi<V)==1){
    return (Z)
   }
  }
}

Algorithme.1.bis <- function(N=100,T=1,x0=X0){
  while(TRUE){
    N <- rpois(n=1,lambda=T*Max.phi(0))
    U <- runif(n=N,min=0,max=T)
    V <- runif(n=N,min=0,max=Max.phi(0))
    Z.T <- sim.h()
    
    if(prod(phi(Z)-min.phi<V)==1){
      return (Z)
    }
  }
}

Algorithme.2.bis <- function(N=100,NT=100,T=1,x0=X0){
  Traj <- c()
  res <- c()
  Max=0
  for(i in 1:NT){
    Traj <- sim.path.brow.t0.t1(X=x0,Y=Z.T,t0=0,t1=T,pas=T/N)
    Max <- Max + max(Traj)
  }
  M <- Max-min.phi
  while(TRUE){
    N <- rpois(n=1,lambda=T*M)
    U <- runif(n=N,min=0,max=T)
    V <- runif(n=N,min=0,max=M)
    Z.T <- sim.h()
    Z <- sim.brow.C(x=X0,y=Z.T,Time=U)
    if(prod(phi(Z)-min.phi<V)==1){
      res$Z <- Z
      res$U <- U
      res$V <- V
      res$Max <- M
      return (res)
    }
  }
}


### Calcul de h
Rlambert <- function(x,N){
  w0=1
  N=0
  for(i in 1:N){
    w0=w0-(w0*exp(w0)-x)/((1+w0)*exp(w0))
  } 
  return (w0)
}
#Calcul de C
u.etoile=(gamma * T + Rlambert( T * beta * S0 * exp(-sigma * X0  - gamma * T ))) / sigma + X0


h <- function(u,t=T,x0=X0){
  P0*exp(A(u)-(u-x0)^2/(2*t))
}

f.A <- function(z,u=u.etoile,T){
  return (1/sqrt(2*pi*T)*exp(-(z-u)^2/T))
}

g.A <- function(z){
  return ((beta*S0)/sigma*(1-exp(sigma*z)-sigma*z*exp(-sigma*u.etoile)))
}

C=f(0,u.etoile,T) * exp(g.A(0) - g.A(u.etoile)) / h(0)

#algorithme de rejet
sim.h <- function(t=T,x0=X0){
  X=0
  repeat{
    U <- runif(1)
    Y <- rnorm(1,mean=u.etoile,sd=t)
    if (C*h(Y)/f.A(Y,u.etoile,t)>U) { X=Y;
                                    break()}
  } 
  return (X)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title  
##' @param u 
##' @param T  
##' @return 
##' @author Ettabib Mohammad
## h <- function(u,T=T){
##   return (A(u)-1/sqrt(2*pi)*exp(-(u-x)^2/(2*T)))
## }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param n the number of realisations 
##' @return a realisation of the function h
##' @author Ettabib Mohammad
## sim.h <- function(n=1){
##     return (rnorm(n=n,mean=x,sd=sqrt(T)))
## }

