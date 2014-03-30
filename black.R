r <- 0.05
sigma <- 0.3
delta <- 0
T <- 1
alpha <- 0.6
beta <- 0.4
K <- 100
S0 <- 100
X0 <- S0/sigma
A <- function(u){
    return (1/2*u^2*(r-delta))
}
a <- function(u){
    return (u*(r-delta))
}
phi <- function(u){
    return (1 / 2 * u ^ 2 * (r - delta) ^ 2  + 1 / 2 * (r  - delta))
}
min.phi=1 / 2 * (r - delta)
