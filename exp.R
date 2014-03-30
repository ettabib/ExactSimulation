A <- function(u){
    return (exp(-u ^ 2 / 2))
}
a <- function(u){
    return (-u * exp(-u ^ 2 / 0.2e1))
}
phi <- function(u){
    return (u ^ 2 * exp(-u ^ 2 / 0.2e1) ^ 2 / 0.2e1 - exp(-u ^ 2 / 0.2e1) / 0.2e1 + u ^ 2 * exp(-u ^ 2 / 0.2e1) / 0.2e1)
}
min.phi=-0.5
