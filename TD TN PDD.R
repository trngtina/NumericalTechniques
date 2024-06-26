# Methode naive

T <- 1
N <- 200
delta_T <- T/N
S0 <- 100
sigma <- 0.45
r <- 0.015
Sa <- 120


S_next <- function(S) {
  X <- rnorm(1)
  return(S*(r*delta_T + sigma*sqrt(delta_T)*X + 1)) 
}

l = 1:200

S = S0
for (i in 1:N){
  l[i] = S
  S = S_next(S)
}

PDD <- exp(-r*1)*(Sa-l[N])*(l[N]<=0.8*Sa)*(max(l[N/2:N])<=Sa)

-------
sim = 150
PDD = 0
for (j in 1:sim){
  
  l = 1:200
  
  S = S0
  for (i in 1:N){
    l[i] = S
    S = S_next(S)
  }
  
  PDD <- PDD + exp(-r*1)*(Sa-l[N])*(l[N]<=0.8*Sa)*(max(l[N/2:N])<=Sa)
  
}
PDD_Monte <- PDD/sim