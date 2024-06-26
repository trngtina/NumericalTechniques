
# Code - Tina -------------------------------------------------------------

library(orthopolynom)
library(dplyr)

# Initialisation -------------------------------------------------------------

# volatilites
s_imp <- 0.25; s_hist <- 0.15;

# rendements
r <- 0.01; mu <- 0.02;

# Parametres du contrat
S0 <- 1; VM0 <- 110; x =0.3; PM0 = 100; R = PM0/VM0
rg <- 0.005

# Parametres de simulations 
# (N = nb simulations primaires, M = nb simulations secondaires)
N <- 3000; M = 1000;

# evolution de l'actif risque en t ayant une volatilite "vol" et drift "d"
# S0 valeur initiale de l'actif
St <- function(S0, d, vol, t, N){
  S0*exp((d-0.5*vol^2)*t+vol*sqrt(t)*rnorm(N))
}

# evolution du fonds VMt a la date t
VMt <- function(St, t, x, VM0, r){
  (1-x)*VM0*exp(r*t) + x*St*VM0
}

BEt <- function(St, t, x, VM0, r){
  R*VM0*( (1-x)*exp(r*t) + x*St )
}

Ot <- function(St, t, x, VM0, r, rg){
  pmax(R*VM0*(exp(rg*t)-(1-x)*exp(r*t)) - R*x*VM0*St,0)
}

# NAV0 -------------------------------------------------------------

# Calcul NAV0_k
N <- 50000;
K <- 50000; # (parametre MC)

## Calculs des NAV0 par MC
S1 <- St(S0, mu, s_hist, 1, N)
S2 <- St(S1, r, s_imp, 1, 1)

VM1 <- VMt(S1, 1, x, VM0, r)
VM2 <- VMt(S2, 2, x, VM0, r)
O2 <- Ot(S2, 2, x, VM0, r, rg)

# NAV0_k
NAV0_k <- VM0 - PM0 - mean(exp(-r*2)*O2)
print(NAV0_k)


# Methode MC pour NAV0
  # Methode rapide

MC_NAV0 <- function(N, K, r, rg){
  NAV0 = c()
  for(k in 1:K){
    S1 <- St(S0, mu, s_hist, 1, N)
    S2 <- St(S1, r, s_imp, 1, 1)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV0_k <- VM0 - PM0 - mean(exp(-r*2)*O2)
    NAV0[k] = NAV0_k
  }
  mean(NAV0)
}

u = sapply(seq(10,500,10),FUN=MC_NAV0, K=3000, r = r, rg= rg)
plot(seq(10,500,10),u,ylim=c(0,10), col = "deepskyblue4", pch='.', main= "Simulation de Monte Carlo",
     xlab = "Nombre de simulations Monte Carlo",
     ylab = "NAV0 avec 3000 simulations")

u = sapply(seq(0,4000,100),FUN=MC_NAV0, N=3000, r= r, rg=rg)
plot(seq(0,4000,100),u,ylim=c(0,10), col = "deepskyblue4", pch='o', main= "Simulation de Monte Carlo",
     xlab = "Nombre de simulations Monte Carlo",
     ylab = "NAV0 avec 3000 simulations")

NAV0 = MC_NAV0(30000,3000,r,rg)
print(NAV0)


# NAV1 -------------------------------------------------------------

NAV1 = c()

S1 <- St(S0, mu, s_hist, 1, N)

for(i in 1:N){
  
  S2 <- St(S1[i], r, s_imp, 1, M)
  VM2 <- VMt(S2, 2, x, VM0, r)
  BE2 <- BEt(S2, 2, x, VM0, r)
  O2 <- Ot(S2, 2, x, VM0, r, rg)
  NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
}

NAV1_act <- NAV1/(1+r)
SCR <- quantile(NAV0 - NAV1_act, .995) #SCR <- NAV0 - quantile(NAV1_act, .005)
SCR

#hist(NAV1) 

plot(S1,NAV1, col = "deepskyblue4", pch='.')

# LSMC -------------------------------------------------------------

model <- lm(NAV1 ~ poly(S1, 3, raw = FALSE)) # Raw = FALSE, les polynomes forment une base orthonormee
summary(model)
points(S1,model$fitted.values, col = "red", pch='.')

test = data.frame(c(0.8,1,1.2,1.4))
colnames(test)=c("S1") # Il faut nommer le vecteur par S1
predictions <- model %>% predict(test)

points(test$S1,predictions, col = "yellow")

NAV1_fitted = model$fitted.values

SCR_poly = NAV0 - quantile(NAV1_fitted, p = .005)/(1+r)
SCR_poly
AIC(model)

fit_polynom <- function(d, vect){
  # NAV1 et S1 sont des variables globales 
  model <- lm(NAV1 ~ poly(S1, d, raw = FALSE))
  vect = data.frame(vect)
  colnames(vect)=c("S1")
  predictions <- model %>% predict(vect)
  predictions
}

fit_polynom(3,c(0.8,1,1.2,1.4))

# AIC en fonction du degre d du modele
# et SCR en fonction du degre d du modele

temp_ = c()
temp__ = c()
for(d in 1:5){
  model <- lm(NAV1 ~ poly(S1, d, raw = FALSE))
  #temp_ = c(temp_,AIC(model)) # AIC
  NAV1_fitted = model$fitted.values
  
  temp_ = c(temp_,mean((NAV1_fitted - NAV1)^2) %>% sqrt()) # RMSE
  
  SCR_poly = NAV0 - quantile(NAV1_fitted,p = .005)/(1+r)
  temp__ = c(temp__,SCR_poly)
  # points(S1,model$fitted.values, col = "red", pch='.')
}

plot(1:5, temp_, col = "deepskyblue4", pch=".", cex= 7, xlim = c(1,5),
     main="RMSE en fonction du degré du polynôme",
     xlab = "Degré du polynôme (base polynomiale canonique)",
     ylab = "RMSE")

plot(1:5, temp__, col = "deepskyblue4", pch=".", cex= 7, xlim = c(1,5), 
     main="SCR en fonction du degré du polynome",
     xlab = "Degre du polynome (base polynomiale  canonique)",
     ylab = "SCR")




# Questions -------------------------------------------------------------

MC_SCR <- function(N, M, r, rg){
  
  # set.seed(4)
  # NAV0 = MC_NAV0(30000,3000,r,rg)
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N)
  
  for(i in 1:N){
    S2 <- St(S1[i], r, s_imp, 1, M)
    VM2 <- VMt(S2, 2, x, VM0, r)
    BE2 <- BEt(S2, 2, x, VM0, r)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
  }
  
  NAV1_act <- NAV1/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
  
}

set.seed(4)
T1<-Sys.time()
MC_SCR(3000,10,r,rg)
T2<-Sys.time()
difftime(T2, T1)

LSMC_SCR <- function(N,M,d,r,rg){
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N) # Methode LSMC pour avoir NAV1 tout de suite pour chaque S1
  
  NAV1 = fit_polynom(d,S1)
  
  NAV1_act <- NAV1/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
}

set.seed(4)
T1<-Sys.time()
LSMC_SCR(3000,30000,3,r,rg)
T2<-Sys.time()
difftime(T2, T1)

LSMC_SCR_dynamique <- function(N,M,d,r,rg){
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N)
  
  for(i in 1:N){
    S2 <- St(S1[i], r, s_imp, 1, M)
    VM2 <- VMt(S2, 2, x, VM0, r)
    BE2 <- BEt(S2, 2, x, VM0, r)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
  }
  
  model <- lm(NAV1 ~ poly(S1, d, raw = FALSE))
  NAV1_fitted = model$fitted.values
  
  NAV1_act <- NAV1_fitted/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
}

LSMC_SCR_dynamique(3000,3,3,r,rg)

LSMC_SCR_dynamique_chebyshev <- function(N,M,d,r,rg){
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N)
  
  for(i in 1:N){
    S2 <- St(S1[i], r, s_imp, 1, M)
    VM2 <- VMt(S2, 2, x, VM0, r)
    BE2 <- BEt(S2, 2, x, VM0, r)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
  }
  
  poly_chebyshev = chebyshev.c.polynomials(n <- d, normalized=TRUE)
  chebyshev = as.matrix(as.data.frame(polynomial.values(polynomials <- poly_chebyshev ,x=scaleX(S1, u=-1, v=1))))[,2:(n+1)] # On ne prend pas la constante et on prend la derniere valeur
  model_chebyshev = lm(NAV1 ~ chebyshev)
  NAV1_chebyshev = model_chebyshev$fitted.values
  
  NAV1_act <- NAV1_chebyshev/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
}

LSMC_SCR_dynamique_hermite <- function(N,M,d,r,rg){
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N)
  
  for(i in 1:N){
    S2 <- St(S1[i], r, s_imp, 1, M)
    VM2 <- VMt(S2, 2, x, VM0, r)
    BE2 <- BEt(S2, 2, x, VM0, r)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
  }
  
  poly_hermite = hermite.h.polynomials(n <- d, normalized=TRUE)
  hermite = as.matrix(as.data.frame(polynomial.values(polynomials <- poly_hermite ,x=scaleX(S1, u=-1, v=1))))[,2:(n+1)] # On ne prend pas la constante et on prend la derniere valeur
  model_hermite = lm(NAV1 ~ hermite)
  NAV1_hermite = model_hermite$fitted.values
  
  NAV1_act <- NAV1_hermite/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
}


# Question 1 -------------------------------------------------------------
# Etudier les estimations precedentes du SCR en fonction de N, M et d (nombre de polynômes).
# N: simulations primaires
# M: simulations secondaires

  # LSMC: Quand d, le degre du polynome augmente, le modele est moins robuste mais la justesse diminue

par(mfrow=c(1,1))

set.seed(1)
LSMC_SCR_N_variable = sapply(seq(50,6000,5), FUN = LSMC_SCR, M = 500, d = 3, r=r, rg=rg)
plot(LSMC_SCR_N_variable, col = "deepskyblue4", ylim = c(0,10), pch='.', cex = 4, type ="l",
     main = "Sensibilite du SCR au nombre de degres du polyn?me",
     ylab = "SCR",
     xlab = "Nombre de simulations")

set.seed(1)
LSMC_SCR_N_variable = sapply(seq(50,6000,5), FUN = LSMC_SCR, M = 500, d = 2, r=r, rg=rg)
points(LSMC_SCR_N_variable, col = "deepskyblue3", ylim = c(0,10), pch='.', cex = 4, type ="l")

T1<-Sys.time()
set.seed(1)
LSMC_SCR_N_variable = sapply(seq(50,6000,5), FUN = LSMC_SCR, M = 500, d = 1, r=r, rg=rg)
T2<-Sys.time()
difftime(T2, T1)
points(LSMC_SCR_N_variable, col = "deepskyblue2", ylim = c(0,10), pch='.', cex = 4, type ="l")
legend(x="bottomright", legend=c("d=3", "d=2", "d=1"), col=c("deepskyblue4","deepskyblue3","deepskyblue2"), lwd=1, cex=0.8, xjust = 1, yjust = 1, title = "Nb de degré(s)")


  # LSMC: Quand N varie
set.seed(1)
LSMC_SCR_N_variable = sapply(seq(50,10000,10), FUN = LSMC_SCR_dynamique, M = 100, d=3,r=r, rg=rg)
plot(seq(50,10000,10),LSMC_SCR_N_variable, col = "deepskyblue4", ylim = c(0,10), pch='.', cex = 4, type ="l")

  # LSMC: Quand M varie
set.seed(1)
LSMC_SCR_M_variable = sapply(seq(50,10000,50), FUN = LSMC_SCR_dynamique, N = 900, d = 3, r=r,rg = rg)
plot(seq(50,10000,50), LSMC_SCR_M_variable, col = "deepskyblue4", ylim = c(0,10), pch='.', cex = 4, type ="l")

# Question 2 -------------------------------------------------------------

# Intervalles de confiance

# d variable 
  # N = 350, M = 150

N = 350
M = 150
IC_LSMC = data.frame('id'=1:600)
IC_LSMC_dyn = data.frame('id'=1:600)

iter <- seq(1,5,1)

for (i in iter) {
  
  set.seed(i)
  IC_SCR = sapply(rep(N,600), FUN = LSMC_SCR, d = i, M=M, r=r,rg=rg)
  IC_LSMC[toString(i)] = IC_SCR
  
  set.seed(i)
  IC_SCR = sapply(rep(N,600), FUN = LSMC_SCR_dynamique, d = i, M=M, r=r,rg=rg)
  IC_LSMC_dyn[toString(i)] = IC_SCR
  print(i)
  
}

boxplot(IC_LSMC[,2:6], col="brown1", ylim=c(4,11),main="Distribution du SCR et degré de polynôme \n LSMC indépendante des simulations secondaires", xlab = "d", ylab = "600 simulations (N = 350, M = 150)")
boxplot(IC_LSMC_dyn[,2:6], col="deepskyblue2", ylim=c(4,11),main="Distribution du SCR et degré de polynôme \n LSMC régression dynamique", xlab = "d", ylab = "600 simulations (N = 350, M = 150)")


# M variable

N = 100
IC_MC = c()
IC_LSMC = c()
IC_LSMC_dyn = c()

IC_MC_sup =c()
IC_MC_inf =c()
IC_LSMC_sup =c()
IC_LSMC_inf =c()
IC_LSMC_dyn_sup =c()
IC_LSMC_dyn_inf =c()
iter <- seq(5,100,2) #seq(5,350,10)
for (i in iter) {
  
  set.seed(i+1000)
  IC_SCR = sapply(rep(N,200), FUN = MC_SCR, M=i, r=r,rg=rg)
  IC_MC_sup = c(IC_MC_sup,t.test(IC_SCR)$conf.int[2])
  IC_MC_inf = c(IC_MC_inf,t.test(IC_SCR)$conf.int[1])
  IC_MC=c(IC_MC,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  set.seed(i+1000)
  IC_SCR = sapply(rep(N,200), FUN = LSMC_SCR, d = 3, M=i, r=r,rg=rg)
  IC_LSMC_sup = c(IC_LSMC_sup,t.test(IC_SCR)$conf.int[2])
  IC_LSMC_inf = c(IC_LSMC_inf,t.test(IC_SCR)$conf.int[1])
  IC_LSMC=c(IC_LSMC,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  set.seed(i+1000)
  IC_SCR = sapply(rep(N,200), FUN = LSMC_SCR_dynamique, d = 3, M=i, r=r,rg=rg)
  IC_LSMC_dyn_sup = c(IC_LSMC_dyn_sup,t.test(IC_SCR)$conf.int[2])
  IC_LSMC_dyn_inf = c(IC_LSMC_dyn_inf,t.test(IC_SCR)$conf.int[1])
  IC_LSMC_dyn=c(IC_LSMC_dyn,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  print(i)
  
}

plot(iter,IC_MC_sup, cex = 4, ylim = c(4,10), type ="l", col = "black",main="Intervalle de confiance à 95% \n (N = 100 et 200 simulations du SCR)", xlab ="M",ylab = "IC à 95%")
points(iter,IC_MC_inf, cex = 4, type ="l", col = "black",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_sup, cex = 4, type ="l", col = "brown1",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_inf, cex = 4, type ="l", col = "brown1",main="Différence entre la borne supérieure et inférieure \n de l'intervalle de confiance à 95%", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
legend(x="bottomright", legend=c("SdS", "LSMC régression indépendante \n des simulations secondaire"), col=c("black","brown1"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

plot(iter,IC_MC_sup, cex = 4, ylim = c(4,10),type ="l", col = "black",main="Intervalle de confiance à 95% \n (N = 100 et 200 simulations du SCR)", xlab ="M",ylab = "IC à 95%")
points(iter,IC_MC_inf, cex = 4, type ="l", col = "black", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_sup, cex = 4, type ="l", col = "deepskyblue3",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_inf, cex = 4, type ="l", col = "deepskyblue3",main="Différence entre la borne supérieure et inférieure \n de l'intervalle de confiance à 95%")
legend(x="bottomright", legend=c("SdS", "LSMC régression dynamique"), col=c("black","deepskyblue3"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

# N variable
M = 100
IC_MC = c()
IC_LSMC = c()
IC_LSMC_dyn = c()

IC_MC_sup =c()
IC_MC_inf =c()
IC_LSMC_sup =c()
IC_LSMC_inf =c()
IC_LSMC_dyn_sup =c()
IC_LSMC_dyn_inf =c()
iter <- seq(5,350,10)
for (i in iter) {
  
  set.seed(i)
  IC_SCR = sapply(rep(M,200), FUN = MC_SCR, N=i, r=r,rg=rg)
  IC_MC_sup = c(IC_MC_sup,t.test(IC_SCR)$conf.int[2])
  IC_MC_inf = c(IC_MC_inf,t.test(IC_SCR)$conf.int[1])
  IC_MC=c(IC_MC,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  set.seed(i)
  IC_SCR = sapply(rep(M,200), FUN = LSMC_SCR, d = 3, N=i, r=r,rg=rg)
  IC_LSMC_sup = c(IC_LSMC_sup,t.test(IC_SCR)$conf.int[2])
  IC_LSMC_inf = c(IC_LSMC_inf,t.test(IC_SCR)$conf.int[1])
  IC_LSMC=c(IC_LSMC,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  set.seed(i)
  IC_SCR = sapply(rep(M,200), FUN = LSMC_SCR_dynamique, d = 3, N=i, r=r,rg=rg)
  IC_LSMC_dyn_sup = c(IC_LSMC_dyn_sup,t.test(IC_SCR)$conf.int[2])
  IC_LSMC_dyn_inf = c(IC_LSMC_dyn_inf,t.test(IC_SCR)$conf.int[1])
  IC_LSMC_dyn=c(IC_LSMC_dyn,t.test(IC_SCR)$conf.int[2]-t.test(IC_SCR)$conf.int[1])
  
  print(i)
  
}

plot(iter,IC_MC_sup, cex = 4, type ="l", col = "black",main="Intervalle de confiance à 95% \n (M = 100 et 200 simulations du SCR)", xlab ="N",ylab = "IC à 95%")
points(iter,IC_MC_inf, cex = 4, type ="l", col = "black", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_sup, cex = 4, type ="l", col = "brown1",xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_inf, cex = 4, type ="l", col = "brown1",main="Différence entre la borne supérieure et inférieure \n de l'intervalle de confiance à 95%", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
legend(x="bottomright", legend=c("SdS", "LSMC régression indépendante \n des simulations secondaire"), col=c("black","brown1"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

plot(iter,IC_MC_sup, cex = 4, type ="l", col = "black",main="Intervalle de confiance à 95% \n (M = 100 et 200 simulations du SCR)", xlab ="N",ylab = "IC à 95%")
points(iter,IC_MC_inf, cex = 4, type ="l", col = "black", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_sup, cex = 4, type ="l", col = "deepskyblue3",xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_inf, cex = 4, type ="l", col = "deepskyblue3",main="Différence entre la borne supérieure et inférieure \n de l'intervalle de confiance à 95%", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
legend(x="bottomright", legend=c("SdS", "LSMC régression dynamique"), col=c("black","deepskyblue3"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

## Application numerique
iter <- seq(5,120,2)
IC_MC_sup =c()
IC_MC_inf =c()
IC_LSMC_dyn_sup =c()
IC_LSMC_dyn_inf =c()
for (i in iter) {
  set.seed(i)
  IC_SCR = sapply(rep(i,800), FUN = MC_SCR, N=100, r=r,rg=rg)
  IC_MC_sup = c(IC_MC_sup,t.test(IC_SCR)$conf.int[2])
  IC_MC_inf = c(IC_MC_inf,t.test(IC_SCR)$conf.int[1])
  
  set.seed(i)
  IC_SCR = sapply(rep(i,800), FUN = LSMC_SCR_dynamique, d = 3, N=100, r=r,rg=rg)
  IC_LSMC_dyn_sup = c(IC_LSMC_dyn_sup,t.test(IC_SCR)$conf.int[2])
  IC_LSMC_dyn_inf = c(IC_LSMC_dyn_inf,t.test(IC_SCR)$conf.int[1])
  
  print(i)
}


plot(iter,IC_MC_sup, cex = 4, ylim =c(3,10),type ="l", col = "black",main="Intervalle de confiance à 95% \n (N = 100 et 800 simulations du SCR)", xlab ="M",ylab = "IC à 95%")
points(iter,IC_MC_inf, cex = 4, type ="l", col = "black", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_sup, cex = 4, type ="l", col = "deepskyblue3",xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
points(iter,IC_LSMC_dyn_inf, cex = 4, type ="l", col = "deepskyblue3",main="Différence entre la borne supérieure et inférieure \n de l'intervalle de confiance à 95%", xlab ="N",ylab = "Différence entre la borne supérieure et inférieure")
legend(x="bottomright", legend=c("SdS", "LSMC régression dynamique"), col=c("black","deepskyblue3"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

set.seed(4)
t_deb <- Sys.time()
IC_SCR = sapply(rep(120,800), FUN = MC_SCR, N=100, r=r,rg=rg)
t_fin <- Sys.time()
t_fin-t_deb

set.seed(4)
t_deb <- Sys.time()
IC_SCR = sapply(rep(30,800), FUN = LSMC_SCR_dynamique, d = 3, N=100, r=r,rg=rg)
t_fin <- Sys.time()
t_fin-t_deb

t_deb <- Sys.time()
set.seed(4)
IC_SCR = sapply(rep(30,800), FUN = LSMC_SCR, d = 3, N=100, r=r,rg=rg)
t_fin <- Sys.time()
t_fin-t_deb 
# Bien trop long comme temps de calcul, possibilite de l'optimiser eventuellement ...

# Temps de calcul (suite)
  # M variable 
set.seed(1)
duree_tab = c()
iter <- seq(1000, 10000, 50)
N <- 1000
M <- 1000

for (i in iter) {
  t_deb <- Sys.time()
  MC_SCR(N,i,r,rg)
  t_fin <- Sys.time()
  duree_tab <- c(duree_tab, t_fin-t_deb)
  print(i)
}
plot(iter,duree_tab,ylim=c(0,0.7), cex = 4, type ="l", col = "black", xlab = "M", ylab = "Temps de calcul (s)", main="Temps de calcul en fonction de M \n N = 1000")

set.seed(1)
duree_tab = c()
iter <- seq(1000, 10000, 50)
N <- 1000
M <- 1000

for (i in iter) {
  
  t_deb <- Sys.time()
  LSMC_SCR_dynamique(N,i,3,r,rg)
  t_fin <- Sys.time()
  duree_tab <- c(duree_tab, t_fin-t_deb)
  
  print(i)
}

points(iter,duree_tab,ylim=c(0,0.7),col="deepskyblue3", cex = 4, type ="l")
legend(x="bottomright", legend=c("SdS", "LSMC régression dynamique"), col=c("black","deepskyblue3"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")

# N variable 
set.seed(1)
duree_tab = c()
iter <- seq(1000, 10000, 50)
N <- 1000
M <- 1000

for (i in iter) {
  t_deb <- Sys.time()
  MC_SCR(i,M,r,rg)
  t_fin <- Sys.time()
  duree_tab <- c(duree_tab, t_fin-t_deb)
  print(i)
}
plot(iter,duree_tab,ylim=c(0,0.7), cex = 4, type ="l", col = "black", xlab = "N", ylab = "Temps de calcul (s)", main="Temps de calcul en fonction de N \n M = 1000")

set.seed(1)
duree_tab = c()
iter <- seq(1000, 10000, 50)
N <- 1000
M <- 1000

for (i in iter) {
  
  t_deb <- Sys.time()
  LSMC_SCR_dynamique(i,M,3,r,rg)
  t_fin <- Sys.time()
  duree_tab <- c(duree_tab, t_fin-t_deb)
  
  print(i)
}

points(iter,duree_tab,ylim=c(0,0.7),col="deepskyblue3", cex = 4, type ="l")
legend(x="bottomright", legend=c("SdS", "LSMC régression dynamique"), col=c("black","deepskyblue3"), lwd=1, cex=0.8, xjust = 1, yjust = 1, 
       title = "Méthode")


# d variable 
set.seed(1)
duree_tab = c()
iter <- seq(1,10,1)
N <- 10000
M <- 10000

for (i in iter) {
  
  t_deb <- Sys.time()
  LSMC_SCR_dynamique(N,M,i,r,rg)
  t_fin <- Sys.time()
  duree_tab <- c(duree_tab, t_fin-t_deb)
  
  print(i)
}

plot(iter,duree_tab,col="deepskyblue3", cex = 4, type ="l", xlab = "d", ylab = "Temps de calcul (s)", main="Temps de calcul en fonction de d")

# Question 3 -------------------------------------------------------------
# En se basant sur la méthode LSMC, calculer la sensibilite du SCR par rapport a rg et r

LSMC_SCR_dynamique_chebyshev <- function(N,M,d,r,rg){
  
  # set.seed(4)
  NAV0 = MC_NAV0(30000,3000,r,rg)
  
  NAV1 = c()
  S1 <- St(S0, mu, s_hist, 1, N)
  
  for(i in 1:N){
    S2 <- St(S1[i], r, s_imp, 1, M)
    VM2 <- VMt(S2, 2, x, VM0, r)
    BE2 <- BEt(S2, 2, x, VM0, r)
    O2 <- Ot(S2, 2, x, VM0, r, rg)
    NAV1[i] <- exp(-r*1)*mean(VM2-BE2-O2)
  }
  
  poly_chebyshev = chebyshev.c.polynomials(n <- d, normalized=TRUE)
  chebyshev = as.matrix(as.data.frame(polynomial.values(polynomials <- poly_chebyshev ,x=scaleX(S1, u=-1, v=1))))[,2:(n+1)] # On ne prend pas la constante et on prend la derniere valeur
  model_chebyshev = lm(NAV1 ~ chebyshev)
  NAV1_chebyshev = model_chebyshev$fitted.values
  
  NAV1_act <- NAV1_chebyshev/(1+r)
  SCR <- quantile(NAV0 - NAV1_act, .995)
  SCR
}

r_sensibilite = c()
i=4
set.seed(i)
a = LSMC_SCR_dynamique_chebyshev(6000,200,3,r+0.0001*r,rg)
set.seed(i)
b = LSMC_SCR_dynamique_chebyshev(6000,200,3,r,rg)
r_sensibilite = (a-b)/(0.0001*r)
r_sensibilite
# sensibilite par rapport a r negative, de -80
# quand r, le taux sans risque augmente, le SCR diminue
# de plus, quand r augmente d'une unite, SCR diminue de 80 unites


rg_sensibilite = c()
i=3
set.seed(i)
a = LSMC_SCR_dynamique_chebyshev(6000,200,3,r,rg+0.0001*rg)
set.seed(i)
b = LSMC_SCR_dynamique_chebyshev(6000,200,3,r,rg)
rg_sensibilite = (a-b)/(0.0001*rg)

# sensibilite par rapport a rg positive, de 90
# quand rg, le taux technique augmente, le SCR augmente aussi
# de plus, quand rg augmente d'une unite, SCR augmente de 90 unites