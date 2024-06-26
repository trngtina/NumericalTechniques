### Parametres
sigma = 0.2 # Coef de volatilite
T = 1 # Maturite
K = 800 # Strike
r = 0.05 # Taux sans risque
N = 300 # Nombre de pas | Temps
delta_tau = T/N # Pas de discretisation | Temps
M = 300 # Nombre de pas | Espace
y_max = log(2000) # Valeur "max" du sous-jacent en T, X_T = 2000
y_min = log(10) # Valeur "min" du sous-jacent en T, X_T = 10
delta_y = (y_max - y_min)/(M+1) # Pas de discretisation | Espace


#v0 = -1:M+1
#v0 = y_min+v0*delta_y
#y = exp(v0)
#v0 = K-exp(v0)
#v0[v0<=0] =0
y = seq(y_min, y_max, by=delta_y)
V0 = pmax(0, K-exp(seq(y_min, y_max, by=delta_y)))

a1 = (r-(sigma^2)/2)*(delta_tau/(delta_y*2))+((sigma**2)*(delta_tau))/(2*(delta_y**2))
a2 = 1-r*delta_tau-((sigma**2)*(delta_tau))/(delta_y**2)
a3 = -(r-(sigma^2)/2)*(delta_tau/(delta_y*2))+((sigma**2)*(delta_tau))/(2*(delta_y**2))
A = matrix(0, nrow = M+2, ncol =M+2)
diag(A) = a2
A[row(A)-col(A) == 1] = a1
A[col(A)-row(A) == 1] = a3

for (i in 1:300){
  V0 = A%*%V0
}

plot(exp(y),V0)
