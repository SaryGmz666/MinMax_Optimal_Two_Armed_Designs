####################################################
# Script en R para probar nuevos algoritmos.
####################################################

# Cargar paquetes, bibliotecas y funciones necesarias.
rm(list=ls())

library(gurobi) # El solucionador gurobi se utiliza en la función objetivo. 
source("evaluation.r") # Este script contiene la función objetivo.
source("opt.design.r")
set.seed(10)

#====================================================================
# CARGAR UNA INSTANCIA
#====================================================================
n <- 10 # Número de sujetos.
p <- 6 # Número de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)

# Este objeto contiene cinco matrices de covariables de los pacientes de la muestra.
# Estas matrices fueron generadas al azar e incluyen dos niveles solamente, -1 y +1.
# El algoritmo debe ser ejecutado en cada matriz.
print(H.mat)

# Select one matrix
ii <- 1
H <- H.mat[[ii]]
dim(H) # Imprime el tamaño de la matriz.
View(H)
#====================================================================
#EVALUAR UN DISEÑO DE DOS BRAZOS
#====================================================================

# Con fines ilustrativos, generamos un diseño de dos brazos al azar.
# Los tratamientos del diseño se codifican como -1 o +1.
X <- sample(c(-1,1), size = n, replace = TRUE)
X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.

# La función objetivo es 'obj_fun' en el archivo "evaluation.r".
# Esta función calcula la varianza máxima para todos los pacientes posibles.

# Calcula el criterio de optimalidad G: Varianza máxima.
obj.val <- obj_fun(H, X, max.time = 300, approx = TRUE, print.output = 1)
print(obj.val)

# Calcular el criterio de optimalidad I: Varianza media.
Ioptimality(H, X, approx = TRUE)

#====================================================================
# UTILIZAR UN ALGORITMO TRADICIONAL PARA CONSTRUIR DISEÑOS ÓPTIMOS.
#====================================================================
set.seed(10)
ii <- 1
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]

# Implementamos un algoritmo de intercambio de coordenadas, que es un enfoque de optimización tradicional
# de optimización tradicional en el campo del diseño de experimentos. Se trata de un algoritmo de búsqueda local
# con una estrategia de primera mejora.
# Establecer el número máximo de iteraciones.
max.iter <- 5
# Esta tolerancia es diferente a la de Zhang et al. (2020), que utilizaron 0,0001.
tol <- 0.0000001  #0.0000001 

# Primero, generamos una solución inicial al azar.
set.seed(10)
X <- sample(c(-1,1), size = n*max.iter, replace = TRUE)
X <- matrix(X, ncol = max.iter, nrow = n) # La función objetivo funciona con matrices.

my.time <- proc.time()
for (j in 1:max.iter){
  obj.curr <- obj_fun(H, as.matrix(X[,j]), max.time = 300, print.output = 0)
  npass <- n*20
  obj.pass <- 10*3
  citer <- 1
  while (citer <= npass){
    for (i in 1:n){
      # One move: Flip the sign of a coordinate.
      X[i,j] <- -1*X[i,j] 
      
      # Evaluate the resulting design.
      obj.test <- obj_fun(H, as.matrix(X[,j]), max.time = 300, print.output = 0)
      
      if ( (obj.test - obj.curr) < -1*tol  ){
        # If improvement. Update best objective value
        obj.curr <- obj.test
      } else {
        # Flip back the coordinate.
        X[i,j] <- -1*X[i,j] 
      }
      
    } # end for (i in 1:n)
    
    # Check if we need to go over all coordinates again.
    if ( abs(obj.pass - obj.curr) > tol ){
      obj.pass <- obj.curr
    } else {
      break
    }
    citer <- citer + 1
  } # end while
} # end for (j in 1:max.iter)
print(proc.time() - my.time)

# Evaluate optimal designs as Zhang et al. (2021)
evaluate.designs.X <- est_lp(H, X, Zs,max.time = 300,print.output = 0)
min(evaluate.designs.X[,1])

# n = 100 and p = 6
# 0.07824436 H1
# 0.08774702 H2
# 0.07517031 H3
# 0.1012448 H4
# 0.07965166 H5

# n = 100 and p = 10
# [1] 0.1633114 H1
# [1] 0.1525554
# [1] 0.1548249
# [1] 0.164883
# [1] 0.1350688

# n = 150 and p = 6
# [1] 0.05631823
# [1] 0.05313066
# [1] 0.05585053
# [1] 0.04897471
# [1] 0.05505756