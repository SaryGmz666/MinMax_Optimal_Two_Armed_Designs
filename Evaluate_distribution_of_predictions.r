rm(list=ls())

library(gurobi) # El solucionador gurobi se utiliza en la función objetivo. 
source("evaluation.r") # Este script contiene la función objetivo.
source("opt.design.r")
source("plot_functions.r")
set.seed(10)

n <- 100 # Número de sujetos.
p <- 6 # Número de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
#file.results.instance <- matrix(sample(c(-1,1), n*p, replace = T), ncol = p)
load(file.results.instance)


ii <- 4
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]
X.opt.iter <- opt.design.iter(H, balance=TRUE, max.time = 20) # 500
X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = 300) 
X <- data.frame("Exact" = X.opt.iter, "Approx" = X.opt.approx)

#JMP1 <- read.table("Ioptdesign.txt", sep = ',', header = TRUE)
#JMP2 <- read.table("Doptdesign.txt", sep = ',', header = TRUE)
#X <- data.frame("Exact" = X.opt.iter, "Approx" = X.opt.approx, "Iopt" = JMP1[,1],"Dopt" = JMP2[,1])

fds.plot(X, H, Zs)

#============================================================
# Construcción de diseños I-óptimos.
#============================================================


n <- 100 # Número de sujetos.
p <- 6 # Número de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)


ii <- 3
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]
X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = 300) 


# Implementamos un algoritmo de intercambio de coordenadas, que es un enfoque de optimización tradicional
# de optimización tradicional en el campo del diseño de experimentos. Se trata de un algoritmo de búsqueda local
# con una estrategia de primera mejora.
# Establecer el número máximo de iteraciones.
max.iter <- 100

# Esta tolerancia es diferente a la de Zhang et al. (2020), que utilizaron 0,0001.
tol <- 0.0000001  #0.0000001 

# Primero, generamos una solución inicial al azar.
set.seed(10)
X <- sample(c(-1,1), size = n*max.iter, replace = TRUE)
Y <- matrix(X, ncol = max.iter, nrow = n) # La función objetivo funciona con matrices.

my.time <- proc.time()
for (j in 1:max.iter){
  obj.curr <- Ioptimality(H, as.matrix(Y[,j])) + obj_fun(H, as.matrix(Y[,j]), print.output = 0)
  npass <- n*20
  obj.pass <- 10*3
  citer <- 1
  while (citer <= npass){
    for (i in 1:n){
      # Un movimiento: Invertir el signo de una coordenada.
      Y[i,j] <- -1*Y[i,j] 
      
      # Evalúa el diseño resultante.
      obj.test <- Ioptimality(H, as.matrix(Y[,j])) + obj_fun(H, as.matrix(Y[,j]), print.output = 0)
      
      if ( (obj.test - obj.curr) < -1*tol  ){
        # Si la mejora. Actualizar el mejor valor objetivo
        obj.curr <- obj.test
      } else {
        # Voltea la coordenada.
        Y[i,j] <- -1*Y[i,j] 
      }
      
    } # end for (i in 1:n)
    
    # Comprueba si hay que volver a repasar todas las coordenadas.
    if ( abs(obj.pass - obj.curr) > tol ){
      obj.pass <- obj.curr
    } else {
      break
    }
    citer <- citer + 1
  } # end while
} # end for (j in 1:max.iter)
print(proc.time() - my.time)
evaluate.designs.X <- Ioptimality(H, Y) + obj_fun(H, as.matrix(Y[,j]), print.output = 0)
best.des <- which.min(evaluate.designs.X)


#ALL.DESIGNS <- data.frame("Exact" = X.opt.iter, "Approx" = X.opt.approx, "Iopt" = JMP1[,1],"Dopt" = JMP2[,1], "Alan" = X[,best.des])
ALL.DESIGNS <- data.frame("Approx" = X.opt.approx,  "Alan" = Y[,best.des])
fds.plot(ALL.DESIGNS, H, Zs)


Umat <- matrix(runif(1000*p), ncol = ncol(Zs), nrow = nrow(Zs))
Ymat <- Zs*Umat

fds.plot(ALL.DESIGNS, H, Ymat)

Ioptimality(H, ALL.DESIGNS)


est_lp(H, ALL.DESIGNS, Zs,max.time = 300,print.output = 0)

