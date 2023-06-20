rm(list=ls())

library(MOEADr)
library(Matrix)
library(gurobi) # El solucionador gurobi se utiliza en la funci?n objetivo. 
source("evaluation.r") # Este script contiene la funci?n objetivo.
source("opt.design.r")
source("plot_functions.r")


################ Creacion de la instancia ################
n <- 100 # N?mero de sujetos.
p <- 6 # N?mero de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)

################## Inicio y busqueda de la solucion optima  ###################

ii <- 1
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]
#X.opt.iter <- opt.design.iter(H, balance=TRUE, max.time = 20) # Exacto
#t1 <- proc.time()
#X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = 300) #Aproximado
#proc.time()-t1
# Calcula el criterio de optimalidad G: Varianza m?xima y guardamos
#obj.val.opt <- obj_fun(H, X.opt.approx, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
#obj.val.opt

# Calcular el criterio de optimalidad I: Varianza media y guardamos
#Iopt <- Ioptimality(H, X.opt.approx, approx = TRUE)
#Iopt


################################################### Inicio de la metaheuristica ###################################################

#Control de variables
# Cardinalidad de la poblacion elite.
N_new         <- 1000    
# Numero de soluciones a tomar en cuenta para cada solucion en el frente de Pareto.
N_elite       <- 50      
# Numero maximo de iteraciones
iter          <- 40      
mostrar_grafica <- TRUE # Mostrar la grafica de evolucion de mejores soluciones.


###### Alta de Matrices ######
#Mariz que contiene las X, los valores de las funciones obj y si pertenece a las soluciones no dominadas
X.pob.E           <- matrix(0, ncol = N_new, nrow = (n+2))

# Alta de matrices y valores constantes que se utilizan en las funciones objetivo.
p <- ncol(H)
n <- nrow(H)
R <- t(H)%*%H # H'H
Rinv <- solve(R) # (H'H)^{-1}
M <- moments.matrix(p-1, n.bin = p-1)
# NOTA: Necesitamos aniadir el factor 2^{p-1} a la matrix M.
# Funcionalidad extra:
modelo.gurobi <- generate_gurobi_model(p, n.bin = p-1)

############################## Creacion de la poblacion inicial ##############################
set.seed(5739584) # EPCEA= 5739584, OLD = 4338663
t <- proc.time()
X <- sample(c(-1,1), size = n*N_new, replace = TRUE)
X.pob.E[1:n,] <- matrix(X, ncol = N_new)
X.pob.E[n+1,] <- Goptimality_fast(H, X.pob.E[1:n,], p, n, R, Rinv, approx = TRUE, modelo.gurobi) # Solucion del nivel inferior
X.pob.E[n+2,] <- Ioptimality_fast(H, X.pob.E[1:n,], p, n, R, Rinv, M, approx = TRUE)

############################## Creacion de la poblacion elite ##############################

# Encontrar las soluciones del frente de pareto
sel_nondominated <- find_nondominated_points( t(X.pob.E[c(n+1,n+2),]))
X.Sol <- X.pob.E[,sel_nondominated]  
# Remover soluciones repetidas.
X.Sol <- X.Sol[,!duplicated(t(X.Sol[c(n+1,n+2),]))]
if(mostrar_grafica){
  plot(X.pob.E[n+2,], X.pob.E[n+1,], main = "Frente de Pareto Inicial",
       xlab = "I-optimalidad", ylab = "G-optimalidad")
  points(X.Sol[n+2,], X.Sol[n+1,], pch = 19, col = 'red')
}

#Crear la poblacion elite para el frente de pareto.------------
X.Elite <- list()
n_sol_pareto <- ncol(X.Sol)
for (i in 1:n_sol_pareto){
  # Para cada solucion en frente de pareto, encontrar las N_elite soluciones mas 
  # cercanas.
  dist_sol_to_i <- colSums((X.pob.E[c(n+1, n+2),] - X.pob.E[c(n+1, n+2),i])^2) 
  aux.i <- order(dist_sol_to_i, decreasing = FALSE)
  # Guardar las N_elite soluciones mas cercanas a cada solucion en frente de pareto.
  X.Elite[[i]] <- X.pob.E[, aux.i[1:N_elite]] 
}

############ Creacion de nuevas soluciones con respecto a la frecuencia calculada ############
for( y in 1:iter ){
  cat("Iteracion", y, "\n")
  X.pob.E           <- matrix(0, ncol = N_new, nrow = (n+2))
  # La siguiente funcion identifica las entradas que son iguales a 1 en la 
  # matriz de soluciones X.E[1:n,]. Despues, calcula la proporcion de 1s 
  # por columna en esta matrix.
  Fre.Par.Sol <- matrix(0, ncol = n_sol_pareto, nrow = n)
  for (i in 1:n_sol_pareto){
    Fre.Par.Sol[,i] <- rowMeans(X.Elite[[i]][1:n,] == 1) # De acuerdo a G-optimalidad
  }

  # Creacion de una nueva poblacion, con respecto a las frecuencias.-----------
  
  #Se guardan las soluciones en el frente de pareto actual.
  X.pob.E[,1:n_sol_pareto] <- X.Sol
  
  # Generar "N_new - N_mejores" soluciones utilizando las probabilidades en Fre.
  # La funcion sapply es basicamente un for loop de 1 a (N_new-N_mejores).
  # En cada iteracion ejecuta el comando 2*(prob >= runif(n))-1, que crea
  # un vector de 0 y 1 de acuerdo a las probabilidades en prob = Fre. Despues
  # transforma ese vector a -1s y 1s: 2*x - 1 con x = 0 o 1.
  
  # Numero de soluciones a generar en la vecindad de cada solucion en el frente
  # de pareto.
  nsols_por_pareto_sol <- ceiling((N_new-n_sol_pareto)/n_sol_pareto)
  
  for (i in 1:n_sol_pareto){
    aux.vec <- (nsols_por_pareto_sol*(i-1) + 1):(i*nsols_por_pareto_sol) + n_sol_pareto
    max.aux.vec <- max(aux.vec)
    if(max.aux.vec > N_new){aux.vec <- head(aux.vec, -(max.aux.vec - N_new)); 
                            nsols_por_pareto_sol <- length(aux.vec)}
    X.pob.E[1:n,aux.vec] <- sapply(1:nsols_por_pareto_sol,
                                               FUN = function(x, prob, n){2*(prob >= runif(n))-1}, 
                                               prob = Fre.Par.Sol[,i], n = n)
  }
  
  
  # Evaluar la nueva poblacion.-------------------------------------------------
  # Ahora, evaluar las "N_new - N_mejores" soluciones.
  # G-optimalidad.
  X.pob.E[n+1,(n_sol_pareto+1):N_new] <- Goptimality_fast(H, X.pob.E[1:n,(n_sol_pareto+1):N_new], 
                                                       p, n, R, Rinv, approx = TRUE, 
                                                       modelo.gurobi)
  # I-optimalidad.
  X.pob.E[n+2,(n_sol_pareto+1):N_new] <- Ioptimality_fast(H, X.pob.E[1:n,(n_sol_pareto+1):N_new], 
                                                       p, n, R, Rinv, M, approx = TRUE)
  
  # Encontrar las soluciones en el frente de Pareto.----------------------------
  sel_nondominated <- find_nondominated_points( t(X.pob.E[c(n+1,n+2),]))
  X.Sol <- X.pob.E[,sel_nondominated]
  if (sum(sel_nondominated) == 1){ # Si solo hay una solucion,
    X.Sol <- matrix(X.Sol, ncol = 1) # Convertir a matriz
  } else { # Si hay mas de una solucion,
    X.Sol <- X.Sol[,!duplicated(t(X.Sol[c(n+1,n+2),]))] # Remover duplicados  
  }
  
  if(mostrar_grafica){
    plot(X.pob.E[n+2,], X.pob.E[n+1,], main = paste("Frente de Pareto, Iteracion", y),
         xlab = "I-optimalidad", ylab = "G-optimalidad")
    points(X.Sol[n+2,], X.Sol[n+1,], pch = 19, col = 'red')
  }
  
  #Actualizacion de las poblacion elite.----------------------------------------
  X.Elite <- list()
  n_sol_pareto <- ncol(X.Sol)
  for (i in 1:n_sol_pareto){
    # Para cada solucion en frente de pareto, encontrar las N_elite soluciones mas 
    # cercanas.
    dist_sol_to_i <- colSums((X.pob.E[c(n+1, n+2),] - X.Sol[c(n+1, n+2),i])^2) 
    aux.i <- order(dist_sol_to_i, decreasing = FALSE)
    # Guardar las N_elite soluciones mas cercanas a cada solucion en frente de pareto.
    X.Elite[[i]] <- X.pob.E[, aux.i[1:N_elite]] 
  }
  
}
save.file <- paste("PXE_nObs_", n, "_nCovar_", p, "_Inst_", ii, ".RData", sep = '')
pxe.time <- proc.time()-t
save(X.Sol, pxe.time, file = save.file)


########## Ya obtenida las soluciones seleccionar la mejor ##########

#Soluciones no dominadas
#indicador <- find_nondominated_points(t(X.Sol[c(n + 1,n + 2),]))
#X.Sol[c(n + 1,n + 2),indicador]
#

