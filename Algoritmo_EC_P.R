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
for(ii in 3:5){
#ii <- 2
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]
#X.opt.iter <- opt.design.iter(H, balance=TRUE, max.time = 20) # Exacto
t1 <- proc.time()
X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = 300) #Aproximado
proc.time()-t1
# Calcula el criterio de optimalidad G: Varianza m?xima y guardamos
obj.val.opt <- obj_fun(H, X.opt.approx, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
obj.val.opt

# Calcular el criterio de optimalidad I: Varianza media y guardamos
Iopt <- Ioptimality(H, X.opt.approx, approx = TRUE)
Iopt

################################################### Inicio de la metaheuristica ###################################################
#Control de variables
N_pob         <- 500     #Cardinalidad de la poblacion
N_elite       <- 50      #Cardinalidad de la poblacion elite
N_new         <- 400     #Cardinalidad de las creaciones nuevas
iter          <- 40      #Numero de iteraciones
#objetivo.ref     <- 2    #1 es para entropia cruzada para G-opt, 2 es para I-opt
N_mejores <- 10  # Numero de soluciones a guardar de la poblacion elite.
mostrar_grafica <- TRUE # Mostrar la grafica de evolucion de mejores soluciones.
muestra <- 5
for(objetivo.ref in 1:2){

tiempo <- list()
X.Soluciones <- list()

for(muestra in 1:muestra){
  cat("(Muestra", muestra, ")\n")
  t <- proc.time()
  ###### Alta de Matrices ######
  #Mariz que contiene las X, los valores de las funciones obj y si pertenece a las soluciones no dominadas
  X.pob           <- matrix(0, ncol = N_pob, nrow = (n+2))
  
  #Creacion de las poblaciones elite y matrices solucion
  X.E           <- matrix(0, ncol = N_elite, nrow = (n+2))
  X.Sol           <- matrix(0, ncol = (iter + 1), nrow = (n+2))
  # Nota: No creo que necesitemos los nombres de las columnas y filas.
  # En cualquier caso, los podemos agregar al final del codigo.
  
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
  
  X <- sample(c(-1,1), size = n*N_pob, replace = TRUE)
  X.pob[1:n,] <- matrix(X, ncol = N_pob)
  X.pob[n+1,] <- Goptimality_fast(H, X.pob[1:n,], p, n, R, Rinv, approx = TRUE, modelo.gurobi) # Solucion del nivel inferior
  X.pob[n+2,] <- Ioptimality_fast(H, X.pob[1:n,], p, n, R, Rinv, M, approx = TRUE)

  ############################## Creacion de la poblacion elite ##############################
  cont <- 1
  #orden de creciente para G optimalidad (objetivo.ref=1) o
  # I optimalidad (objetivo.ref=2)
  aux <- order(X.pob[n+objetivo.ref,], decreasing = FALSE)
  
  #Crear la poblacion elite
  X.Sol[,cont] <- X.pob[, aux[1]] # Guardar la mejor solucion.
  X.E <- X.pob[, aux[1:N_elite]] # Guardar las mejores N_elite soluciones.
  
  if(mostrar_grafica){
    plot(cont, X.Sol[n+objetivo.ref,cont], main = "Evolucion de soluciones",
         xlab = "Iteracion", ylab = "Valor Objetivo (Individual)", 
         xlim = c(1, iter+1))
  }
  
  cont  <- cont + 1
  ############ Creacion de nuevas soluciones con respecto a la frecuencia calculada ############
  for( y in 1:iter ){
    cat("Iteracion", y, "\n")
    X.pob.E           <- matrix(0, ncol = N_new, nrow = (n+2))
    # La siguiente funcion identifica las entradas que son iguales a 1 en la 
    # matriz de soluciones X.E[1:n,]. Despues, calcula la proporcion de 1s 
    # por columna en esta matrix.
    Fre <- rowMeans(X.E[1:n,] == 1)
    
    # Creacion de una nueva poblacion, con respecto a las frecuencias.-----------
    
    #Se guardan las primeras N_mejores=10 soluciones de X.E
    X.pob.E[,1:N_mejores] <- X.E[,1:N_mejores]
    # Generar "N_new - N_mejores" soluciones utilizando las probabilidades en Fre.
    # La funcion sapply es basicamente un for loop de 1 a (N_new-N_mejores).
    # En cada iteracion ejecuta el comando 2*(prob >= runif(n))-1, que crea
    # un vector de 0 y 1 de acuerdo a las probabilidades en prob = Fre. Despues
    # transforma ese vector a -1s y 1s: 2*x - 1 con x = 0 o 1.
    X.pob.E[1:n,(N_mejores+1):N_new] <- sapply(1:(N_new-N_mejores),
                                            FUN = function(x, prob, n){2*(prob >= runif(n))-1}, 
                                            prob = Fre, n = n)
    # Ahora, evaluar las "N_new - N_mejores" soluciones.
    # G-optimalidad.
    X.pob.E[n+1,(N_mejores+1):N_new] <- Goptimality_fast(H, X.pob.E[1:n,(N_mejores+1):N_new], 
                                                         p, n, R, Rinv, approx = TRUE, 
                                                         modelo.gurobi)
    # I-optimalidad.
    X.pob.E[n+2,(N_mejores+1):N_new] <- Ioptimality_fast(H, X.pob.E[1:n,(N_mejores+1):N_new], 
                                                         p, n, R, Rinv, M, approx = TRUE)
    
    #Actualizacion de las poblaciones elite
    #orden de creciente para G optimalidad (objetivo.ref=1) o
    # I optimalidad (objetivo.ref=2)
    aux <- order(X.pob.E[n+objetivo.ref,], decreasing = FALSE)
    
    #Crear la poblacion elite
    X.Sol[,cont] <- X.pob.E[, aux[1]] # Guardar la mejor solucion.
    X.E <- X.pob.E[, aux[1:N_elite]] # Guardar las mejores N_elite soluciones.
    
    if (mostrar_grafica){
      points(cont, X.Sol[n+objetivo.ref,cont])  
    }
    
    cont  <- cont + 1
  }

  ########## Ya obtenida las soluciones seleccionar la mejor ##########
  
  indicador <- find_nondominated_points(t(X.Sol[c(n + 1,n + 2),]))
  X.Sol <- X.Sol[,indicador]
  X.Soluciones[[muestra]] <- X.Sol
  
  tiempo[[muestra]] <- proc.time()-t

}

save.file <- paste("Extremos_EC/AEC_nObs_", n, "_nCovar_", p, "_Inst_", ii,"_Opt_", objetivo.ref, "_", muestra, "_muestras_.RData", sep = '')
save(X.Soluciones, tiempo, file = save.file)
}
}
