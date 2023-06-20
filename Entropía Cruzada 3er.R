rm(list=ls())

library(MOEADr)
library(gurobi) # El solucionador gurobi se utiliza en la función objetivo. 
source("evaluation.r") # Este script contiene la función objetivo.
source("opt.design.r")
source("plot_functions.r")


################ Creacion de la instancia ################
n <- 100 # Número de sujetos.
p <- 6 # Número de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)

################## Inicio y busqueda de la solucion optima  ###################
t1 <- proc.time()

ii <- 5
H <- H.mat[[ii]]
Zs <- Zs.mat[[ii]]
#X.opt.iter <- opt.design.iter(H, balance=TRUE, max.time = 20) # Exacto
X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = 300) #Aproximado
# Calcula el criterio de optimalidad G: Varianza máxima y guardamos
obj.val.opt <- obj_fun(H, X.opt.approx, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
obj.val.opt

# Calcular el criterio de optimalidad I: Varianza media y guardamos
Iopt <- Ioptimality(H, X.opt.approx, approx = TRUE)
Iopt

proc.time()-t1

################################################### Inicio de la metaheuristica ###################################################
t <- proc.time()

#Control de variables
N_pob         <- 500     #Cardinalidad de la poblacion
N_elite       <- 50      #Cardinalidad de la poblacion elite
N_new         <- 400     #Cardinalidad de las creaciones nuevas
iter          <- 40      #Numero de iteraciones
EntroCruz     <- FALSE    #True es para entropia cruzada para I-opt, False es para G-opt

###### Alta de Matrices ######
#Mariz que contiene las X, los valores de las funciones obj y si pertenece a las soluciones no dominadas
X.pob           <- matrix(rep(0, times = N_pob*(n+2)), ncol = N_pob, nrow = (n+2))
rownames(X.pob) <- c(paste0("Per ", 1:(nrow(X.pob)-2)), "G opt", "I opt")
colnames(X.pob) <- paste0("Sol. ", 1:N_pob)

#Creacion de las poblaciones elite y matrices solucion
X.E           <- matrix(rep(0, times = N_elite*(n+2)), ncol = N_elite, nrow = (n+2))
rownames(X.E) <- c(paste0("Per ", 1:(nrow(X.E)-2)), "G opt", "I opt")
colnames(X.E) <- paste0("Sol. ", 1:N_elite)

X.Sol           <- matrix(rep(0, times = (iter + 1)*(n+2)), ncol = (iter + 1), nrow = (n+2))
rownames(X.Sol) <- c(paste0("Per ", 1:(nrow(X.Sol)-2)), "G opt", "I opt")
colnames(X.Sol) <- paste0("Sol. ", 1:(iter + 1))


############################## Creacion de la poblacion inicial ##############################
for(i in 1:N_pob){
  X <- sample(c(-1,1), size = n, replace = TRUE)
  X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
  #Guardamos la solucion
  for(k in 1:(n)){
    if(k<=n){
      X.pob[k,i] <- X[k,]
    }
  }
  
  # Calcula el criterio de optimalidad G: Varianza máxima y guardamos
  obj.val <- obj_fun(H, X, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
  X.pob[(n + 1),i] <- obj.val
  
  # Calcular el criterio de optimalidad I: Varianza media y guardamos
  Iopt <- Ioptimality(H, X, approx = TRUE)
  X.pob[(n + 2),i] <- Iopt
}

############################## Creacion de la poblacion elite ##############################
cont <- 1

#orden de creciente para G optimalidad e I optimalidad
if(EntroCruz == TRUE){
  aux <- order(X.pob[102,], decreasing = FALSE)
}else{
  aux <- order(X.pob[101,], decreasing = FALSE)
}

#Crear la poblacion elite
for(i in 1:N_elite){
  for(j in 1:(n+2)){
    if(i == 1){
      X.Sol[j,cont]  <- X.pob[j,aux[i]]
      X.E[j,i]       <- X.pob[j,aux[i]]
    }else
      X.E[j,i]       <- X.pob[j,aux[i]]
  }
}
cont  <- cont + 1


############ Creacion de nuevas soluciones con respecto a la frecuencia calculada ############
for( y in 1:iter ){
  
  X.pob.E           <- matrix(rep(0, times = N_new*(n+2)), ncol = N_new, nrow = (n+2))
  rownames(X.pob.E) <- c(paste0("Per ", 1:(nrow(X.pob.E)-2)), "G opt", "I opt")
  colnames(X.pob.E) <- paste0("Sol. ", 1:N_new)
  
  Fre <- matrix(rep(0, times = n), ncol = 1, nrow = n)
  
  #Calculo de la frecuencia
  for(i in 1:n){
    for(j in 1:N_elite){
      
      if(X.E[i,j] == 1){
        Fre[i,1] <- Fre[i,1] + 1
      }
      
    }
  }
  
  for(i in 1:n){
    Fre[i,1] <- Fre[i,1] / N_elite
  }
  
  # Creacion de una nueva poblacion, con respecto a las frecuencias
  for(i in 1:(N_new)){
    
    if(i <= 10){ #Se guardan las primeras 10 soluciones de X.E
      
      for(j in 1:(n+2)){
        X.pob.E[j,i] <- X.E[j,i]
      }
      
    }else{
      
      #Calculo para la primera frecuencia
      X <- sample(runif(n, min=0, max=1), size = n, replace = TRUE)
      X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
      #creacion a traves de la frecuencia
      for(j in 1:n){
        if(Fre[j,1] >= X[j,]){
          X[j,] <- 1
        }else{
          X[j,] <- -1
        }
      }
      
      #Guardamos la solucion
      for(k in 1:n){
        X.pob.E[k,i] <- X[k,]
      }
      # Calcula el criterio de optimalidad G: Varianza máxima y guardamos
      obj.val <- obj_fun(H, X, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
      X.pob.E[(n + 1),i] <- obj.val
      # Calcular el criterio de optimalidad I: Varianza media y guardamos
      Iopt <- Ioptimality(H, X, approx = TRUE)
      X.pob.E[(n + 2),i] <- Iopt
      
    }
  }
  
  #Actualizacion de las poblaciones elite
  #orden de creciente para G optimalidad e I optimalidad
  if(EntroCruz == TRUE){
    aux <- order(X.pob.E[102,], decreasing = FALSE)
  }else{
    aux <- order(X.pob.E[101,], decreasing = FALSE)
  }
  
   
  #Crear la poblacion elite
  for(i in 1:N_elite){
    for(j in 1:(n+2)){
      if(i == 1){
        X.Sol[j,cont]  <- X.pob.E[j,aux[i]]
        X.E[j,i]       <- X.pob.E[j,aux[i]]
      }else
        X.E[j,i]       <- X.pob.E[j,aux[i]]
    }
  }
  cont  <- cont + 1
  
}

########## Ya obtenida las soluciones seleccionar la mejor ##########

#Soluciones no dominadas
indicador <- find_nondominated_points(t(X.Sol[c(n + 1,n + 2),]))
X.Sol[c(n + 1,n + 2),indicador]

proc.time()-t
