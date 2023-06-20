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
ii <- 1
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

################################## Inicio de la metaheuristica ##################################

#Control de variables
N_pob         <- 1000            #Cardinalidad de la poblacion
N_elite       <-100  #Cardinalidad de la poblacion elite
#Mariz que contiene las X, los valores de las funciones obj y si pertenece a las soluciones no dominadas
X.pob         <- matrix(rep(0, times = N_pob*(n+3)), ncol = N_pob, nrow = (n+3))
rownames(X.pob) <- c(paste0("Per ", 1:(nrow(X.pob)-3)), "G opt", "I opt", "No Dom")
colnames(X.pob) <- paste0("Sol. ", 1:N_pob)


# Creacion de la poblacion inicial
for(i in 1:N_pob){
  X <- sample(c(-1,1), size = n, replace = TRUE)
  X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
  #Guardamos la solucion
  for(k in 1:(n+3)){
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

#Creacion de la poblacion elite

X.pob.E           <- matrix(rep(0, times = N_elite*(n+3)), ncol = N_elite, nrow = (n+3))
rownames(X.pob.E) <- c(paste0("Per ", 1:(nrow(X.pob.E)-3)), "G opt", "I opt", "Frontera")
colnames(X.pob.E) <- paste0("Sol. ", 1:N_elite)

cont <- 1
dom <- 1

while (cont < N_elite) {
  
  indicador <- find_nondominated_points(t(X.pob[c(n + 1,n + 2),]))
  for(r in 1:N_pob){
    if(indicador[r] == TRUE){
      X.pob[n+3,r] <- 1
    }
  }
  
  for(i in 1:N_pob){
    if(X.pob[n+3, i] == 1){
      for(k in 1:(n+3)){
        X.pob.E[k,cont] <- X.pob[k,i]
        X.pob[k,i] <- 100
      }
      X.pob.E[n+3,cont]<- dom
      cont <- cont + 1
      
      if(cont == N_elite +1){
        break
      }
    } 
  } 
  
  dom <- dom +1
  
}

########## Creacion de nuevas soluciones con respecto a la frecuencia calculada ########## 

sol               <- matrix(rep(0, times = n+3), ncol = 1, nrow = (n+3))
Solucion          <- matrix(rep(0, times = n+3), ncol = 1, nrow = (n+3))
iter <- 3


for(y in 1:iter) {
  X.pob           <- matrix(rep(0, times = N_pob*(n+3)), ncol = N_pob, nrow = (n+3))
  rownames(X.pob) <- c(paste0("Per ", 1:(nrow(X.pob)-3)), "G opt", "I opt", "No Dom")
  colnames(X.pob) <- paste0("Sol. ", 1:N_pob)
  
  Fre <- matrix(rep(0, times = n), ncol = 1, nrow = n)
  
  #Calculo de la frecuencia
  for(i in 1:n){
    for(j in 1:N_elite){
      
      if(X.pob.E[i,j] == 1){
        Fre[i,1] <- Fre[i,1] +1
      }
      
    }
  }
  
  for(i in 1:n){
    Fre[i,1] <- Fre[i,1] / N_elite
  }
  
  # Creacion de una nueva poblacion, con respecto a las frecuencias
  for(i in 1:N_pob){
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
      X.pob[k,i] <- X[k,]
    }
    
    # Calcula el criterio de optimalidad G: Varianza máxima y guardamos
    obj.val <- obj_fun(H, X, max.time = 300, approx = TRUE, print.output = 1) # Solucion del nivel inferior
    X.pob[(n + 1),i] <- obj.val
    
    # Calcular el criterio de optimalidad I: Varianza media y guardamos
    Iopt <- Ioptimality(H, X, approx = TRUE)
    X.pob[(n + 2),i] <- Iopt
  }
  
  #Actualizacion de la poblacion elite
  cont <- 1
  dom <- 1
  
  X.pob.E           <- matrix(rep(0, times = N_elite*(n+3)), ncol = N_elite, nrow = (n+3))
  rownames(X.pob.E) <- c(paste0("Per ", 1:(nrow(X.pob.E)-3)), "G opt", "I opt", "Frontera")
  colnames(X.pob.E) <- paste0("Sol. ", 1:N_elite)
  
  while (cont < N_elite) {
    
    indicador <- find_nondominated_points(t(X.pob[c(n + 1,n + 2),]))
    for(r in 1:N_pob){
      if(indicador[r] == TRUE){
        X.pob[n+3,r] <- 1
      }
    }
    
    for(i in 1:N_pob){
      if(X.pob[n+3, i] == 1){
        for(k in 1:(n+3)){
          if(dom == 1){
            X.pob.E[k,cont] <- X.pob[k,i]
            sol[k,1] <- X.pob[k,i]
            X.pob[k,i] <- 100
          }else{
            X.pob.E[k,cont] <- X.pob[k,i]
            X.pob[k,i] <- 100
          }
        }
        
        X.pob.E[n+3,cont]<- dom
        if(dom == 1){
          Solucion <- cbind(Solucion, sol)
        }
        
        cont <- cont + 1
        
        if(cont == N_elite +1){
          break
        }
      } 
    } 
    
    dom <- dom +1
  }
}

########## Ya obtenida las soluciones no dominadas las volvemos a calcular el frente de pareto 

Solucion <- Solucion[, -1]
indicador <- find_nondominated_points(t(Solucion[c(n + 1,n + 2),]))
Solucion[,indicador]



#Otra idea; agarrar el 5% de las mejores soluciones con respecto a la I optimalidad, y el 5% con respecto a G optimalidad; despues de calculan las frecuencias y
#continuamos con el ciclo