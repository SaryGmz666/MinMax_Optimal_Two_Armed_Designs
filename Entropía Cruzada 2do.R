rm(list=ls())

library(MOEADr)
library(gurobi) # El solucionador gurobi se utiliza en la función objetivo. 
source("evaluation.r") # Este script contiene la función objetivo.
source("opt.design.r")
source("plot_functions.r")


################ Creacion de la instancia ################
n <- 150 # Número de sujetos.
p <- 10 # Número de covariables (incluyendo la columna de unos).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)

################## Inicio y busqueda de la solucion optima  ###################
t1 <- proc.time()

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

proc.time()-t1
################################## Inicio de la metaheuristica ##################################
t <- proc.time()

#Control de variables
N_pob         <- 500     #Cardinalidad de la poblacion
N_elite       <- 50      #Cardinalidad de la poblacion elite
N_new         <- 300     #Cardinalidad de las creaciones nuevas
iter          <- 20      #Numero de iteraciones

#Mariz que contiene las X, los valores de las funciones obj y si pertenece a las soluciones no dominadas
X.pob         <- matrix(rep(0, times = N_pob*(n+3)), ncol = N_pob, nrow = (n+3))
rownames(X.pob) <- c(paste0("Per ", 1:(nrow(X.pob)-3)), "G opt", "I opt", "Prop")
colnames(X.pob) <- paste0("Sol. ", 1:N_pob)

#Creacion de las poblaciones elite y matrices solucion
X.E1           <- matrix(rep(0, times = N_elite*(n+2)), ncol = N_elite, nrow = (n+2))
rownames(X.E1) <- c(paste0("Per ", 1:(nrow(X.E1)-2)), "G opt", "I opt")
colnames(X.E1) <- paste0("Sol. ", 1:N_elite)

X.E2           <- matrix(rep(0, times = N_elite*(n+2)), ncol = N_elite, nrow = (n+2))
rownames(X.E2) <- c(paste0("Per ", 1:(nrow(X.E2)-2)), "G opt", "I opt")
colnames(X.E2) <- paste0("Sol. ", 1:N_elite)

X.E3           <- matrix(rep(0, times = N_elite*(n+2)), ncol = N_elite, nrow = (n+2))
rownames(X.E3) <- c(paste0("Per ", 1:(nrow(X.E3)-2)), "G opt", "I opt")
colnames(X.E3) <- paste0("Sol. ", 1:N_elite)

X.Sol           <- matrix(rep(0, times = (3*(iter + 1))*(n+2)), ncol = (3*(iter + 1)), nrow = (n+2))
rownames(X.Sol) <- c(paste0("Per ", 1:(nrow(X.E1)-2)), "G opt", "I opt")
colnames(X.Sol) <- paste0("Sol. ", 1:(3*(iter + 1)))


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
  
  #Calculo de la proporcion
  X.pob[(n + 3),i] <- (Iopt + obj.val)/2
}


############################## Creacion de la poblacion elite ##############################
cont <- 1

#orden de creciente para G optimalidad e I optimalidad
aux.G <- order(X.pob[101,], decreasing = FALSE)
aux.I <- order(X.pob[102,], decreasing = FALSE)
aux.P <- order(X.pob[103,], decreasing = FALSE)

#Crear la poblacion elite (poblacion para Gopt)
for(i in 1:N_elite){
  for(j in 1:(n+2)){
    if(i == 1){
      X.Sol[j,cont] <- X.pob[j,aux.G[i]]
      X.E1[j,i]       <- X.pob[j,aux.G[i]]
    }else
      X.E1[j,i]       <- X.pob[j,aux.G[i]]
  }
}
cont  <- cont + 1

#Crear la poblacion elite (poblacion para Iopt)
for(i in 1:N_elite){
  for(j in 1:(n+2)){
    if(i == 1){
      X.Sol[j,cont] <- X.pob[j,aux.I[i]]
      X.E2[j,i]       <- X.pob[j,aux.I[i]]
    }else
      X.E2[j,i]       <- X.pob[j,aux.I[i]]
  }
}
cont  <- cont + 1

#Crear la poblacion elite (poblacion para Panza)
for(i in 1:N_elite){
  for(j in 1:(n+2)){
    if(i == 1){
      X.Sol[j,cont] <- X.pob[j,aux.P[i]]
      X.E3[j,i]       <- X.pob[j,aux.P[i]]
    }else
      X.E3[j,i]       <- X.pob[j,aux.P[i]]
  }
}
cont  <- cont + 1


############ Creacion de nuevas soluciones con respecto a la frecuencia calculada ############
for( y in 1:iter ){
  
  X.pob.E           <- matrix(rep(0, times = (3*N_new)*(n+3)), ncol = (3*N_new), nrow = (n+3))
  rownames(X.pob.E) <- c(paste0("Per ", 1:(nrow(X.pob.E)-3)), "G opt", "I opt", "Prop")
  colnames(X.pob.E) <- paste0("Sol. ", 1:(3*N_new))
  
  Fre.E1 <- matrix(rep(0, times = n), ncol = 1, nrow = n)
  Fre.E2 <- matrix(rep(0, times = n), ncol = 1, nrow = n)
  Fre.E3 <- matrix(rep(0, times = n), ncol = 1, nrow = n)
  
  #Calculo de la frecuencia para E1
  for(i in 1:n){
    for(j in 1:N_elite){
      
      if(X.E1[i,j] == 1){
        Fre.E1[i,1] <- Fre.E1[i,1] + 1
      }
      
    }
  }
  
  for(i in 1:n){
    Fre.E1[i,1] <- Fre.E1[i,1] / N_elite
  }
  
  #Calculo de la frecuencia para E2
  for(i in 1:n){
    for(j in 1:N_elite){
      
      if(X.E2[i,j] == 1){
        Fre.E2[i,1] <- Fre.E2[i,1] + 1
      }
      
    }
  }
  
  for(i in 1:n){
    Fre.E2[i,1] <- Fre.E2[i,1] / N_elite
  }
  
  #Calculo de la frecuencia para E3
  for(i in 1:n){
    for(j in 1:N_elite){
      
      if(X.E3[i,j] == 1){
        Fre.E3[i,1] <- Fre.E3[i,1] + 1
      }
      
    }
  }
  
  for(i in 1:n){
    Fre.E3[i,1] <- Fre.E3[i,1] / N_elite
  }
  
  # Creacion de una nueva poblacion, con respecto a las frecuencias
  for(i in 1:(3*N_new)){
    
    if(i < 301){
      
      #Calculo para la primera frecuencia
      X <- sample(runif(n, min=0, max=1), size = n, replace = TRUE)
      X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
      #creacion a traves de la frecuencia
      for(j in 1:n){
        if(Fre.E1[j,1] >= X[j,]){
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
      X.pob.E[(n + 3),i] <- (obj.val+Iopt)/2
      
    }
    if(i>300 & i<601){
      
      #Calculo para la segunda frecuencia
      X <- sample(runif(n, min=0, max=1), size = n, replace = TRUE)
      X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
      #creacion a traves de la frecuencia
      for(j in 1:n){
        if(Fre.E2[j,1] >= X[j,]){
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
      X.pob.E[(n + 3),i] <- (obj.val+Iopt)/2
      
    }
    else{
      
      #Calculo para la segunda frecuencia
      X <- sample(runif(n, min=0, max=1), size = n, replace = TRUE)
      X <- matrix(X, ncol = 1) # La función objetivo funciona con matrices.
      #creacion a traves de la frecuencia
      for(j in 1:n){
        if(Fre.E3[j,1] >= X[j,]){
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
      X.pob.E[(n + 3),i] <- (obj.val+Iopt)/2
      
    }
    
    
  }
  
  #Actualizacion de las poblaciones elite
  #orden de creciente para G optimalidad e I optimalidad
  aux.G <- order(X.pob.E[101,], decreasing = FALSE)
  aux.I <- order(X.pob.E[102,], decreasing = FALSE)
  aux.P <- order(X.pob.E[103,], decreasing = FALSE)
  
  #Crear las poblaciones elite (poblacion para Gopt)
  for(i in 1:N_elite){
    for(j in 1:(n+2)){
      if(i == 1){
        X.Sol[j,cont] <- X.pob.E[j,aux.G[i]]
        X.E1[j,i]       <- X.pob.E[j,aux.G[i]]
      }else
        X.E1[j,i]       <- X.pob.E[j,aux.G[i]]
    }
  }
  cont  <- cont + 1
  #Crear las poblaciones elite (poblacion para Iopt)
  for(i in 1:N_elite){
    for(j in 1:(n+2)){
      if(i == 1){
        X.Sol[j,cont] <- X.pob.E[j,aux.I[i]]
        X.E2[j,i]       <- X.pob.E[j,aux.I[i]]
      }else
        X.E2[j,i]       <- X.pob.E[j,aux.I[i]]
    }
  }
  cont  <- cont + 1
  #Crear la poblacion elite (poblacion para Panza)
  for(i in 1:N_elite){
    for(j in 1:(n+2)){
      if(i == 1){
        X.Sol[j,cont] <- X.pob.E[j,aux.P[i]]
        X.E3[j,i]       <- X.pob.E[j,aux.P[i]]
      }else
        X.E3[j,i]       <- X.pob.E[j,aux.P[i]]
    }
  }
  cont  <- cont + 1
  
}

########## Ya obtenida las soluciones seleccionar la mejor ##########

#Soluciones no dominadas
indicador <- find_nondominated_points(t(X.Sol[c(n + 1,n + 2),]))
X.Sol[c(n + 1,n + 2),indicador]

proc.time()-t
