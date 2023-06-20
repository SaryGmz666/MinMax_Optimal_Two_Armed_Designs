# Load required packages, libraries and functions.
rm(list=ls())
library(gurobi)
source("evaluation.r") # This script contains the objective function.
source("opt.design.r")
source("EPCEA_functions.R")

#====================================================================
# CARGAR UNA INSTANCIA
#====================================================================
n <- 100 # Numero de sujetos.
p <- 6 # Numero de covariables (including the a column of ones).
file.results.instance <- paste("instances/Results_", n,"_obs_",p, "_covariates.RData",sep='')
load(file.results.instance)

n.instances <- length(H.mat)
pareto.fronts <- list()
pareto.design.sets <- list()
epcea.times <- rep(NA, n.instances)
for (ii in 1:n.instances){
  cat("Wroking on instance", ii, "\n")
  H <- H.mat[[ii]]
  R <- t(H)%*%H # H'H
  Rinv <- solve(R) # (H'H)^{-1}
  M <- moments.matrix(p-1, n.bin = p-1)
  modelo.gurobi <- generate_gurobi_model(p, n.bin = p-1)
  
  init.time <- proc.time()
  design.epcea <- coexchpf(level = c(-1L, 1L), nrun = n, nfactor = 1L, ns = 5, 
                           H, p, R, Rinv, M, modelo.gurobi) 
  final.time <- proc.time() - init.time
  epcea.times[ii] <- final.time[3]
  pareto.fronts[[ii]] <- -1*design.epcea[[1]]
  pareto.design.sets[[ii]] <- design.epcea[[2]]
  file.name <- paste("Frentes_Pareto_EPCEA/nObs_", n, "_nCovar_", p, "_Inst_", ii, ".csv", sep = '')
  write.csv(x = pareto.fronts[[ii]], file = file.name)
}

main.file.name <- paste("Frentes_Pareto_EPCEA/nObs_", n, "_nCovar_", p, ".RData", sep = '')
save(pareto.fronts, pareto.design.sets, file.results.instance, n, p, H.mat, epcea.times,
     file = main.file.name)

#---

load("Frentes_Pareto_EPCEA/nObs_100_nCovar_6.RData")
epcea.times
pareto.fronts[[1]]
