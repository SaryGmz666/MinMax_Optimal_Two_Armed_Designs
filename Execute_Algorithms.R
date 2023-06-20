####################################################
# Computational experiments to re-produce the results 
# from Zhang et al. (2021).
####################################################

rm(list=ls())

library(gurobi)
source("opt.design.r")
source("evaluation.r")
set.seed(10)

#############################################################
# SECTION A. USE ALGORITHMS TO CONSTRUCT TWO-ARMED DESIGNS
#
#############################################################
mirco.rep <- 5# number of mirco.replication

#ns <- c(60, 100, 120, 150) # number of patients
#ps <- c(4, 6, 10, 15, 20) # number of covariates
ns <- c(100, 150) # number of patients
ps <- c(6, 10, 20) # number of covariates

# Some extra parameters.
max.time.exact <- 300
max.time.approx <- 300
max.time.evaluate <- 300


for(i in 1:length(ns)) {
  for(j in 1:length(ps)) {
    
    # Section 0. Set the number of patients and the number of covariates
    n <- ns[i] # Number of patients
    p <- ps[j] # Number of covariates
    re.all <- NULL # Object to save the results for a single combination (n,p)
    
    # The following objects are explained as follows:
    # H: saves the matrices of covariates of the patients in the sample, 
    #     generated at random.
    # Zs: saves the matrices of covariates of the patients out of the sample, 
    #     generated at random.
    H.mat <- list()
    Zs.mat <- list()
    
    # Section 1. Repeat the experiments 'mirco.rep' times for each combination (n,p) 
    for(it in 1: mirco.rep) {
      # Generate the covariate matrix H for patients in sample. This matrix is
      # used in the objective function.
      H <- matrix(sample(c(-1, 1), n*(p-1), replace=TRUE), n, p-1) 
      H <- cbind(1, H)
      H.mat[[it]] <- H # Save matrix.
      
      
      # Generate the covariate matrix Zs for patients out of sample. This matrix 
      # is used to predict the variances of patients not used to construct the 
      # optimal design x. 
      Zs <- matrix(sample(c(-1, 1), 1000*(p-1), replace=TRUE), 1000, p-1) 
      Zs <- cbind(1, Zs)
      Zs.mat[[it]] <- Zs # Save matrix.
      
      #=========================================================================
      # Lower bound approximation algorithm in Section 3.3 of Zhang et al. (2021)
      #=========================================================================
      
      current.time <- proc.time()
      # Construct design using covariate matrix H as input.
      # Function 'opt.design.approx' is in file "opt.design.R"
      X.opt.approx <- opt.design.approx(H, balance=TRUE, max.time = max.time.approx) 
      # Evaluate design using using covariate matrix Zs as input. 
      # Function 'est_lp' is in file "evaluation.r"
      # The object 're' is the maximum variance for a sample of subjects, not 
      # necessarily in the sample, whose covariate vectors are in matrix Zs.
      # The object 're' is a vector. The first element uses the original matrix
      # in Equation (7) of the Zhang et al. (2021). The second element uses the surrogate 
      # matrix in Equation (9) of this paper. 
      
      re <- est_lp(H, X.opt.approx, Zs, max.time.evaluate)
      # Remark: If the number of covariates is smaller than 50, the function 
      # outputs the maximum variance over the whole feasible space of covariates.
      # Otherwise, it outputs the maximum variance over the covariate vectors
      # in matrix Zs.
      if(nrow(re)>0) {
        time.duration <- proc.time()-current.time
        re$method <- "LB_APPROX"
        re$IT <- it
        re$time <- time.duration[3]
      }
      # Save the results of this iteration for this algorithm.
      re.all <- rbind(re.all, re) 
      
    
      #=========================================================================
      # Exact algorithm in Section 3.2 of Zhang et al. (2021)
      #=========================================================================
      
      current.time <- proc.time()
      # Construct design using covariate matrix H as input.
      # Function 'opt.design.iter' is in file "opt.design.R"
      X.opt.iter <- opt.design.iter(H, balance=TRUE, max.time = max.time.exact) # 500
      # Evaluate design using using covariate matrix Zs as input. 
      # Function 'est_lp' is in file "evaluation.r"
      # The object 're' is the maximum variance for a sample of subjects, not 
      # necessarily in the sample, whose covariate vectors are in matrix Zs.
      # The object 're' is a vector. The first element uses the original matrix
      # in Equation (7) of the Zhang et al. (2021). The second element uses the surrogate 
      # matrix in Equation (9) of this paper. 
      
      re <- est_lp(H, X.opt.iter, Zs, max.time.evaluate)
      # Remark: If the number of covariates is smaller than 50, the function 
      # outputs the maximum variance over the whole feasible space of covariates.
      # Otherwise, it outputs the maximum variance over the covariate vectors
      # in matrix Zs.
      if(nrow(re)>0) {
        time.duration <- proc.time()-current.time
        re$method <- "EXACT"
        re$IT <- it
        re$time <- time.duration[3]
      }
      # Save the results of this iteration for this algorithm.
      re.all <- rbind(re.all, re)
      
      #=========================================================================
      # Generate two-treatment designs at random.
      #=========================================================================
      
      nrep <- 100 # Number of designs to generate.
      
      # random balanced design
      current.time <- proc.time()
      # Vector of n/2 -1s and n/2 +1s, where n is the number of patients.
      x <- rep(c(-1, 1), each=n/2) 
      # Generate a random sample of 100 vectors.
      # Matrix X.bl is ns x nrep, where ns is the number of patients in sample.
      X.bl <- replicate(nrep, sample(x))
      # Evaluate each of the randomly generated designs.
      re <- est_lp(H, X.bl, Zs, max.time.evaluate)
      # Remark: If the number of covariates is smaller than 50, the function 
      # outputs the maximum variance over the whole feasible space of covariates.
      # Otherwise, it outputs the maximum variance over the covariate vectors
      # in matrix Zs.
      time.duration <- proc.time()-current.time
      re$method <- "RAND"
      re$IT <- it
      re$time <- time.duration[3]
      # Identify the best design in terms of the original objective function.
      best.design <- which.min(re$obj)
      re.all <- rbind(re.all, re[best.design,])
    }
    
    # ==========================================================================
    # Save results and important objects.
    # ==========================================================================
    re.all$n <- n
    re.all$p <- p
    
    file.results.instance <- paste("Results_", n,"_obs_",p, "_covariates.RData",sep='')
    save(re.all, H.mat, Zs.mat, file = file.results.instance)
    
    print(c(n,p))
    
  }
}

###################################################################
# SECTION B. COMPILE THE RESULTS FROM ALL INSTANCES AND ALGORTHMS.
#
###################################################################

collect.re <- NULL
for(i in 1:length(ns)) {
  for(j in 1:length(ps)) {
    n <- ns[i]
    p <- ps[j]
    file.results.instance <- paste("Results_", n,"_obs_",p, "_covariates.RData",sep='')
    load(file.results.instance)
    collect.re <- rbind(collect.re, re.all)
  }
}

write.csv(x = collect.re,file = "allresults.csv")



