est <- function(H, X, Zs) {
  Sigmainvs <- list()
  for(j in 1:ncol(X)) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    Sigmainvs[[j]] <- solve(Sigma)
  }
  zs <- matrix(0, 0, ncol(X))
  for(i in 1:nrow(Zs)) {
      z <- Zs[i, ]
      zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
      
  }
  return(apply(zs, 1, mean))
}



est_lp <- function(H, X, Zs,max.time = 300,print.output = 0) {
#############################################
# Objective function.
# If the number of covariates is smaller than 50,
# the function outputs the maximum variance over the
# whole feasible space of covariates.
# Otherwise, it outputs the maximum variance 
# over the covariate vectors in matrix Zs.
#############################################

 p <- ncol(H)
 ss <- ncol(X)
 max.obj <- rep(999, ss)
 max.obj.approx <- rep(999, ss)
 if(p <= 50) {
  for(j in 1:ss) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
    Sigmainv <- try(solve(Sigma), TRUE) 
    # If matrix is invertible, then use it to find z which maximizes the variance
    # over the whole feasible covariate space. 
    if(is.matrix(Sigmainv)) {
     zast <- lp_sdp(Sigmainv,max.time,print.output)
     max.obj[j] <- t(zast)%*% Sigmainv %*% zast
    } 
    Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
    zast <- lp_sdp(Sigma.approx,max.time,print.output)
    max.obj.approx[j] <- t(zast)%*% Sigma.approx %*% zast
  }
 re <- data.frame(obj=max.obj, obj.approx=max.obj.approx)
 }
 if(p > 50) {
    Sigmainvs <- list()
    Sigma.approx <- list()
    for(j in 1:ncol(X)) {
      x <- X[,j]
      n <- length(x)
      D <- diag(x)
      R <- t(H)%*%H
      Rinv <- solve(R)
      Rx <- t(H)%*%D%*%H
      Sigma <- R-Rx%*%Rinv%*%Rx
      Sigmainvs[[j]] <- solve(Sigma)
      Sigma.approx[[j]] <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
   }
   zs <- matrix(0, 0, ncol(X))
   for(i in 1:nrow(Zs)) {
     z <- Zs[i, ]
     zs <- rbind(zs, sapply(Sigmainvs, function(U) t(z)%*%U%*%z))
   }
   max.obj <- apply(zs, 2, max)
   zs <- matrix(0, 0, ncol(X))
   for(i in 1:nrow(Zs)) {
     z <- Zs[i, ]
     zs <- rbind(zs, sapply(Sigma.approx, function(U) t(z)%*%U%*%z))
   }
   max.obj.approx <- apply(zs, 2, max)
   re <- data.frame(obj=max.obj, obj.approx=max.obj.approx)
 }
 return(re)
}

############################
# Custom objective function
############################

obj_fun <- function(H, X, max.time = 300, approx = FALSE, print.output=1) {
#############################################
# Objective function.
# If the number of covariates is smaller than 50,
# the function outputs the maximum variance over the
# whole feasible space of covariates.
# Otherwise, it outputs the maximum variance 
# over the covariate vectors in matrix Zs.
#############################################

 p <- ncol(H)
 ss <- ncol(X)
 max.obj <- rep(999, ss)

  for(j in 1:ss) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H # H'H
    Rinv <- solve(R) # (H'H)^{-1}
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    
    if ( approx ){
      Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv # (H'H)^{-1} + phi.
      zast <- lp_sdp(Sigma.approx,max.time,print.output) # lp_sdp is in "opt.design.r"
      max.obj[j] <- t(zast)%*% Sigma.approx %*% zast
    } else {
      # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
      Sigmainv <- try(solve(Sigma), TRUE) 
      # If matrix is invertible, then use it to find z which maximizes the variance
      # over the whole feasble covariate space. 
      if(is.matrix(Sigmainv)) {
        zast <- lp_sdp(Sigmainv,max.time,print.output)
        max.obj[j] <- t(zast)%*% Sigmainv %*% zast
      }
    }

  }
 
 return(max.obj)
}

moments.matrix <- function(m, n.bin = 0, type = 'l'){
  #=============================================================================
  # moments.matrix: Compute the moments matrix for linear and quadratic models.
  #                 For linear models, the function can output the moments matrix
  #                 for a model containing continuous and binary covariates in 
  #                 {-1, +1}.
  # Inputs: 
  #       p: number of covariates.
  #       n.bin: number of binary covariates.
  #       type: 'l' for linear model and 'q' for quadratic model.
  #
  # Output:
  #       M: Moments matrix.
  #=============================================================================
  if (type == 'l'){
    M <- 1/3*diag(m+1)
    M[1,1] <- 1
    
    # If there are binary variables, assign them at the end.
    if (n.bin > 0){
      las.elem <- tail(1:(m+1),n.bin)
      diag(M)[las.elem] <- 1
    }
    return(M)
  } else if (type == 'q' && n.bin == 0){
    p <- (m+1)*(m+2)/2 
    M <- diag(p)
    diag(M) <- c(1, rep(1/3, m), rep(1/9, m*(m-1)/2), rep(1/5-1/9, m))
    las.elem <- tail(1:p,m)
    M[las.elem,las.elem] <- M[las.elem,las.elem] + 1/9
    M[1,las.elem] <- 1/3
    M[las.elem,1] <- 1/3
    return(M)
  } else {print('Undefined model')}
  
}

Ioptimality <- function(H, X, n.bin = 0, approx = FALSE) {
  #############################################
  # Objective function.
  # If the number of covariates is smaller than 50,
  # the function outputs the maximum variance over the
  # whole feasible space of covariates.
  # Otherwise, it outputs the maximum variance 
  # over the covariate vectors in matrix Zs.
  #############################################
  
  p <- ncol(H)
  ss <- ncol(X)
  ave.obj <- rep(999, ss)
  M <- moments.matrix(p-1, n.bin)
  
  for(j in 1:ss) {
    x <- X[,j]
    n <- length(x)
    D <- diag(x)
    R <- t(H)%*%H
    Rinv <- solve(R)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    
    if ( approx ){
      Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
      ave.obj[j] <- sum(diag(Sigma.approx %*% M))
    } else {
      # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
      Sigmainv <- try(solve(Sigma), TRUE) 
      # If matrix is invertible, then use it to find z which maximizes the variance
      # over the whole feasble covariate space. 
      if(is.matrix(Sigmainv)) {
        ave.obj[j] <- sum(diag(Sigmainv %*% M))
      }
    }
    
  }
  return(ave.obj)
}

sing_Iopt <- function(x, p, n, ss, H, R, Rinv, M, approx = TRUE){
  D <- diag(x)
  Rx <- t(H)%*%D%*%H
  Sigma <- R-Rx%*%Rinv%*%Rx
  
  if ( approx ){
    Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv
    obj.val <- sum(diag(Sigma.approx %*% M))
  } else {
    # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
    Sigmainv <- try(solve(Sigma), TRUE) 
    # If matrix is invertible, then use it to find z which maximizes the variance
    # over the whole feasble covariate space. 
    if(is.matrix(Sigmainv)) {
      obj.val <- sum(diag(Sigmainv %*% M))
    }
  }
  return(obj.val)
}

Ioptimality_fast <- function(H, X, p, n, R, Rinv, M, approx = FALSE) {
  #############################################
  # Objective function.
  # If the number of covariates is smaller than 50,
  # the function outputs the maximum variance over the
  # whole feasible space of covariates.
  # Otherwise, it outputs the maximum variance 
  # over the covariate vectors in matrix Zs.
  #############################################
  
  ss <- ncol(X)
  ave.obj <- apply(X, MARGIN = 2, FUN = sing_Iopt, p=p, n=n, ss=ss, H=H, R=R, 
                  Rinv=Rinv, M=M, approx = approx)
  
  return(ave.obj)
}

solve.sub.problem <- function(R, p, gurobi.model){
  R22 <- R[2:p, 2:p]
  R21 <- R[1, 2:p]
  gurobi.model$model$obj   <- c(as.vector(R22), R21-rowSums(R22))
  result <- gurobi(gurobi.model$model, gurobi.model$params)
  z <- 2*result$x[-(1:(p-1)^2)]-1
  z <- c(1, z)
  return(z)
}

Goptimality_fast <- function(H, X, p, n, R, Rinv, approx = FALSE, gurobi.model) {
  #############################################
  # Objective function.
  # If the number of covariates is smaller than 50,
  # the function outputs the maximum variance over the
  # whole feasible space of covariates.
  # Otherwise, it outputs the maximum variance 
  # over the covariate vectors in matrix Zs.
  #############################################
  
  ss <- ncol(X)
  max.obj <- rep(999, ss)
  
  for(j in 1:ss) {
    x <- X[,j]
    D <- diag(x)
    Rx <- t(H)%*%D%*%H
    Sigma <- R-Rx%*%Rinv%*%Rx
    
    if ( approx ){
      Sigma.approx <- Rinv+Rinv%*%Rx%*%Rinv%*%Rx%*%Rinv # (H'H)^{-1} + phi.
      #zast <- lp_sdp(Sigma.approx,max.time,print.output) # lp_sdp is in "opt.design.r"
      zast <- solve.sub.problem(Sigma.approx,p,gurobi.model) 
      max.obj[j] <- t(zast)%*% Sigma.approx %*% zast
    } else {
      # Compute the original covariance matrix in Equation (7) of Zhang et al. (2021)
      Sigmainv <- try(solve(Sigma), TRUE) 
      # If matrix is invertible, then use it to find z which maximizes the variance
      # over the whole feasble covariate space. 
      if(is.matrix(Sigmainv)) {
        #zast <- lp_sdp(Sigmainv,max.time,print.output)
        zast <- solve.sub.problem(Sigma.approx,p,gurobi.model)
        max.obj[j] <- t(zast)%*% Sigmainv %*% zast
      }
    }
    
  }
  
  return(max.obj)
}

generate_gurobi_model <- function(p, n.bin = 0, max.time = 300, print.output = 0){
  #=============================================================================
  # generate_gurobi_model: Initialize gurobi model.
  #
  # Inputs: 
  #       p: number of covariates.
  #       n.bin: number of binary covariates.
  #       max.time: maximum time for the gurobi solver.
  #       print..output: 1 if yes; 0 otherwise,
  #
  # Output:
  #       result: Gurobi model (list) and parameters for the model (list).
  #=============================================================================
  
  m <- p-1 # Number of covariates.
  n.vars <- p-1 + (p-1)^2 # Number of variables in optimization problem.
  n.bil.vars <- (p-1)^2 # Number of bilinear/auxiliar variables.
  result <- list()
  
  # Intialize model.
  model <- list()
  model$modelsense <- 'max'
  
  # MODELO LINEAR PARA VARIABLES BINARIAS.========================================
  if (m == n.bin){
    
    A0 <- diag(1, n.bil.vars)
    A1 <- rbind(-A0, -A0, A0, A0)
    A22 <- matrix(0, n.bil.vars, m)
    A23 <- -kronecker(diag(1, m), matrix(1, m, 1))
    A24 <- -kronecker(matrix(1, m, 1), diag(1, m))
    A21 <- -A23-A24
    A2 <- rbind(A21, A22, A23, A24)
    model$A <- cbind(A1, A2)
    model$rhs   <- c(rep(1, n.bil.vars), rep(0, 3*n.bil.vars))
    model$sense <- rep('<', nrow(model$A))
    model$vtype <- c(rep('C',n.bil.vars), rep('B', m))
    
    # Define parameters of model.
    params <- list(TimeLimit=max.time,OutputFlag=print.output)
    
  } else {
    # If at least one decision variable is continuous.
    # Define domain of decision variables.
    var.type <- rep('C', m)
    if (n.bin > 0) {var.type[tail(1:m,n.bin)] <- rep('B', n.bin) }
    
    # Define linear constraints of model.
    temp.A <- rbind(diag(1, n.bil.vars), -1*diag(1, n.bil.vars))
    model$A <- cbind(temp.A, matrix(0, ncol = m, nrow = nrow(temp.A)))
    model$rhs   <- c(rep(1, n.bil.vars), rep(0, n.bil.vars))
    model$sense <- rep('<', nrow(model$A))
    model$vtype <- c(rep('C',n.bil.vars), var.type)
    
    # Define bilinear constraints.
    citer <- 1
    all_qcl <- list()
    for(i in 1:(p-1)){
      for (j in 1:(p-1)){
        qcl <- list()
        qcl$Qc <- spMatrix(n.vars, n.vars, c(n.bil.vars+i), c(n.bil.vars+j), c(1.0))
        qcl$rhs <- 0
        qcl$sense <- c('=')
        tmp.vec <- rep(0,n.vars)
        tmp.vec[citer] <- -1
        qcl$q <- tmp.vec
        all_qcl[[citer]] <- qcl
        citer <- citer + 1
      }
    }
    model$quadcon <- all_qcl
    
    # Define parameters of model.
    params <- list(NonConvex=2,OutputFlag=print.output,TimeLimit=max.time)
    
  }
  
  result$model <- model
  result$params <- params
  return(result)
}