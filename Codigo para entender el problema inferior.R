H <- matrix(sample(c(-1,1), 30, replace = T), ncol = 3)
H[,1] <- 1
R <- t(H)%*%H
max.time = 60 
print.output=1

model <- list()
params <- list(MIPGap=0.001, TimeLimit=max.time,OutputFlag=print.output)

n <- nrow(R)
R22 <- R[2:n, 2:n]
R21 <- R[1, 2:n]
model$obj   <- c(as.vector(R22), R21-rowSums(R22))
model$modelsense <- 'max'

A0 <- diag(1, (n-1)^2)
A1 <- rbind(-A0, -A0, A0, A0)
A22 <- matrix(0, (n-1)^2, n-1)
A23 <- -kronecker(diag(1, n-1), matrix(1, n-1, 1))
A24 <- -kronecker(matrix(1, n-1, 1), diag(1, n-1))
A21 <- -A23-A24
A2 <- rbind(A21, A22, A23, A24)
model$A <- cbind(A1, A2)
model$rhs   <- c(rep(1, (n-1)^2), rep(0, 3*(n-1)^2))
model$sense <- rep('<', nrow(model$A))
model$vtype <- c(rep('C',(n-1)^2), rep('B', n-1))
result <-gurobi(model, params)
# The elements in the output vector z can only take 
# values -1 and +1
z <- 2*result$x[-(1:(n-1)^2)]-1
z <- c(1, z)