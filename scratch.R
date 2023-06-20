#=========================================================================
# Additive model algorithm which ignores the patients' covariates
# when constructing the design. This algorithm solves the problem in
# Equation (6) of Zhang et al. (2021) using Gurobi.
#=========================================================================

current.time <- proc.time()
# Construct design using covariate matrix H as input.
# Function 'opt.design.approx' is in file "opt.design.R"
X.opt.approx <- opt.design.approx1(H, balance=TRUE)
# Evaluate design using using covariate matrix Zs as input. 
# Function 'est_lp' is in file "evaluation.r"
# The object 're' is the maximum variance for a sample of subjects, not 
# necessarily in the sample, whose covariate vectors are in matrix Zs.
# The object 're' is a vector. The first element uses the original matrix
# in Equation (7) of the Zhang et al. (2021). The second element uses the surrogate 
# matrix in Equation (9) of this paper. 
re <- est_lp(H, X.opt.approx, Zs)
if(nrow(re)>0) {
  time.duration <- proc.time()-current.time
  re$method <- "ADDITIVE"
  re$IT <- it
  re$time <- time.duration[3]
}
re.all <- rbind(re.all, re)

