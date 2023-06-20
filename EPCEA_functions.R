##########################################################
# In the following example, I am generating 2 criteria 
# Optimal designs for 3 factors (x1, x2, x3) each with 3 levels 
# (-1, 0, 1). The 2 criteria are D-efficiency and Pure error df
# Let me know if you still have any questions ...
##########################################################


library(gurobi)
source("evaluation.r") # This script contains the objective function.
source("opt.design.r")


##########################################################
# Evalfunc: a function for calculating the criteria values
##########################################################
evalfunc = function(des, H, p, n, R, Rinv, M, gurobi.model) {
  # Changed to reflect our objectives
  Iopt <- Ioptimality_fast(H, des, p, n, R, Rinv, M, approx = TRUE)
  Gopt <- Goptimality_fast(H, des, p, n, R, Rinv, approx = TRUE, gurobi.model)
  -1*c(Iopt, Gopt)
  #c(0, 0)
}

compare = function(newpt, newdes, curpf, curpfdes, nfactor)
{
  # Comparison 2. If a design is not dropped in comparison 1, then it will be
  # compared to the designs which exist in the pareto front P or curpf.
  # Arguments:
  # newpt: objective criteria values of test design (vector).
  # newdes: test design (matrix).
  # curpf: objective criteria values of pareto designs (matrix).
  # curpfdes: designs in pareto front P. (matrix).
  
  # Compare new design with designs in pareto front.
  # g1 is a list
  g1 = round(newpt[1L], 8L) > round(curpf[,1L], 8L) 
  g2 = round(newpt[2L], 8L) > round(curpf[,2L], 8L) 
  
  ge1 = round(newpt[1L], 8L) >= round(curpf[,1L], 8L)
  ge2 = round(newpt[2L], 8L) >= round(curpf[,2L], 8L)
  
  l1 = round(newpt[1L], 8L) < round(curpf[,1L], 8L)
  l2 = round(newpt[2L], 8L) < round(curpf[,2L], 8L)
  
  le1 = round(newpt[1L], 8L) <= round(curpf[,1L], 8L)
  le2 = round(newpt[2L], 8L) <= round(curpf[,2L], 8L)
  
  eq1 = round(newpt[1L], 8L) == round(curpf[,1L], 8L)
  eq2 = round(newpt[2L], 8L) == round(curpf[,2L], 8L)
  
  # cond1 is a vector. If an element is True, then new design does not dominate 
  # the corresponding design. Otherwise, it dominates a design in pareto front.
  cond1 = g1*ge2 + g2*ge1 == 0L
  # cond2 is 0 when new design neither dominates nor is dominated by a design in 
  # the pareto front.
  cond2 = sum(l1*le2 + l2*le1 + eq1*eq2)
  # Condition3 tells me the designs in the pareto front that dominates new design.
  cond3 = seq(1L, nrow(curpf))[cond1] 
  # Initialize cond4. Vector which will tell me which columns to keep from 
  # the matrix containing all pareto optimal designs. 
  # The following command depends on number of factors.
  cond4 = rep(0L, nfactor*length(cond3))
  
  # If there is at least one design in pareto front that dominates the new design.
  if(length(cond3) > 0L)
  {
    for(i in 1L:length(cond3))
    {
      # The following command depends on the number of factors.
      cond4[i*nfactor - nfactor + 1L:nfactor]=cond3[i]*nfactor - nfactor + 1L:nfactor
    }
  }	
  
  # From pareto front, select the designs that dominate the new design.
  # Or, in other words, remove the designs that are dominated by new design.
  newpf = curpf[cond1,]
  # From pareto front designs, select the designs that dominate the new design.
  # Remember, curpfdes will be a matrix with nfactor x npareto columns.
  newpfdes = as.matrix(curpfdes[,cond4])
  
  # If new design neither dominates nor is dominated by a design in 
  # the pareto front.
  if(cond2 == 0L)
  {
    # Append objective values of new design to newpf.
    newpf = rbind(newpf, newpt) 
    # Append new design to newpfdes.
    newpfdes = cbind(newpfdes, newdes)
  }
  
  # Pareto front and Pareto front designs.
  list(matrix(newpf, ncol = 2L), newpfdes)
}


cnt = function(A, B){
  ind = rep(TRUE, nrow(A))
  for (i in 1L:nrow(A)){
    for (j in 1L:nrow(B)){
      if(all(round(A[i,], 4L) == round(B[j,], 4L))){ 
        ind[i] = FALSE
      } 
    }
  }
  ind
}

###########################################
# Pareto based coordinate exchange operator 
###########################################

coexch1 = function(des, level, H, p, n, R, Rinv, M, gurobi.model) {	
  
  cdesmat = as.matrix(des) # copy design.
  nfactor = ncol(cdesmat)
  cbestvals = evalfunc(des, H, p, n, R, Rinv, M, gurobi.model) # compute two objective function values.
  # Make the previous results as a matrix. This is because cpfvals will
  # contain the objective values of all designs in the pareto front P.
  cpfvals = matrix(cbestvals, ncol = 2L) 
  # cchange: controls the number of times to pass through all coordinates.
  cchange = TRUE 
  cpfdes = as.matrix(des)
  
  # We stop when there is no new design that dominates the previously generated
  # good design.
  while(cchange)  
  {
    cchange = FALSE
    for (i in 1L:length(cdesmat)) # proceed through all nrun x nfactor coordinates 
    {
      for (j in level) # test the levels in vector "level"
      { 
        cnewdesmat = cdesmat # Copy design.
        if(cnewdesmat[i] == j) next # If same coordinate, do nohing and proceed to next coordinate
        cnewdesmat[i] = j # exchange coordinate
        cnewbestval = evalfunc(cnewdesmat, H, p, n, R, Rinv, M, gurobi.model) # evaluate new design.
        
        # Comparison 1. If new design dominates the previous design...
        if (any(round(cnewbestval, 8L) > round(cbestvals, 8L)) & all(round(cnewbestval, 8L) >= round(cbestvals, 8L)))
        {
          cbestvals = cnewbestval # Replace objective values
          cchange = TRUE					# Keep switching coordinates.
          cdesmat = cnewdesmat		# Replace the current best design.	
        }					
        # Q: What if current design is dominated by new design? In this case,
        # shouldn't I skip comparison 2?
        
        # Comparison 2. Proceed to second comparison using newly obtained design.
        ctemppf = compare(t(cnewbestval), cnewdesmat, cpfvals, cpfdes, nfactor)
        cpfvals = ctemppf[[1L]]
        cpfdes = ctemppf[[2L]]
        
      } # end for (j in level) 
    }  # end for (i in 1L:length(cdesmat))
  } # end while (cchange)
  
  ctemppf
}

################################################
# EPCEA for a single random start
################################################

coexch2 = function(level, nrun, nfactor, H, p, R, Rinv, M, gurobi.model) {
  
  # (1) Randomly generate a design.
  des = matrix(sample(level, nrun*nfactor, replace = TRUE), ncol = nfactor)
  desfr = as.data.frame(des)
  
  # (2) Perform the Pareto-based coordinate exchange operator in Section 3.2.
  temp = coexch1(des, level, H, p, n = nrun, R, Rinv, M, gurobi.model) 
  psold = psnew = temp[[2L]] # Save new pareto optimal designs.
  pfold = pfnew = temp[[1L]] # Save objective values in pareto set.
  
  # (3) Enhanced elitism operator.
  # For each design in pareto front obtained from des.
  for (i in 1L:nrow(pfold))
  {
    
    des = as.matrix(psold[,(i-1L)*nfactor + 1L:nfactor])
    temp = coexch1(des, level, H, p, n = nrun, R, Rinv, M, gurobi.model) # Apply PBCE operator to this design.
    # Question: We do not skip original design "des" if it is included in psold. Does it matter?
    for(j in 1L:nrow(temp[[1L]]))
    {
      # Compare each pareto front with the pareto front obtained from the inital design.
      temppf = compare(matrix(temp[[1L]][j,], ncol = 2L), as.matrix(temp[[2L]][,(j-1)*nfactor+1L:nfactor]), pfnew, psnew, nfactor)
      pfnew = temppf[[1L]] # Objective values
      psnew = temppf[[2L]]	# Set of designs.		
      
    }
    
  } # end for (i in 1L:nrow(pfold))
  
  # Step (3b). Repeat step (3a) by performing the pareto-based coordinate exchange operator.
  # on the designs in the current Pareto optimal set.
  ind = cnt(pfnew, pfold)
  while(sum(ind) > 0)
  {
    dind = rep(ind, each = nfactor)
    pfdiff = matrix(pfnew[ind,], ncol = 2L)
    psdiff = as.matrix(psnew[,dind]) # Reshape set of designs.
    
    pfold = pfnew
    psold = psnew
    
    for (i in 1L:nrow(pfdiff))
    {
      desmat = as.matrix(psdiff[,(i-1L)*nfactor + 1L:nfactor])
      temp = coexch1(desmat, level, H, p, n = nrun, R, Rinv, M, gurobi.model)
      for(j in 1L:nrow(temp[[1L]]))
      {
        temppf = compare(matrix(temp[[1L]][j,], ncol = 2L), as.matrix(temp[[2L]][,(j-1)*nfactor+1L:nfactor]), pfnew, psnew, nfactor)
        pfnew = temppf[[1L]]
        psnew = temppf[[2L]]			
        
      }
      
    }
    
    ind = cnt(pfnew, pfold)
    
  } # end while(sum(ind) > 0)
  
  list(pfnew, psnew)
}

##################################################
# EPCEA with multiple random starts 
##################################################

coexchpf = function(level, nrun, nfactor, ns, H, p, R, Rinv, M, gurobi.model, my.seed = 5739584) 
{
  set.seed(5739584)
  temp = coexch2(level, nrun, nfactor, H, p, R, Rinv, M, gurobi.model)
  paretodes = temp[[2]]
  paretodescrit = temp[[1]]
  
  for (icount in 2L:ns) 
  {
    
    temp = coexch2(level, nrun, nfactor, H, p, R, Rinv, M, gurobi.model)
    for(i in 1L:nrow(temp[[1L]]))
    {
      temppf = compare(matrix(temp[[1L]][i,], ncol = 2L),as.matrix(temp[[2L]][,(i-1)*nfactor + 1L:nfactor]), paretodescrit, paretodes, nfactor)
      paretodescrit=temppf[[1L]]
      paretodes=temppf[[2L]]			
      
    }
    
  }
  
  list(paretodescrit, paretodes)
}

