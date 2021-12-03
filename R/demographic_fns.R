## Demographic functions -------------------------------------------------------

generation_time<- function(Amat, Fmat){
  # This uses the Bienvenu and Legendre (Am Nat 2015) definition of generation time.
  # T = lambda*v*w/(v*F*w)

  lambda<- Re(eigen(Amat, only.values = T)$values[1])
  w<- Re(eigen(Amat)$vectors[,1]) # right eigenvector of the mean matrix
  v<- Re(eigen(t(Amat))$vectors[,1]) # left eigenvector of the mean matrix

  gentime<- lambda*v%*%w/(v%*%Fmat%*%w)
}

r_nought<- function(Amat, Fmat){
  Umat<- Amat-Fmat
  Nmat<- fundamental_matrix(Umat)
  R0<- Re(eigen(Fmat%*%Nmat, only.values = T)$values[1])
  return(R0)
}

lifespan<- function(Umat, all_ages=T){
  # expected lifespan (in timesteps of the model, not years!) is t(e)*Nmat
  # if all.ages=T, then return the vector of expected lifespan remaining
  # if all.ages=F, then return only the first entry (expected lifespan at the beginning of life)
  Nmat<- fundamental_matrix(Umat)
  eT<- t(rep(1, dim(Umat)[1])) # column of ones
  eta<- eT%*%Nmat
  if (all_ages==F){
    eta<- eta[1]
  }
  return(eta)
}

# Calculate the mean matrix from a list of matrices:
mean_matrix<- function(Aobj){
  # Aobj is a list of matrices. We want to calculate the mean of all indices:

  Aobj_flat<- collapse_mat_list(Aobj)
  meanmat_flat<- apply(Aobj_flat, MARGIN=2, FUN=mean)
  meanmat<- reMat(meanmat_flat)
  return(meanmat)
}
