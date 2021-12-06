## Auxiliary functions ---------------------------------------------------------
# A function for variance assuming a complete sample:
variance_complete <- function(x) {
  xx<- x[!is.na(x)]
  mean((xx - mean(xx))^2)
}

# a function for the covariance object:
cov_matrix<- function(Aobj){
  # covariance matrix is: C = Expected[vec(A)*t(vec(A))] - vec(Amean)*t(vec(Amean))
  Amean<- apply(Aobj, 2, mean) # calculate the mean matrix
  second_term<- Amean%*%t(Amean)
  first_term_elements<- list() # initialize the holder for the part inside square brackets
  for (i in 1:dim(Aobj)[1]){
    first_term_elements[[i]]<- Aobj[i,]%*%t(Aobj[i,])
  }
  # first_term_elements is a list of matrices. We want the average matrix from those:
  first_term<- apply(simplify2array(first_term_elements), 1:2, mean)
  covmat<- first_term - second_term
  return(covmat)
}

# A function to reconstruct a matrix that was collapsed column-wise into a row vector:
reMat<- function(vecM,j=NULL){
  # vecM is a matrix where each row is vec(Aj), and j is the row index to be re-constructed

  # if vecM is a single row, we need to write it a little differently:
  if (is.null(dim(vecM))){
    matdim<- sqrt(length(vecM)) # this is only for making square matrices, so
    # the matrix dimensions are the square root of the length of the row vector.
    matted<- matrix(vecM, matdim, matdim)
  } else {
    if(is.null(j)){
      warning("You have not specified which row to reMat, assuming j=1.")
      j<- 1
    }
    matdim<- sqrt(length(vecM[j,]))
    matted<- matrix(vecM[j,], matdim, matdim)
  }
  return(matted)
}

# A function to collapse a list of matrices into rows:
collapse_mat_list<- function(Alist){
  # check what size the Aobj matrix needs to be:
  matdim<- dim(Alist[[1]])[1]
  # initialize the Aobj matrix:
  Aobj<-matrix(NA, nrow=length(Alist), ncol=matdim^2)

  for (i in 1:length(Alist)){
    Aobj[i,]<- as.vector(Alist[[i]])
  }
  return(Aobj)
}

# A function to compute the variance in lambda across a set of matrices, with some parameters fixed at their mean value:
lamVar<- function(Aobj, which.fixed=NULL) {

  if (is.list(Aobj)){
    vecM<- collapse_mat_list(Aobj)
  } else {
    vecM<- Aobj
  }

  Mfix<- vecM
  # for each of the parameters to be held fixed, set it equal to the mean for
  # all matrices in the dataset
  if (length(which.fixed)>0){
    for(j in which.fixed) {
      Mfix[,j]<- mean(vecM[,j])
    }
  }

  lambdas<- numeric(nrow(Mfix)) # initialize lambdas
  # calculate the lambda of each population matrix (with the fixed parameters
  # held constant at their means)
  for(j in 1:nrow(Mfix)) {
    Mj<- reMat(Mfix,j)
    lambdas[j]<- Re(eigen(Mj, only.values = T)$values[1])
  }
  # calculate the variance in lambda
  return(variance_complete(lambdas))
}

# A function to compute the difference in lambda across a set of matrices, with some parameters fixed at their mean value:
lamDiff<- function(Aobj, which.fixed=NULL) {
  # Aobj can either be a list of matrices, where Aobj[[1]] is Aref and Aobj[[2]] is Atest
  # or have 2 rows that are the vec(Aref) and vec(Atest)

  if (is.list(Aobj)){
    if (length(Aobj)>2){
      warning("Aobj contains more than 2 matrices. lamDiff assumes Aobj[[1]] is Aref and Aobj[[2]] is Atest.")
    }
    else if (length(Aobj)<2){
      stop("Aobj contains fewer than 2 matrices, so we cannot compute the difference in lambda.")
    }
    Mtest<- rbind(as.vector(Aobj[[1]]), as.vector(Aobj[[2]])) # convert to the vec'ed format
  } else { # if Aobj isn't a list, then we assume it's a matrix.
    if (dim(Aobj)[1]>2){
      warning("Aobj contains more than 2 matrices. lamDiff assumes Aobj[1,] is vec(Aref) and Aobj[,2] is vec(Atest).")
    } else if (dim(Aobj)[1]<2){
      stop("Aobj contains fewer than 2 matrices, so we cannot compute the difference in lambda.")
    }
    Mtest<- Aobj[1:2,]
  }

  # for each of the parameters to be held fixed, set it equal to the corresponding
  # value in the reference matrix
  if (length(which.fixed)>0){
    for(j in which.fixed) {
      Mtest[2,j]<- Aobj[1,j]
    }
  }

  lambdas<- numeric(nrow(Mtest)) # initialize lambdas
  # calculate the lambda of each population matrix (with the fixed parameters
  # held constant at their means)
  for(j in 1:nrow(Mtest)) {
    Mj<- reMat(Mtest,j)
    lambdas[j]<- Re(eigen(Mj, only.values = T)$values[1])
  }
  # calculate the difference in lambda as Atest-Aref:
  return(lambdas[2]-lambdas[1])
}

# recursive function for making a Gmatrix (from Poelwijk, Krishna, Ranganathan 2016 paper):
make.Gmatrix<-function(n) {
  G<- 1;
  for(k in 1:n) {
    topMat<- cbind(G, 0*G)
    botMat<- cbind(-G, G)
    G<- rbind(topMat, botMat)}
  return(G)
}

# A lower-level function for calculating the responses of a matrix, according to some other function:
calc_matrix_responses<- function(Aobj, ind_vary, FUN, maxint="all"){
  # count how many indices vary:
  n_vary<- length(ind_vary)

  if (maxint=="all"){ # use expand.grid:
    # find all the combinations of ind_vary for generating variances (or differences) in lambda (the nu's):
    bin_order<- expand.grid(replicate(n_vary, c(0,1), simplify=FALSE), KEEP.OUT.ATTRS = FALSE)
  } else {
    if (!is.numeric(maxint)){
      stop("Input variable maxint, the maximum interaction order to calculate, must either be ''all'' or a numeric variable.")
    } else if (maxint>=n_vary){
      warning("You have requested an interaction order that is greater than or equal to the number of parameters that vary, defaulting to maxint='all'.")
      # find all the combinations of ind_vary for generating variances (or differences) in lambda (the nu's):
      bin_order<- expand.grid(replicate(n_vary, c(0,1), simplify=FALSE), KEEP.OUT.ATTRS = FALSE)
      maxint<- "all"
    } else {
      list_ind_vary<- list()
      list_ind_vary[[1]]<- vector() # row 1 is empty, fix all of the variables
      for (i in 1:maxint){
        newrows<- t(utils::combn(ind_vary, i))
        newrows<- lapply(seq_len(nrow(newrows)), function(i) newrows[i,])
        list_ind_vary<- c(list_ind_vary, newrows)
      }
    }
  }
  # okay, now we have either bin_order or list_ind_vary, both of which have the
  # necessary information to calculate the nu vector

  if (maxint=="all"){ # this is the bin_order method.
    # convert bin_order to list_ind_vary:
    list_ind_vary<- list()
    list_ind_vary[[1]]<- vector()
    for (i in 2:dim(bin_order)[1]){
      newlistitem<- ind_vary[which(bin_order[i,]==1)]
      list_ind_vary[[i]]<- newlistitem
    }
    rm(bin_order)
  }

  # the number of nus is the number of list entries in list_ind_vary
  n_nus<- length(list_ind_vary)
  # initialize the holder for the nu values:
  nus<- numeric(length=n_nus)

  # loop for calculating the nus:
  for (i in 1:n_nus){
    ind_i<- list_ind_vary[[i]]
    # fix.these should be the indices that vary that are not in ind_i:
    fix.these<- ind_vary[!(ind_vary %in% ind_i)]
    nus[i]<- FUN(Aobj, fix.these)
  }
  # need to return the nus and the list_ind_vary:
  results<- list(nus=nus, list_ind_vary=list_ind_vary)
  return(results)

}

# An alternative to using the Gmatrix approach, this converts from responses to effects:
calc_matrix_effects<- function(responses, list_ind_vary){
  effects<- vector(length=length(list_ind_vary))
  effects[1]<- responses[1]
  term_order<- sapply(list_ind_vary, length)

  for (i in 2:length(effects)){# go row-wise, starting with row 2
    ind_i<- list_ind_vary[[i]] # which indices were varying for this row?

    thisrow<- rep(0, length=length(effects))
    thisrow[i]<- 1 # always need to include the response of the term of interest

    # we only need to check the entries that are lower order than term_order[i]
    tocheck<- which(term_order<term_order[i])

    for (j in tocheck){
      ind_j<- list_ind_vary[[j]]

      if (sum(ind_j %in% ind_i)==length(ind_j)){
        thisrow[j]<- ifelse(term_order[i]%%2 == term_order[j]%%2, 1, -1)
      }
    }
    # multiply this row by the response vector to get the effects entry
    effects[i]<- thisrow%*%responses

  }
  return(effects)
}

# Run matrix checks on an Aobj array (each row is a matrix, collapsed column-wise):
run_matrix_checks<- function(Aobj){
  # Check 1: are they all square matrices?
  n_elements<- length(Aobj[1,])
  if (n_elements%%1!=0){
    stop("Your matrices do not appear to be square.", call.=FALSE)
  }
  # Check 2: Are all elements non-negative?
  if (sum(Aobj<0)>0){
    stop("Population matrices must be strictly non-negative.", call.=FALSE)
  }

  # A few checks that we have to run after reshaping each matrix:
  for (i in 1:dim(Aobj)[1]){
    this_matrix<- reMat(Aobj, i)

    # Check 3: Is the matrix ergodic?
    if (popdemo::isErgodic(this_matrix)==FALSE){
      stop(paste("Matrix", i, "is not ergodic."), call.=FALSE)
    }

    # Check 4: is the matrix irreducible?
    if (popdemo::isIrreducible(this_matrix)==FALSE){
      stop(paste("Matrix", i, "is reducible."), call.=FALSE)
    }
    # Check 5: is the matrix primitive?
    if (popdemo::isPrimitive(this_matrix)==FALSE){
      stop(paste("Matrix", i, "is imprimitive."), call.=FALSE)
    }

    # Check 5: Check survival sums are not greater than 1.
    # assume that the first row is fecundity and all other entries represent survival/growth
    survival<- this_matrix[-1,]
    # if this_matrix is 2x2, the survival will be a vector. need to treat that differently than a matrix.
    if (is.null(dim(survival))){
      surv_sums<- survival
    } else { # if survival is still a vector, then we take the column sums:
      surv_sums<- apply(survival, 2, FUN=sum, na.rm=TRUE)
    }
    if (sum(surv_sums>1)>0){
      warning(paste("Survival terms may be misspecified. Matrix", i, "has column sums greater than 1 (assuming row 1 and only row 1 represents fecundity)."),
              call.=FALSE)
    }
  }
}
