# Functions for exactLTRE analyses:

## Main LTRE functions ---------------------------------------------------------
# a wrapper function for standard LTRE:
#' Approximate LTRE analysis
#'
#' @param Aobj An object containing all the population projection matrices to be
#'  included in the analysis. It should either be a list, or a matrix where each
#'   row is the column-wise vectorization of a matrix.
#'
#'  For fixed design, exactly 2 matrices must be provided, ordered as
#'  `[reference matrix, test matrix`]. For random design, any set of 2 or more
#'  matrices can be provided. The set of matrices passed in must all have the
#'  same dimensions.
#'
#' @param method Either "random" or "fixed." The default behavior is "random."
#' See details for more information.
#'
#' @return A matrix of contributions to variance (random design) or difference
#' (fixed design) in lambda. Lambda is the asymptotic population growth rate,
#' defined as the largest eigenvalue of the population projection matrix.
#'
#' @export
#'
#'@details  Lambda is the asymptotic population growth rate, defined as the
#'  largest eigenvalue of the population projection matrix. A fixed design LTRE
#'  decomposes the difference in lambda due to differences at each position of
#'  the matrices. For a fixed design LTRE, exactly 2 matrices must be provided,
#'  ordered as `[reference matrix, test matrix`]. The matrix of contributions
#'  returned from a standard method fixed design LTRE will have the same shape
#'  as the provided matrices.
#'
#'  A random design LTRE decomposes the variance in lambda due to variance and
#'  covariance in the entries at each position in the matrices. For a random
#'  design LTRE, at least 2 matrices must be provided. The matrix of
#'  contributions returned from a standard method random design LTRE will
#'  include both first-order terms (due to variance) and interaction terms (due
#'  to covariance). Therefore, if the provided matrix is 3x3, the matrix of
#'  contributions will be 9x9 (the size of the variance-covariance matrix is the
#'  square of the size of the original matrix). The contributions of variances
#'  are found on the diagonal of the contribution matrix, and the contributions
#'  of covariances are symmetric. So the contribution of covariance between two
#'  vital rate parameters is the sum of the two corresponding off-diagonal terms.
#'
#'  The standard methods of LTRE analysis are based on approximation methods.
#'  The equations and descriptions can be found in Caswell's 2001 textbook.
#'
#' @examples
#' A1<- matrix(data=c(0,0.8,0, 0,0,0.7, 5,0,0.2), nrow=3, ncol=3)
#' A2<- matrix(data=c(0,0.9,0, 0,0,0.5, 4,0,0.3), nrow=3, ncol=3)
#' A3<- matrix(data=c(0,0.4,0, 0,0,0.6, 6,0,0.25), nrow=3, ncol=3)
#' cont_diff<- standardLTRE(list(A1,A2), method='fixed') # contributions to the difference in lambda
#' cont_var<- standardLTRE(list(A1,A2,A3), method='random') # contributions to the variance of lambda
standardLTRE<- function(Aobj, method="random"){

  # if Aobj is passed in as a list, collapse into the row-format:
  if (is.list(Aobj)){
    Aobj<- collapse_mat_list(Aobj)
  }

  if (method=="random"){
    output<- standardLTRE_random(Aobj)
  } else if (method=="fixed"){
    # For fixed LTRE, it is important that Aobj is 2 rows. Row 1 contains vec(Aref), and Row 2 contains vec(Atest)
    Aref<- reMat(Aobj,1)
    Atest<- reMat(Aobj,2)
    output<- standardLTRE_fixed(Aref, Atest)
  }
  return(output)
}

# A wrapper function for exact LTRE:
exactLTRE<- function(Aobj, method="random", maxint="all"){

  # if Aobj is passed in as a list, collapse into the row-format:
  if (is.list(Aobj)){
    Aobj<- collapse_mat_list(Aobj)
  }

  if (method=="random"){
    output<- exactLTRE_random(Aobj, maxint)
  } else if (method=="fixed"){
    output<- exactLTRE_fixed(Aobj, maxint)
  }
  return(output)
}

# A function for standard LTRE, fixed design:
standardLTRE_fixed<- function(Aref, Atest){

  # run some matrix checks and return warnings as needed:
  run_matrix_checks(rbind(as.vector(Aref), as.vector(Atest)))

  Amean<- (Atest+Aref)/2 # define the mean matrix
  wmean<- Re(eigen(Amean)$vectors[,1]) # right eigenvector of the mean matrix
  vmean<- Re(eigen(t(Amean))$vectors[,1]) # left eigenvector of the mean matrix
  sensmat<- vmean%*%t(wmean)/as.vector(vmean%*%wmean)
  diffmat<- Atest-Aref
  C_m<- diffmat*sensmat
  return(C_m)
}
# A function for standard LTRE, random design:
standardLTRE_random<- function(Aobj){
  # each row of Aobj should be the elements of an individual population matrix, collapsed column-wise
  # if Aobj is passed in as a list, collapse into the row-format:
  if (is.list(Aobj)){
    Aobj<- collapse_mat_list(Aobj)
  }

  # run some matrix checks and return warnings as needed:
  run_matrix_checks(Aobj)

  # covariance matrix for the parameters in Aobj:
  Cmat<- mycov(Aobj)

  # calculate the mean parameter values and reshape into a square matrix:
  Amean<- matrix(apply(Aobj,2,mean),sqrt(dim(Aobj)[2]),sqrt(dim(Aobj)[2]))

  # calculate the sensitivity, evaluated at the mean matrix:
  wmean<- Re(eigen(Amean)$vectors[,1]) # right eigenvector of the mean matrix
  vmean<- Re(eigen(t(Amean))$vectors[,1]) # left eigenvector of the mean matrix
  sensmat<- vmean%*%t(wmean)/as.vector(vmean%*%wmean)

  # calculate contributions according to: Cov(aij,akl)*sij*skl
  s<- as.vector(sensmat) # this is vec
  contmat<- Cmat*(s%*%t(s)) # Cmat *hadamard* vec(s)*t(vec(s))
  return(contmat)
}


# A function for exact LTRE, random design:
exactLTRE_random<- function(Aobj, maxint="all"){
  # each row of Aobj should be the elements of an individual population matrix, collapsed column-wise
  # if Aobj is passed in as a list, collapse into the row-format:
  if (is.list(Aobj)){
    Aobj<- collapse_mat_list(Aobj)
  }

  # run some matrix checks and return warnings as needed:
  run_matrix_checks(Aobj)

  #find the indices of the columns in Aobj that have non-zero variance. In other
  # words, which matrix entries vary among populations?
  ind_vary<- which(apply(Aobj, MARGIN=2, FUN=Vc)>0)

  # calculate the responses, the vector nu: aka the variance in lambda,
  # given each possible combination of parameters varying or held at their mean:
  responses<- calc_matrix_responses(Aobj, ind_vary, FUN=lamVar, maxint)

  # Calculate the effects or "epsilons," depending on maxint.
  if (maxint=="all"){ # use Poelwijk approach

    Gmatrix<- make.Gmatrix(length(ind_vary)) # order of Gmatrix is the number of parameters that vary
    # calculate the epsilon values:
    epsilons<- Gmatrix%*%responses$nus

  } else { # directly calculate the epsilons:

    epsilons<- calc_matrix_effects(responses$nus, responses$list_ind_vary)

  }

  # make the data to return:
  output<- list(indices.varying = ind_vary,
                varying.indices.list = responses$list_ind_vary,
                epsilons = epsilons)
  return(output)

}
# A function for exact LTRE, fixed design
exactLTRE_fixed<- function(Aobj, maxint="all"){
  # each row of Aobj should be the elements of an individual population matrix, collapsed column-wise
  # It is important that Aobj is 2 rows. Row 1 contains vec(Aref), and Row 2 contains vec(Atest)

  # if Aobj is passed in as a list, collapse into the row-format:
  if (is.list(Aobj)){
    Aobj<- collapse_mat_list(Aobj)
  }

  # run some matrix checks and return warnings as needed:
  run_matrix_checks(Aobj)

  # get rid of column and row names, just in case:
  colnames(Aobj)<- NULL; rownames(Aobj)<- NULL

  #find the indices of the columns in Aobj that have non-zero difference. In other
  # words, which matrix entries vary between populations?
  ind_vary<- which(abs(Aobj[1,]-Aobj[2,])>0)

  # checks on maxint:
  if (maxint!="all"){
    if (!is.numeric(maxint)){
      stop("Input variable maxint, the maximum interaction order to calculate, must either be ''all'' or a numeric variable.")
    } else if (maxint>=length(ind_vary)){
      warning("You have requested an interaction order that is greater than or equal to the number of parameters that vary, defaulting to maxint='all'.")
      maxint="all"
    }
  }

  # calculate the matrix responses:
  responses<- calc_matrix_responses(Aobj, ind_vary, FUN=lamDiff, maxint)

  # Calculate the effects or "epsilons," depending on maxint.
  if (maxint=="all"){ # use Poelwijk approach

    Gmatrix<- make.Gmatrix(length(ind_vary)) # order of Gmatrix is the number of parameters that vary
    # calculate the epsilon values:
    epsilons<- Gmatrix%*%responses$nus

  } else { # directly calculate the epsilons:

    epsilons<- calc_matrix_effects(responses$nus, responses$list_ind_vary)

  }

  # make the data to return:
  output<- list(indices.varying = ind_vary,
                varying.indices.list = responses$list_ind_vary,
                epsilons = epsilons)
  return(output)

}

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

## Auxiliary functions ---------------------------------------------------------
# A function for variance assuming a complete sample:
Vc <- function(x) {
  xx<- x[!is.na(x)]
  mean((xx - mean(xx))^2)
}

# a function for the covariance object:
mycov<- function(Aobj){
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
  return(Vc(lambdas))
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

# function for generating "my" Gmatrix from list_ind_vary:
my.Gmatrix<- function(list_ind_vary){
  matdim<- length(list_ind_vary)
  Gmatrix<- diag(nrow=matdim, ncol=matdim) # start with a diagonal matrix

  for (i in 2:matdim){# go row-wise, starting with row 2
    ind_i<- list_ind_vary[[i]] # which indices were varying for this row?
    order_i<- length(list_ind_vary[[i]]) # what order interaction is represented by row i?

    for (j in 1:(i-1)){# we only want to check the sub-diagonal entries

      if (j==1){
        # the first column is the multiplicative factor for the baseline/constant term
        # (note, this is not relevant for LTRE because nu_0=eps_0=0 ALWAYS in LTRE)
        Gmatrix[i,j]<- (-1)^(order_i)
        next
      }
      # which indices were varying for the entry of nu corresponding to this column?
      ind_j<- list_ind_vary[[j]]
      order_j<- length(list_ind_vary[[j]])

      if (sum(ind_j %in% ind_i)==length(ind_j)){ # this is the row-entries that will be non-zero
        if (order_i%%2 == order_j%%2){
          Gmatrix[i,j]<- 1
        } else {
        Gmatrix[i,j]<- (-1)
        }
      }
    }
  }
  return(Gmatrix)
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

# Calculate the fundamental matrix
fundamental_matrix<- function(Umat){
  # Umat contains all the transitions for individuals (so Amat-Fmat)

  # the fundamental matrix, Nmat, is inv(I-U))
  Nmat<- matrixcalc::matrix.inverse(diag(dim(Umat)[1])-Umat)
  return(Nmat)
}

