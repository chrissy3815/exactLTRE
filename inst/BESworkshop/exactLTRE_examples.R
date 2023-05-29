#################################################################
##
## Code for tutorial/workshop on exactLTRE package
##
#################################################################

# Other libraries needed for this tutorial code:
# If you don't have them, install them from CRAN
library(Rcompadre)

# Install exactLTRE from CRAN:
install.packages("exactLTRE")

# Load the exactLTRE library:
library(exactLTRE)

## Examples from the package documentation:

# Build some example matrices:
A1<- matrix(data=c(0,0.8,0, 0,0,0.7, 5,0,0.2), nrow=3, ncol=3)
A2<- matrix(data=c(0,0.9,0, 0,0,0.5, 4,0,0.3), nrow=3, ncol=3)
A3<- matrix(data=c(0,0.4,0, 0,0,0.6, 6,0,0.25), nrow=3, ncol=3)

# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- classicalLTRE(list(A1,A2), method='fixed')
cont_diff; # matrix of contributions
sum(cont_diff); # sum of approximate contributions
lamDiff(list(A1,A2)); # true difference in lambda

# contributions to the variance of lambda
cont_var<- classicalLTRE(list(A1,A2,A3), method='random')
round(cont_var,digits=5); # matrix of contributions
sum(cont_var); # sum of contributions
lamVar(list(A1,A2,A3)); # true variance in lambda

# Exact LTRE:
# contributions to the difference in lambda
cont_diff<- exactLTRE(list(A1,A2), method='fixed')
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions
lamDiff(list(A1,A2)); # true difference in lambda

# contributions to the variance of lambda
cont_var<- exactLTRE(list(A1,A2,A3), method='random')
cbind(as.character(cont_var$varying.indices.list),round(cont_var$epsilons,digits=5));
sum(cont_var$epsilons); # sum of contributions
lamVar(list(A1,A2,A3)); # true variance in lambda

# Generation Time:
F1<- matrix(0, nrow=3, ncol=3)
F1[1,3]<- A1[1,3]
#F1 is all zeros, except the upper right corner which matches A1 for adult fertility
gen_time <- generation_time(A1, F1)

# R0, the expected lifetime reproductive output for an individual
R0<- r_nought(A1, F1)

# expected lifespan
U1<- A1
U1[1,3]<- 0
# the upper right corner represents adult fertility in this model. U1, the
# survival matrix, contains all the transitions *except* for fertility.
eta<- lifespan(U1, all_ages=TRUE)
eta_1<- lifespan(U1, all_ages=FALSE) # eta_1 should match the first entry of eta

#################################################################
## Examples with COMPADRE
#################################################################

# Load comadre and compadre databases:
comadre<- cdb_fetch("comadre")
compadre<- cdb_fetch("compadre")

## Note: if you can't install Rcompadre because of OS incompatibility, then
# uncomment this line, make sure the file path is correct, and run it to load
# the matrices and metadata for the following few examples.
# This file contains the following four objects:
# geocrinia, geocrinia_mats, spermophilus, spermophilus_mats,
load('inst/BESworkshop/example_matrices.Rdata')

# We'll first work through a simple example of fixed design LTREs with white-bellied frogs.
geocrinia<- cdb_flatten(comadre[comadre$SpeciesAuthor=="Geocrinia_alba",])
# cdb_flatten extracts the metadata only, without the `mat` list object.

# You can look through the metadata about these matrices in the object "geocrinia"
# For a full list of the available metadata, try:
names(geocrinia)
# You can learn more about the variables available in metadata here:
# https://jonesor.github.io/CompadreGuides/user-guide.html#variables-in-metadata

# Here's an example of pulling out some metadata:
geocrinia[,c("SpeciesAccepted", "Country", "MatrixPopulation", "MatrixTreatment", "MatrixID")]

# pull out the two matrices of interest:
geocrinia_mats<- c(matA(comadre[comadre$MatrixID==248239,]),matA(comadre[comadre$MatrixID==248238, ]))
lapply(geocrinia_mats, eigen, only.values=TRUE)

# lambda for Forest Grove South population:
eigen(geocrinia_mats[[1]], only.values=TRUE)$values[1]
# lambda for Bruce Road population:
eigen(geocrinia_mats[[2]], only.values=TRUE)$values[1]
# Forest Grove South population is growing, Bruce Rd population is shrinking.
# difference in lambda (Bruce Road - Forest Grove South):
eigen(geocrinia_mats[[2]], only.values=TRUE)$values[1]-eigen(geocrinia_mats[[1]], only.values=TRUE)$values[1]

# How do the matrix elements differ?
differences<- as.vector(geocrinia_mats[[2]] - geocrinia_mats[[1]])
barplot(differences, main='Differences', names.arg=c('', 'sJ', 'f', 'sA'))

# Which matrix elements are driving the difference?
# Evaluate a fixed symmetric design LTRE:
result<- exactLTRE(geocrinia_mats, method='fixed', fixed.directional = FALSE)
result$varying.indices.list
barplot(t(result$epsilons[2:length(result$epsilons)]), main='Contributions',
        names.arg=c("sJ", "f","sJ, f", "sA", "sJ, sA", "f, sA", "sJ, f, sA"),
        las=2)

#################################################################
## Ground squirrel example:
#################################################################

## Spermophilus armatus, ground squirrels.
spermophilus<- cdb_flatten(comadre[comadre$MatrixID %in% c(249840,249844),])
# cdb_flatten extracts the metadata only.

# pull out the matrices:
spermophilus_mats<- c(matA(comadre[comadre$MatrixID==249840]), matA(comadre[comadre$MatrixID==249844]))

# lambda for unmanipulated population:
eigen(spermophilus_mats[[1]], only.values=TRUE)$values[1]
# lambda for density reduction population:
eigen(spermophilus_mats[[2]], only.values=TRUE)$values[1]
# These populations had very similar population growth rates!
# difference in lambda (treatment - reference):
eigen(spermophilus_mats[[2]], only.values=TRUE)$values[1]-eigen(spermophilus_mats[[1]], only.values=TRUE)$values[1]

# How do the matrix elements differ?
differences<- as.vector(spermophilus_mats[[2]] - spermophilus_mats[[1]])
differences<- differences[abs(differences)>0]
barplot(differences, main='Differences', names.arg=c("F1", "P1", "F2", "P2", "F3", "P3"))

# Should we do a symmetric or directional LTRE?
spermophilus[,c("MatrixID", "MatrixPopulation", "MatrixTreatment")]

# We want to compare a treatment with a control, so this should be DIRECTIONAL
spermophilus_dir<- exactLTRE(spermophilus_mats, method='fixed', maxint=3, fixed.directional = TRUE)
xlabels<- c("F1", "P1", "F2", "P2", "F3", "P3",
            "F1+P1", "F1+F2", "F1+P2", "F1+F3", "F1+P3", "P1+F2", "P1+P2",
            "P1+F3", "P1+P3", "F2+P2", "F2+F3", "F2+P3", "P2+F3", "P2+P3",
            "F3+P3", "F1+P1+F2", "F1+P1+P2", "F1+P1+F3", "F1+P1+P3",
            "F1+F2+P2", "F1+F2+F3", "F1+F2+P3", "F1+P2+F3", "F1+P2+P3",
            "F1+F3+P3", "P1+F2+P2", "P1+F2+F3", "P1+F2+P3", "P1+P2+F3",
            "P1+P2+P3", "P1+F3+P3", "F2+P2+F3", "F2+P2+P3", "F2+F3+P3", "P2+F3+P3")
barplot(spermophilus_dir$epsilons[-1], beside = TRUE, names.arg = xlabels, las=2)

# Let's also see how the interpretation changes when we switch between symmetric and directional.
spermophilus_symm<- exactLTRE(spermophilus_mats, method='fixed', maxint=3, fixed.directional = FALSE)

# plot the comparison:
toplot<- rbind(spermophilus_symm$epsilons, spermophilus_dir$epsilons)
toplot<- toplot[,-1]
xlabels<- c("F1", "P1", "F2", "P2", "F3", "P3",
            "F1+P1", "F1+F2", "F1+P2", "F1+F3", "F1+P3", "P1+F2", "P1+P2",
            "P1+F3", "P1+P3", "F2+P2", "F2+F3", "F2+P3", "P2+F3", "P2+P3",
            "F3+P3", "F1+P1+F2", "F1+P1+P2", "F1+P1+F3", "F1+P1+P3",
            "F1+F2+P2", "F1+F2+F3", "F1+F2+P3", "F1+P2+F3", "F1+P2+P3",
            "F1+F3+P3", "P1+F2+P2", "P1+F2+F3", "P1+F2+P3", "P1+P2+F3",
            "P1+P2+P3", "P1+F3+P3", "F2+P2+F3", "F2+P2+P3", "F2+F3+P3", "P2+F3+P3")
barplot(toplot, beside = TRUE, names.arg = xlabels, las=2,
        legend=c("Symmetric", "Directional"), args.legend = list(x="topright"),)

#################################################################
## Random design LTRE example
#################################################################

## Alliaria petiolata, garlic mustard
alliaria<- cdb_flatten(compadre[compadre$MatrixID %in% c(241484,241485, 241486),])
# Note: cdb_flatten extracts the metadata only.
alliaria[,c("MatrixID", "MatrixPopulation", "MatrixStartYear", "MatrixEndYear")]
# We are going to look at the variance in lambda across 3 years at the 'Ives Road' population

# pull out the matrices:
alliaria_mats<- matA(compadre[compadre$MatrixID %in% c(241484,241485, 241486)])
# Note: for random design, order does not matter.

# lambda for first year:
eigen(alliaria_mats[[1]], only.values=TRUE)$values[1]
# lambda for second year:
eigen(alliaria_mats[[2]], only.values=TRUE)$values[1]
# lambda for third year:
Re(eigen(alliaria_mats[[3]], only.values=TRUE)$values[1])
# These years had very similar population growth rates!
# variance in lambda
lamVar(alliaria_mats)

# How do the matrix elements vary?
elems<- vector()
elems<- sapply(alliaria_mats, function(x){cbind(elems, as.vector(x))})
variances<- apply(elems, MARGIN = 1, variance_complete)
variances<- variances[variances>0]
# Plot the variances:
barplot(variances, main='Variances', names.arg=c("Ps", "Gr", "Gf", "Fs", "Fr"))

# Run the exact LTRE:
alliaria_exact<- exactLTRE(alliaria_mats, method='random', maxint = 3)
# Note: that warning is coming up because reproductive individuals can produce two types of offspring!
matF(compadre[compadre$MatrixID==241484])

# Look at some of the results:
alliaria_exact$indices.varying
alliaria_exact$varying.indices.list

# Plot the results of the exact LTRE
xlabels<- c("Ps", "Gr", "Gf", "Fs", "Fr", "Ps+Gr", "Ps+Gf", "Ps+Fs", "Ps+Fr",
            "Gr+Gf", "Gr+Fs", "Gr+Fr", "Gf+Fs", "Gf+Fr", "Fs+Fr", "Ps+Gr+Gf",
            "Ps+Gr+Fs", "Ps+Gr+Fr", "Ps+Gf+Fs", "Ps+Gf+Fr", "Ps+Fs+Fr",
            "Gr+Gf+Fs", "Gr+Gf+Fr", "Gr+Fs+Fr", "Gf+Fs+Fr")
barplot(alliaria_exact$epsilons[-1], beside = TRUE, names.arg = xlabels, las=2)

