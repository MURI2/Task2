rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")

#strains <- c('B', 'C', 'D', 'F', 'J', 'P', 'S')
strains <- c('P')
strains.vec <- vector(mode="list", length=6)
names(strains.vec) <- c('B', 'C', 'D', 'F', 'J', 'P')
strains.vec[[1]] <- 'Bacillus'; strains.vec[[2]] <- 'Caulobacter'
strains.vec[[3]] <- 'Deinococcus'; strains.vec[[4]] <- 'Pedobacter'; 
strains.vec[[5]] <- 'Janthinobacterium';strains.vec[[6]] <- 'Pseudomonas'
package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR', 'indicspecies')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}

