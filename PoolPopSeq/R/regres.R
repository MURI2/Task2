rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")

package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}

# figure out how to do this with time series data.....

df <- read.table("data/gene_by_sample_long/gene_by_sample_long_P.txt", sep = "\t", header = TRUE)

df_sub <- select(df, LocusTag, Line, Muts)
gene_by_pop <- spread(df_sub, Line, Muts)
pop_by_gene <- t(subset(gene_by_pop, select = -c(LocusTag) ))

##Go through each row and determine if a value is zero
#row_sub <- apply(pop_by_gene, 1, function(row) all(row !=0 ))
##Subset as usual
#pop_by_gene.no0 <- pop_by_gene[row_sub,]

pop_by_gene.no0 <- pop_by_gene[rowSums(pop_by_gene[,-1]) != 0,]
pop_by_gene.db <- vegdist(pop_by_gene.no0, method = "bray", upper = TRUE, diag = TRUE)
fish.pcoa <- cmdscale(pop_by_gene.db, eig = TRUE, k = 3) 


explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)
# Initiate Plot
plot(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

pmatch("frequency_L0", c(rownames(fish.pcoa$points)))

cols <- c()
treats <- c()
for (x in rownames(fish.pcoa$points)){
  if (grepl("frequency_L0", x)){
    treats <- c(treats, "1")
    cols <- c(cols, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats <- c(treats, "10")
    cols <- c(cols, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats <- c(treats, "100")
    cols <- c(cols, "#FF6347")
  }
}


times <- c()
for (x in rownames(fish.pcoa$points)){
  if (grepl("D100", x)){
    times <- c(times, "100")
  } else if ( grepl("D200", x)) {
    times <- c(times, "200")
  } else if (grepl("D300", x)) {
    times <- c(times, "300")
  }
}


v <- c('a','b','c','e')

grepl("frequency_L0", rownames(fish.pcoa$points)[1])


# Add Points & Labels
cols_t1 <- c('#87CEEB', '#87CEEB', '#87CEEB', '#87CEEB', 
             '#FFA500',  '#FFA500',  '#FFA500', '#FFA500',
             '#FF6347', '#FF6347', '#FF6347', '#FF6347')
#plot
points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
       pch = 19, cex = 3, bg = "gray", col = cols_t1)

text(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], labels=times, cex= 0.7)

# Create "Factors" vector
quality <- c(rep("1", 4), rep("10", 4), rep("100", 4))
ordiellipse(fish.pcoa, treats, conf = 0.95)



#  Identifying and Visualizing Influential genes in PCoA Basic ordination plots 
fishREL <- pop_by_gene.no0
for(i in 1:nrow(pop_by_gene.no0)){
  fishREL[i, ] = pop_by_gene.no0[i, ] / sum(pop_by_gene.no0[i, ])
} 

# Now, we use this information to calculate and add species scores
fish.pcoa <- add.spec.scores(fish.pcoa,fishREL,method = "pcoa.scores")
text(fish.pcoa$cproj[ ,1], fish.pcoa$cproj[ ,2], 
     labels = row.names(fish.pcoa$cproj), col = "black")

spe.corr <- add.spec.scores(fish.pcoa, fishREL, method = "cor.scores")$cproj
corrcut  <- 0.7       # user defined cutoff
imp.spp  <- spe.corr[abs(spe.corr[, 1]) >= corrcut | abs(spe.corr[, 2]) >= corrcut, ]

# Permutation Test for Species Abundances Across Axes
fit <- envfit(fish.pcoa, fishREL, perm = 999)

which(fit$vectors$pvals < 0.05)
# 305 310 327 398 444

names <- gene_by_pop[,1]
names[305]#DR_2186
names[310]#DR_2226
names[327]#DR_2326
names[398]#DR_A0233 oxidoreductase, iron-sulfur subunit
names[444]#DR_C0029 transposase

pop_by_gene.no0[,305]
pop_by_gene.no0[,310]
pop_by_gene.no0[,327]
pop_by_gene.no0[,398]
pop_by_gene.no0[,444]



# Run PERMANOVA with adonis function
adonis(pop_by_gene.no0 ~ quality, method = "bray", permutations = 999)


# Pseudomonas


df <- read.table("data/gene_by_sample_long/gene_by_sample_long_P.txt", sep = "\t", header = TRUE)

df_sub <- select(df, LocusTag, Line, Muts)
gene_by_pop <- spread(df_sub, Line, Muts)
pop_by_gene <- t(subset(gene_by_pop, select = -c(LocusTag) ))

##Go through each row and determine if a value is zero
#row_sub <- apply(pop_by_gene, 1, function(row) all(row !=0 ))
##Subset as usual
#pop_by_gene.no0 <- pop_by_gene[row_sub,]

pop_by_gene.no0 <- pop_by_gene[rowSums(pop_by_gene[,-1]) != 0,]
pop_by_gene.db <- vegdist(pop_by_gene.no0, method = "bray", upper = TRUE, diag = TRUE)
fish.pcoa <- cmdscale(pop_by_gene.db, eig = TRUE, k = 3) 


explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)
# Initiate Plot
plot(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols_t1 <- c('#87CEEB', '#87CEEB', '#87CEEB', '#87CEEB', '#87CEEB',
             '#FFA500',  '#FFA500',  '#FFA500', '#FFA500', '#FFA500',
             '#FF6347', '#FF6347', '#FF6347', '#FF6347')
#plot
points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
       pch = 19, cex = 3, bg = "gray", col = cols_t1)
# Create "Factors" vector
quality <- c(rep("1", 5), rep("10", 5), rep("100", 4))
ordiellipse(fish.pcoa, quality, conf = 0.95)



#  Identifying and Visualizing Influential genes in PCoA Basic ordination plots 
fishREL <- pop_by_gene.no0
for(i in 1:nrow(pop_by_gene.no0)){
  fishREL[i, ] = pop_by_gene.no0[i, ] / sum(pop_by_gene.no0[i, ])
} 

# Now, we use this information to calculate and add species scores
fish.pcoa <- add.spec.scores(fish.pcoa,fishREL,method = "pcoa.scores")
text(fish.pcoa$cproj[ ,1], fish.pcoa$cproj[ ,2], 
     labels = row.names(fish.pcoa$cproj), col = "black")

spe.corr <- add.spec.scores(fish.pcoa, fishREL, method = "cor.scores")$cproj
corrcut  <- 0.7       # user defined cutoff
imp.spp  <- spe.corr[abs(spe.corr[, 1]) >= corrcut | abs(spe.corr[, 2]) >= corrcut, ]

# Permutation Test for Species Abundances Across Axes
fit <- envfit(fish.pcoa, fishREL, perm = 999)

which(fit$vectors$pvals < 0.05)
# 305 310 327 398 444

names <- gene_by_pop[,1]
names[36]#G_01240
names[97]#G_03466
names[108]#G_03893

pop_by_gene.no0[,36]
pop_by_gene.no0[,97]
pop_by_gene.no0[,108]


# Run PERMANOVA with adonis function
adonis(pop_by_gene.no0 ~ quality, method = "bray", permutations = 999)





####### J
df <- read.table("data/gene_by_sample_long/gene_by_sample_long_F.txt", sep = "\t", header = TRUE)

df_sub <- select(df, LocusTag, Line, Muts)
gene_by_pop <- spread(df_sub, Line, Muts)
pop_by_gene <- t(subset(gene_by_pop, select = -c(LocusTag) ))

##Go through each row and determine if a value is zero
#row_sub <- apply(pop_by_gene, 1, function(row) all(row !=0 ))
##Subset as usual
#pop_by_gene.no0 <- pop_by_gene[row_sub,]

pop_by_gene.no0 <- pop_by_gene[rowSums(pop_by_gene[,-1]) != 0,]
pop_by_gene.db <- vegdist(pop_by_gene.no0, method = "bray", upper = TRUE, diag = TRUE)
fish.pcoa <- cmdscale(pop_by_gene.db, eig = TRUE, k = 3) 


explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)
# Initiate Plot
plot(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols_t1 <- c('#87CEEB', '#87CEEB', '#87CEEB', '#87CEEB', '#87CEEB',
             '#FFA500',  '#FFA500', '#FFA500', 
             '#FF6347', '#FF6347', '#FF6347', '#FF6347', '#FF6347')
#plot
points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
       pch = 19, cex = 3, bg = "gray", col = cols_t1)
# Create "Factors" vector
quality <- c(rep("1", 5), rep("10", 3), rep("100", 5))
ordiellipse(fish.pcoa, quality, conf = 0.95)



#  Identifying and Visualizing Influential genes in PCoA Basic ordination plots 
fishREL <- pop_by_gene.no0
for(i in 1:nrow(pop_by_gene.no0)){
  fishREL[i, ] = pop_by_gene.no0[i, ] / sum(pop_by_gene.no0[i, ])
} 

# Now, we use this information to calculate and add species scores
fish.pcoa <- add.spec.scores(fish.pcoa,fishREL,method = "pcoa.scores")
text(fish.pcoa$cproj[ ,1], fish.pcoa$cproj[ ,2], 
     labels = row.names(fish.pcoa$cproj), col = "black")

spe.corr <- add.spec.scores(fish.pcoa, fishREL, method = "cor.scores")$cproj
corrcut  <- 0.7       # user defined cutoff
imp.spp  <- spe.corr[abs(spe.corr[, 1]) >= corrcut | abs(spe.corr[, 2]) >= corrcut, ]

# Permutation Test for Species Abundances Across Axes
fit <- envfit(fish.pcoa, fishREL, perm = 999)

which(fit$vectors$pvals < 0.05)
# 305 310 327 398 444

names <- gene_by_pop[,1]
names[36]#G_01240
names[97]#G_03466
names[108]#G_03893

pop_by_gene.no0[,36]
pop_by_gene.no0[,97]
pop_by_gene.no0[,108]


# Run PERMANOVA with adonis function
adonis(pop_by_gene.no0 ~ quality, method = "bray", permutations = 999)


