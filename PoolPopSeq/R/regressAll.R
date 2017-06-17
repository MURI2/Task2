rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")

strains <- c('B', 'C', 'D', 'F', 'J', 'P')
#strains <- c('P')

package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}

for (strain in strains) {
  #strain
  #tabele.path <- 
  sdata <- c("data/gene_by_sample_long/gene_by_sample_long_", strain, ".txt")
  #print(paste(sdata, collapse = ''))
  classes <- as.character(c("numeric", "character", "numeric", "numeric", "character", 
               "character", "numeric", "numeric", "character", "numeric"))
  df <- read.table(paste(sdata, collapse = ''), sep = "\t", header = TRUE, colClasses = classes)
  df_sub <- select(df, LocusTag, Line, Muts)
  gene_by_pop <- spread(df_sub, Line, Muts)
  pop_by_gene <- t(subset(gene_by_pop, select = -c(LocusTag) ))
  pop_by_gene.no0 <- pop_by_gene[rowSums(pop_by_gene[,-1]) != 0,]
  pop_by_gene.db <- vegdist(pop_by_gene.no0, method = "bray", upper = TRUE, diag = TRUE)
  fish.pcoa <- cmdscale(pop_by_gene.db, eig = TRUE, k = 3) 
  
  explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
  explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
  explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
  sum.eig <- sum(explainvar1, explainvar2, explainvar3)
  png(filename = paste(c("figs/pcoa/D100_CDS_", strain, ".png"), collapse = ''),
      width = 1200, height = 1200, res = 96*2)
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
  #plot
  points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
         pch = 19, cex = 3, bg = "gray", col = cols)
  ordiellipse(fish.pcoa, treats, conf = 0.95)
  dev.off()
  #dev.off()
  #graphics.off()
  #img.name <- paste(c("figs/pcoa/D100_CDS_", strain, ".txt"), collapse = '')
  #img <- readPNG("./figures/Figure15.png")
  #grid.raster(img)
}



