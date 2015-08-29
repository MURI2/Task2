rm(list=ls())
getwd()
setwd('/Users/WRShoemaker/MURI/MURI_growth_curves/GC_dilution/')
getwd()

library(ggplot2)
library(xlsx)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(Rmisc)

ODs <- read.xlsx("GC_dilution_OD.xlsx", 2)

melted <- melt(ODs, id.vars="Time")

tgc <- summarySE(melted, measurevar="value", groupvars=c("variable", "Time"))

ggplot(tgc, aes(x=Time, y=value, colour=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point() +  
  xlim(0, 28) + 
  ylim(0, 6)

ggplot(data = tolerance, aes(x = time, y = tolerance)) + geom_line() +
  facet_wrap(~id)

ggplot(tgc, aes(x=Time, y=value, colour=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point() + 
  facet_wrap(~variable)



ggplot(tgc, aes(x=Time, y=value, colour=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point() + 
  facet_wrap(~variable) + 
  xlim(0, 28) + 
  ylim(0, 6)

#myColors <- brewer.pal(24,"Set1")
#names(myColors) <- levels(melted$variable)
#colScale <- scale_colour_manual(name = "variable",values = myColors)

#P <- ggplot(data=melted, aes(x=Time, y=value, group=variable, color = variable)) + geom_line() + 
#  geom_point(size = 3) + 
#  theme_gray()
#P1 <-P + colScale



#ggplot(data=melted, aes(x=Time, y=value, group=variable)) + 
#  geom_bar(width=1)+scale_y_continuous(expand = c(0,0))+ theme(axis.text.x=element_text(angle=90))
