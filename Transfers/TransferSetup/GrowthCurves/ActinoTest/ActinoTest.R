rm(list=ls())
getwd()
setwd('/Users/WRShoemaker/MURI/MURI_growth_curves/ActinoTest///')
getwd()

library(ggplot2)
library(xlsx)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(Rmisc)

ODs <- read.xlsx("ActinoTestData.xlsx", 1)

melted <- melt(ODs, id.vars="Time")

tgc <- summarySE(melted, measurevar="value", groupvars=c("variable", "Time"))

ggplot(tgc, aes(x=Time, y=value, colour=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point() +  
  xlim(0, 48) + 
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
  xlim(0, 48) + 
  ylim(0, 8)
