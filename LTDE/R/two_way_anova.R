
rm(list = ls())
getwd()
setwd("~/github/Task2/LTDE/")

## Load Data
obs <- read.csv("data/Final/MAPGD_Evol_iRep.txt", 
                header = TRUE, stringsAsFactors = FALSE, sep = '\t')


int <- aov(obs$T_D ~ obs$model_type+obs$iRep+  obs$model_type:obs$iRep)
summary(int)


summary(aov(newdata$iRep ~ newdata$Genus+newdata$T_D))


int1 <- lm(T_D~evolvability+Genus + evolvability:Genus, data=newdata)
summary(int1)


###### Tajima's D vs evol, controlling for slope 
quads <- obs[ which(obs$model_type==1), ]
int1 <- lm(T_D~decay  , data=quads)

int2 <- lm(T_D~decay + evolvability , data=quads)
summary(int1)


AIC(int1, int2)
