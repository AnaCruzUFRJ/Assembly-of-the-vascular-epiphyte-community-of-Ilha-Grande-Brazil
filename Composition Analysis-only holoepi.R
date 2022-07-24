###############################################
#                                             #
# Date: July-2022
# Name: Canonical analysis of principal coordinates
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description: analyze the variation in the composition of epiphyte species in the RAPELD plots of Ilha Grande, Brazil

# ONLY HOLOEPIPHYTES

## Dataset referring to the manuscript: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#                                           ################################################
#                                       #
##### load directory
setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/Manuscrito 2/Repositorio_Scripts")
################################################
#                                       #
# load packages
#
#############################################################
library(tidyverse) # data manipulation
library(tidyselect) # data manipulation
library (vegan)# for the hellinger function
library(dplyr) # data manipulation
library(ape) # to make pcoa
library(BiodiversityR) # to make pcoa
library(permute) # to do ca and cca
library(lattice) # ditto above
library(ggplot2) # graphics
library(ggthemes)# graphics
library(RColorBrewer)#graphics
###############################################################
################################################################################################################################### canonical analysis of principal coordinates#####################################################################################################################################################################################
#############
# Only CAP gives the flexibility to use any measure of distance, similarity or dissimilarity desired , such as the Bray-Curtis dissimilarity.# 
##############
setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/Manuscrito 2/Repositorio_Scripts")

# matriz de especies
especiescompleto<-read.table("MATRIZGERAL2HOLOEPI.csv",header=T,sep=";")

#View(especiescompleto)
especies<- especiescompleto %>% select(-Segmento, -Plots)
#View(especies)
head(especies)
str(especies)
rowSums(especies) # no line sum zero
# environmental matrix
ambientecompleto<- read.table("Dados_conjunto3.csv",header=T,sep=";", dec=",")
head(ambientecompleto)
#View(ambientecompleto)
ambiente <- ambientecompleto %>% select(Cobertura_dossel,Inclinacao, Densidade_For, DAP, HF, H, PROF_COPA,DIAM_COPA)
head(ambiente)
ambiente <- data.frame(scale(ambiente))
View(ambiente)
cap <- capscale(especies ~ ., data=ambiente, distance="bray", add=T, scale = T)
cap
summary(cap) #
plot(cap, choices=c(1, 2))
anova(cap)
anova(cap, by="axis")
anova(cap, by="terms")
# checking self-correlated variables #
vif.cca(cap)
# Remember that VIF values > 4 indicate that the variable is redundant and should be removed from the analysis.
ambiente$H <- NULL
ambiente$PROF_COPA <- NULL
ambiente$HF <- NULL
head(ambiente)
cap2 <- capscale(especies ~ ., data=ambiente, distance="bray", add=T, scale=T)

vif.cca(cap2)
cap2 # the total variation is 47.53 and the one explained by the variables is 7.35 or 15.46%
summary(cap2)
plot(cap2, choices=c(1, 2),display=c("species", "cn"))
anova(cap2)
anova(cap2, by="axis")
anova(cap2, by="terms")

# are significant but explain very little

RsquareAdj(cap2)

cap0 <- capscale(especies ~ 1, data=ambiente, distance="bray", add=T)
cap0
cap.step <- ordiR2step(cap0, cap2)
v11()
par(mfrow=c(1, 2))
plot(cap.step, choices=c(1, 2), display=c("sites", "cn"))
plot(cap.step, choices=c(1, 2), display=c("species", "cn"), ylim = c(-0.7, 0.7), xlim = c(0, 0.15))
plot(cap2, choices=c(1, 2))

