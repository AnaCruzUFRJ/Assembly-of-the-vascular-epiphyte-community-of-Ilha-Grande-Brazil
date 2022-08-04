
###############################################
#                                             #
# Date: 19-August-2021
# Name: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest 
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description:recognize, execute and interpret the accumulation curve, richness estimators and species rarefaction curve for data from Ilha Grande, RJ
# 
## Dataset referring to the manuscript: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#                                           ################################################
#                                       #


# install.packages(BiodiversityR)
# install.packages(vegan)

library(vegan)
library(BiodiversityR)
library(tidyverse)

# grafico da curva de acumulacao de especies para as parcelas na ilha grande  

############# curva de acumulacao da restinga 1

setwd("D:/Arquivos/Documents/ANA 2021/TESE/Boneca 2/Analises capitulo 2")


REST1<-read.table("composição por segmento_restinga_O1.csv",header=T,sep=";", dec=",")
View(REST1)

# curva de acumulacao de especies
acum<-specaccum(REST1)
acum
plot(acum, col="black",xlab="Sample units", ylab="Accumulated wealth of species", main="Plot 1")
# ou

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="REST1", xlim=c(0,30), ylim=c(0,70))

# estimadores de riqueza chao 1, Jackknife 1 e 2 e Bootstrap 

Chao22<-specpool(REST1)
Chao22 # anotar os valores dos estimadores e desvios padroes

Chao2<-poolaccum(REST1, permutations=100)
summary(Chao2, display = "chao")

plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", main="")

# restinga 2


REST2<-read.table("composição por segmento_restinga_O2.csv",header=T,sep=";", dec=",")

rarecurve(REST1, col = "firebrick", cex = 0.6, xlab = "Tamanho amostral",ylab="Espécies")
# curva de acumulacao
acum<-specaccum(REST2)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="REST2", xlim=c(0,30), ylim=c(0,70))

Chao22<-specpool(REST2)
Chao22 # anotar os valores dos estimadores e desvios padroes

Chao2<-poolaccum(REST2, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Restinga 2")

# parnaioca


Parna<-read.table("composição por segmento_parnaioca.csv",header=T,sep=";", dec=",")


Chao22<-specpool(Parna)
Chao22 # anotar os valores dos estimadores e desvios padroes

# curva de acumulacao
acum<-specaccum(Parna)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling affort", ylab="Accumulated species richness", col="black", main="DOLF", xlim=c(0,30), ylim=c(0,70))

#

Chao2<-poolaccum(Parna, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Parnaioca")

##### poco do soldado 

poco<-read.table("composição por segmento_poço.csv",header=T,sep=";", dec=",")

# curva de acumulacao
acum<-specaccum(poco)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="DOSF1", xlim=c(0,30), ylim=c(0,70))

Chao22<-specpool(poco)
Chao22 # anotar os valores dos estimadores e desvios padroes

Chao2<-poolaccum(poco, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Poco do Soldado")

### britador 

brit<-read.table("composição por segmento_britador.csv",header=T,sep=";", dec=",")

# curva de acumulacao
acum<-specaccum(brit)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="DOSF2", xlim=c(0,30), ylim=c(0,70))


# 
Chao22<-specpool(brit)
Chao22 # anotar os valores dos estimadores e desvios padroes

Chao2<-poolaccum(brit, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Britador")


# indice de diversidade 

# restinga 1
divers.all <- diversityresult(REST1, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(REST1, index = "Simpson", method="pooled") 
divers.all

# restinga 2
divers.all <- diversityresult(REST2, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(REST2, index = "Simpson", method="pooled") 
divers.all

# parnaioca
divers.all <- diversityresult(Parna, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(Parna, index = "Simpson", method="pooled") 
divers.all

# poco do soldado

divers.all <- diversityresult(poco, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(poco, index = "Simpson", method="pooled") 
divers.all

# britador


divers.all <- diversityresult(brit, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(brit, index = "Simpson", method="pooled") 
divers.all
