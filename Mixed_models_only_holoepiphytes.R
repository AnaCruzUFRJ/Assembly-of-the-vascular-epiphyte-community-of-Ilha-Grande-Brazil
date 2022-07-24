
###############################################
#                                             #
# Date: 21-July-2022
# Name: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest 
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description: analysis of the effect of environmental and morphometric variations of trees on the richness and abundance of vascular epiphytes in 5 RAPELD plots of Ilha Grande, Brazil
# 
##### ONLY HOLOEPIPHYTES #### 


## Dataset referring to the manuscript: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#                                           ################################################
#                                       #
##### load directory
setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/Manuscrito 2/Repositorio_Scripts")
################################################
#                                       #
# load packages
library(tidyverse) # data manipulation
library(dplyr)# data manipulation
library(psych) # descriptive analysis
library(car) # normality test
library(MASS)# to boxcox
library(ggthemes)#graphics
library(vegan)
library(unholy)
library(AID)# to boxcox
library(broom) #graphics
library(ggplot2) #for plots
library(RColorBrewer)#graphics
library(RVAideMemoire) #to plot normality plots of residuals
library(lme4)# to LMM
library(bbmle) # we use funcao AICctab do package bbmle
library(rJava)
library(glmulti)# select better models
library(ggpubr)# join figures
library(hrbrthemes)# graphic
library(metaphor)
library(jtools)# output models
library(multcomp) # for tukey - glm post hoc test
################################################
#
dados1<-read.table("Dados_conjunto_sem_hemiep.csv",header=T,sep=";", dec=",")
dados2<-read.table("Dados_conjunto3.csv",header=T,sep=";", dec=",")
head(dados1)
head(dados2)
dados.scaled <- as.data.frame(scale(dados2)) #using standard deviation units because it has information in different measurement units
head(dados.scaled)
dados<-data.frame(dados1,dados.scaled) 
#View(dados)
head(dados)
str(dados)

dados$Parcela<-as.factor(dados$Parcela)
str(dados)

###### testing data normality

# Richness
shapiro.test(dados$S)
# there is no normality
mS <- lm(S ~ 1, data = dados)
summary(mS)
shapiro.test(resid(mS)) # testing normality of residuals
plotresid(mS, shapiro = T) # there is no normality
BOX_COX<-dados %>%
  do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(dados$Abd)
# there is no normality
mAbd <- lm(Abd ~ 1, data = dados)
summary(mAbd)
shapiro.test(resid(mAbd)) # testing normality of residuals
plotresid(mAbd, shapiro = T) # # there is no normality
##############################################################################################################
# 
# species richness analysis # 
#
##############################################################################################################################
# 
# GLM versus glmm
# GLM
head(dados)
Modelos.riqueza <- glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)

# no variables with high correlation
# Best model: S~1+Cobertura_dossel+Densidade_For+DAP
gml<-glm(S~1+Densidade_For+DAP, data = dados, family = "poisson")
summary(gml)
# null model
gml_null<-glm(S ~ 1, data = dados, family = "poisson")
summary(gml_null)
######
#  GLMM
######
# analyze random effects and fixed effects
modelo.misto.riq.cheio <- glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados, family = poisson)
summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados, family = "poisson")
glmm2<-glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados, family = "poisson")
glmm3<-glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados, family = "poisson")
glmm4<-glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados, family = "poisson")
glmm5<-glmer(S ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados, family = "poisson")
glmm6<-glmer(S ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados, family = "poisson")
glmm7<-glmer(S ~  Cobertura_dossel + (1 | Parcela), data = dados, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
# the best model is the one with the lowest AIC, so it's m5
summary(glmm4)
summ(glmm4,  confint = TRUE, digits = 3)
# comparing the two best models: with random effect and without #
anova(glmm4,gml, test="Chisq")
summary(glmm4)
summ(glmm5, scale=T) # put the response variable in units of standard deviation
summary(gml) #best model
summ(gml)
####### graphic for what was significant ######
###########################################
# Tree density and species richness # 
ggplot(dados, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(dados, aes(x = Densidade_For, y = S, color = Parcela)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = FALSE) +
  geom_point(aes(x = Densidade_For, y = S, fill = Parcela), shape = 21, size = 1.5, inherit.aes = F) +
  geom_smooth(aes(x = Densidade_For, y = S), method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, inherit.aes = F, size = 1.25, color = "black", linetype = "longdash") +
  labs(x = "Number of trees", y = "Richness") +
  theme_apa() 

ggplot(dados, aes(x = DAP, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at Breast Height (m)", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(dados, aes(x = DAP, y = S, color = Parcela)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = FALSE) +
  geom_point(aes(x = Densidade_For, y = S, fill = Parcela), shape = 21, size = 1.5, inherit.aes = F) +
  geom_smooth(aes(x = Densidade_For, y = S), method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, inherit.aes = F, size = 1.25, color = "black", linetype = "longdash") +
  labs(x = "DAP", y = "Richness") +
  theme_apa() 
################################
################################################
#############################################
# 
# abundance analysis # 
#
#############################################
#############################################
# GLM  and GLMM

head(dados)
Modelos.abd <- glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Abd~1+Densidade_For+DAP
best.model<-glm(Abd~1+Densidade_For+DAP, data=dados, family = "poisson") 
summary(best.model)
summ(best.model)
# null model
glm0_abd<-glm(Abd~1, data = dados, family = poisson)
summary(glm0_abd)
anova(best.model, glm0_abd, refit = F, test = "Chisq")
######################################
#  GLMM
######################################
# analyze random effects and mixed effects
######################################
modelo.misto.abd <- glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados, family = poisson)
summary(modelo.misto.abd)
glmm0<-glmer(Abd~1 + (1 | Parcela), data = dados, family = "poisson")
glmm1 <- glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados, family = "poisson")
glmm2<-glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados, family = "poisson")
glmm3<-glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados, family = "poisson")
glmm4<-glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados, family = "poisson")
glmm5<-glmer(Abd ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados, family = "poisson")
glmm6<-glmer(Abd ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados, family = "poisson")
glmm7<-glmer(Abd ~  Cobertura_dossel + (1 | Parcela), data = dados, family = "poisson")
AICctab(glmm0, glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
# # best model is the one with the lowest AIC
summary(glmm4)
summ(glmm4)
# comparing the two best models: with random effect and without
anova(glmm4,best.model, test="Chisq")
summary(glmm4) # this is best model
summ(glmm4)
AICctab(glmm4,best.model)
######################
##### graphics 
######################

ggplot(dados, aes(x = Densidade_For, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(dados, aes(x = DAP, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at Breast Height (m)", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

