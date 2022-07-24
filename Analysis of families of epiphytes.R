###############################################
#                                             #
# Date: September-2021
# Name: Effect of environmental conditions and tree morphometry on the diversity of each epiphyte family 
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description: analysis of the effect of environmental conditions and tree morphometry on the diversity of each epiphyte family through generalized linear models and mixed GLM
# 
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
library(car) # normality test
library(MASS)# boxcox
library(ggthemes)# graphics 
library(AID)# boxcox
library(broom) # graphics 
library(ggplot2) # graphics
library(RColorBrewer)#graphics 
library(RVAideMemoire) # to plot residual normality graphs
library(lme4)# linear models
library(bbmle) # we use the AICctab function from the bbmle package
library(rJava)
library(glmulti) # select best models
library(ggpubr)#join figures
library(hrbrthemes)#graphics
library(metafor)
library(jtools)# model output
library(multcomp) # for the glm tukey post hoc test

dados.not.scaled<-read.table("Dados_conjunto3.csv",header=T,sep=";", dec=",")
head(dados.not.scaled)
dados <- as.data.frame(scale(dados.not.scaled)) #using standard deviation units because it has information in different measurement units

#View(dados)
head(dados)
str(dados)

# Bromeliaceae: line 50
# Polypodiaceae: line 170
# Orchidaceae: line 400 
# Araceae: line 570
# Cactaceae: line 741
# Other epiphytes: line 886
# Other ferns: line 1027

###########################################
####### first family: Bromeliaceae  ####### 
###########################################

bromelias<-read.table("Riqueza_Bromeliaceae_por_segmento.csv",header=T,sep=";", dec=",")
str(bromelias)
#View(bromelias)

# testing normality #
# Richness
shapiro.test(bromelias$Riqueza)
 # not normal
mS <- lm(Riqueza ~ 1, data = bromelias)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
  # do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(bromelias$Abundancia)
# not normal
mAbd <- lm(Abundancia ~ 1, data = bromelias)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal
# Richness analysis #
dados.bromelias<-data.frame(bromelias, dados)
#View(dados.bromelias)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.bromelias, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Riqueza~1+Densidade_For+HF 
gml.riq.brom<-glm(Riqueza ~ Densidade_For + HF, data = dados.bromelias, family = "poisson")
summ(gml.riq.brom)
summary(gml.riq.brom) # 
# null model 
glm0_S<-glm(Riqueza~1, data = dados.bromelias, family = poisson)
summary(glm0_S)
AICctab(glm0_S,gml.riq.brom)
anova(gml.riq.brom, glm0_S, refit = T, test = "Chisq")

############ GLMM ##########

modelo.misto.riq.cheio <- glmer(Riqueza ~ + Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.bromelias, family = poisson)

summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.bromelias, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm7,gml.riq.brom, test="Chisq")
summ(glmm7)
summary(glmm7)
summary(gml.riq.brom) 
summ(gml.riq.brom) #this is the best model
# graphic 
ggplot(dados.bromelias, aes(x = HF, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Richness of bromeliads") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
#####################################################
# bromeliad abundance analysis now #
#####################################################
head(dados.bromelias)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.bromelias, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+Densidade_For
gml.abd.brom<-glm(Abundancia ~ Densidade_For, data = dados.bromelias, family = "poisson")
summ(gml.abd.brom)
summary(gml.abd.brom)
# null model
glm0_Abd<-glm(Abundancia~1, data = dados.bromelias, family = poisson)
summary(glm0_Abd)
AICctab(glm0_Abd,gml.abd.brom)
anova(gml.abd.brom, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.bromelias, family = poisson)
summary(modelo.misto.abd.cheio)
summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.bromelias, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.bromelias, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm5,gml.abd.brom, test="Chisq")
summ(glmm5)
summary(glmm5) # this is the best model
# with plot influence (mixed effect)
###########################################
# graphic 
###########################################
ggplot(dados.bromelias, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance of bromeliads") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
##########################################
# checking the effect on each plot #
ggplot(dados.bromelias, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Densidade_For", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())
################################################################################################################################################################################################################################################################################################################################################################################################
# POLYPODIACEAE #
################################################################################################################################################################################################################################################################################################################################################################################################

dados.not.scaled<-read.table("Dados_conjunto3.csv",header=T,sep=";", dec=",")
head(dados.not.scaled)
dados <- as.data.frame(scale(dados.not.scaled)) #using standard deviation units because it has information in different measurement units
head(dados)

samambaias<-read.table("Riqueza_Polypodiaceae_por_segmento.CSV",header=T,sep=";", dec=",")
str(samambaias)
#View(samambaias)
# testing normality 
# Richness
shapiro.test(samambaias$Riqueza)
# not normal 
mS <- lm(Riqueza ~ 1, data = samambaias)
summary(mS)
shapiro.test(resid(mS))
#plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(samambaias$Abundancia)
# not normal 
mAbd <- lm(Abundancia ~ 1, data = samambaias)
summary(mAbd)
shapiro.test(resid(mAbd))
# plotresid(mAbd, shapiro = T)
# not normal 
# Richness analisis#
dados.samambaias<-data.frame(samambaias, dados)
#View(dados.samambaias)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.samambaias, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
#  Riqueza~1+Densidade_For+PROF_COPA
gml.riq.samambaias<-glm(Riqueza~1+Densidade_For+PROF_COPA, data = dados.samambaias, family = "poisson")
summ(gml.riq.samambaias)
summary(gml.riq.samambaias)
# null model 
glm0_S<-glm(Riqueza~1, data = dados.samambaias, family = "poisson")
summary(glm0_S)
AICctab(glm0_S,gml.riq.samambaias)
summ(glm0_S)
anova(gml.riq.brom, gml.riq.samambaias, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.samambaias, family = poisson)

#summary(modelo.misto.riq.cheio)
#summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.samambaias, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
summary(glmm5)
summ(glmm5)
anova(glmm5,gml.riq.samambaias, test="Chisq")
summary(gml.riq.samambaias) # this is the best model
# no part effect #
#############################################
ggplot(dados.samambaias, aes(x = Densidade_For, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
#
ggplot(dados.samambaias, aes(x = PROF_COPA, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown depth", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
#############################################
# analise of abundance in the Polypodiceae now #
#############################################
head(dados.samambaias)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.samambaias, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Abundancia~1+Inclinacao+Densidade_For+PROF_COPA
gml.abd.samambaias<-glm(Abundancia~1+Inclinacao+Densidade_For+PROF_COPA, data = dados.samambaias, family = "poisson")
summ(gml.abd.samambaias)
summary(gml.abd.samambaias)
# null model
glm0_Abd<-glm(Abundancia~1, data = dados.samambaias, family = "poisson")
summary(glm0_Abd)
AICctab(glm0_S,gml.abd.samambaias)
summ(glm0_Abd)
anova(gml.riq.brom, gml.abd.samambaias, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.samambaias, family = poisson)
#summary(modelo.misto.abd.cheio)
# summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.samambaias, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.samambaias, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm2 ,gml.abd.samambaias, test="Chisq")
summary(glmm2) # this is best model
summ(glmm2)
###########################################
# graphics # 
###########################################
ggplot(dados.samambaias, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.samambaias, aes(x = Inclinacao, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Slope", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.samambaias, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.samambaias, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.samambaias, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

#####################################################
# checking the plots separately
#####################################################
ggplot(dados.samambaias, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.samambaias, aes(x = Inclinacao, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Slope", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.samambaias, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.samambaias, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.samambaias, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

##########################################################################################################################################################################################################################################################################################################################################################################################
# ORCHIDACEAE #
##########################################################################################################################################################################################################################################################################################################################################################################################
orquideas<-read.table("Riqueza_Orchidaceae_por_segmento.csv",header=T,sep=";", dec=",")
str(orquideas)
#View(orquideas)
# Richness
# testing normality #
shapiro.test(orquideas$Riqueza)
#NOT NORMAL
mS <- lm(Riqueza ~ 1, data = orquideas)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(orquideas$Abundancia)
# not normal
mAbd <- lm(Abundancia ~ 1, data = orquideas)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal 
# Richness analisis #
dados.orquideas<-data.frame(orquideas, dados)
#View(dados.orquideas)
# GLM # 

Modelos.riqueza <- glmulti(Riqueza ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.orquideas, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
#  Riqueza~1+Inclinacao+Densidade_For+DAP+PROF_COPA
gml.riq.orq<-glm( Riqueza~1+Inclinacao+Densidade_For+DAP+PROF_COPA, data = dados.orquideas, family = "poisson")
summ(gml.riq.orq)
summary(gml.riq.orq)
# null model
glm0_S<-glm(Riqueza~1, data = dados.orquideas, family = "poisson")
summary(glm0_S)
AICctab(glm0_S,gml.riq.orq)
summ(glm0_S)
anova(gml.riq.orq, glm0_S, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.orquideas, family = poisson)
#summary(modelo.misto.riq.cheio)
# summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.orquideas, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
summary(glmm5)
anova(glmm5,gml.riq.orq, refit = T, test = "Chisq")
summary(gml.riq.orq) # this is best model 
# not effect plot# 
summ(gml.riq.orq)
####################################################
# GRAPHICS#
####################################################
ggplot(dados.orquideas, aes(x = Inclinacao, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Slope", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.orquideas, aes(x = Densidade_For, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.orquideas, aes(x = DAP, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.orquideas, aes(x = PROF_COPA, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "crown depth", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
####################################################
# analyze the abundance of orchids now #
####################################################
head(dados.orquideas)
# GLM # 
Modelos.riqueza <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.orquideas, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+Densidade_For+HF
gml.abd.orq<-glm(Abundancia~1+Densidade_For+HF, data = dados.orquideas, family = "poisson")
summ(gml.abd.orq)
summary(gml.abd.orq)
# null model # 
glm0_Abd<-glm(Abundancia~1, data = dados.orquideas, family = "poisson")
summary(glm0_Abd)
AICctab(glm0_Abd,gml.riq.orq)
summ(glm0_Abd)
anova(gml.abd.orq, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.orquideas, family = poisson)
#summary(modelo.misto.abd.cheio)
# summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.orquideas, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.orquideas, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm3,gml.abd.orq, test="Chisq")
summary(glmm3) # this is best model #
# effect plots # 
summ(glmm3)
###########################################
# graphics #
###########################################
ggplot(dados.orquideas, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.orquideas, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.orquideas, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
###### isolated plots######
ggplot(dados.orquideas, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.orquideas, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.orquideas, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundancia") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

##########################################################################################################################################################################################################################################################################################################################################################################################
# Araceae #
##########################################################################################################################################################################################################################################################################################################################################################################################
araceae<-read.table("Riqueza_Araceae_por_segmento.csv",header=T,sep=";", dec=",")
str(araceae)
#View(araceae)
# testing normality #
# Richness
shapiro.test(araceae$Riqueza)
# not normal 
mS <- lm(Riqueza ~ 1, data = araceae)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(araceae$Abundancia)
#nao normais 
mAbd <- lm(Abundancia ~ 1, data = araceae)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal
# Richness analysis #
dados.araceae<-data.frame(araceae, dados)
#View(dados.araceae)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.araceae, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Riqueza~1+Cobertura_dossel+DAP+PROF_COPA
gml.riq.araceae<-glm(Riqueza~1+Cobertura_dossel+DAP+PROF_COPA, data = dados.araceae, family = "poisson")
summ(gml.riq.araceae)
summary(gml.riq.araceae) # 
# null model
glm0_S<-glm(Riqueza~1, data = dados.araceae, family = "poisson")
summary(glm0_S)
AICctab(glm0_S,gml.riq.araceae)
summ(glm0_S)
anova(gml.riq.araceae, glm0_S, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~ + Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.araceae, family = poisson)
summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.araceae, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm7,gml.riq.araceae, test="Chisq")
summary(glmm7) # este eh o melhor modelo 
summ(glmm7)
###################################################
# graphic
###################################################
ggplot(dados.araceae, aes(x = Cobertura_dossel, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

# plots isolated #
ggplot(dados.araceae, aes(x = Cobertura_dossel, y = Riqueza)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())
#####################################################
# analisis abundance #
#####################################################
head(dados.araceae)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.araceae, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+Cobertura_dossel+Densidade_For+PROF_COPA
gml.abd.araceae<-glm(Abundancia~1+Cobertura_dossel+Densidade_For+PROF_COPA, data = dados.araceae, family = "poisson")
summ(gml.abd.araceae)
summary(gml.abd.araceae)
# null model 
glm0_Abd<-glm(Abundancia~1, data = dados.araceae, family = "poisson")
summary(glm0_Abd)
AICctab(glm0_Abd,gml.abd.araceae)
summ(glm0_Abd)
anova(gml.abd.araceae, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.araceae, family = poisson)
summary(modelo.misto.abd.cheio)
summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.araceae, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.araceae, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm1,gml.abd.araceae, test="Chisq")
summary(glmm1) # this is best model 
summ(glmm1)
# with effect plot 
###################################################
# graphic
###################################################

ggplot(dados.araceae, aes(x = Cobertura_dossel, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.araceae, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.araceae, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.araceae, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.araceae, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.araceae, aes(x = DIAM_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown diameter", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

# plots isolated #

ggplot(dados.araceae, aes(x = Cobertura_dossel, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.araceae, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.araceae, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.araceae, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.araceae, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.araceae, aes(x = DIAM_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Crown diameter", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())
##########################################################################################################################################################################################################################################################################################################################################################################################
# Cactaceae #
##########################################################################################################################################################################################################################################################################################################################################################################################
cactus<-read.table("Riqueza_cactaceae_por_segmento.csv",header=T,sep=";", dec=",")
str(cactus)
# View(cactus)
# testing normality #
# Richness
shapiro.test(cactus$Riqueza)
# not normal
mS <- lm(Riqueza ~ 1, data = cactus)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(cactus$Abundancia)
# not normal
mAbd <- lm(Abundancia ~ 1, data = cactus)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal
# Richness analisis #
dados.cactus<-data.frame(cactus, dados)
View(dados.cactus)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.cactus, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Riqueza~1+DAP
gml.riq.cactus<-glm(Riqueza~1+DAP, data = dados.cactus, family = "poisson")
summ(gml.riq.cactus)
summary(gml.riq.cactus) 
# null model #
glm0_S<-glm(Riqueza~1, data = dados.cactus, family = "poisson")
AICctab(glm0_S,gml.riq.cactus)
summ(glm0_S)
anova(gml.riq.cactus, glm0_S, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~ + Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.cactus, family = poisson)
summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.cactus, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
summ(glmm7)
anova(glmm7,gml.riq.cactus, test="Chisq")
anova(glmm7, gml.riq.cactus, refit = T, test = "Chisq")
summ(gml.riq.cactus) # this is best model
summary(gml.riq.cactus)
###########
# graphic #
###########
ggplot(dados.cactus, aes(x = DAP, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
#########################
# Abundance analisis
#########################
head(dados.cactus)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.cactus, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+DAP+HF
gml.abd.cactus<-glm(Abundancia~1+DAP+HF, data = dados.cactus, family = "poisson")
summ(gml.abd.cactus)
summary(gml.abd.cactus)
# null model
glm0_Abd<-glm(Abundancia~1, data = dados.cactus, family = "poisson")
AICctab(glm0_Abd,gml.abd.cactus)
summ(glm0_Abd)
anova(gml.abd.cactus, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.cactus, family = poisson)
summary(modelo.misto.abd.cheio)
summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.cactus, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.cactus, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm1,gml.abd.cactus, test="Chisq")
summary(glmm1) # This is best model 
summ(glmm1)
# with effect plot
###########
# graphic #
###########
ggplot(dados.cactus, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.cactus, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.cactus, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
##############
# plots isolated
##############
ggplot(dados.cactus, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.cactus, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(dados.cactus, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = T, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_bw() +
  theme(panel.grid = element_blank())
##########################################################################################################################################################################################################################################################################################################################################################################################
# Other epiphytes #
##########################################################################################################################################################################################################################################################################################################################################################################################
others<-read.table("Riqueza_outras familias_por_segmento.csv",header=T,sep=";", dec=",")
str(others)
#View(others)
# testing normality #
# Richness
shapiro.test(others$Riqueza)
# not normal
mS <- lm(Riqueza ~ 1, data = others)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(others$Abundancia)
# not normal
mAbd <- lm(Abundancia ~ 1, data = others)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal
# Richness analisis #
dados.others<-data.frame(others, dados)
head(dados.others)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.others, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Riqueza~1+DAP
gml.riq.others<-glm(Riqueza~1+DAP, data = dados.others, family = "poisson")
summ(gml.riq.others)
summary(gml.riq.others) # 
# null model
glm0_S<-glm(Riqueza~1, data =dados.others, family = "poisson")
AICctab(glm0_S,gml.riq.others)
summ(glm0_S)
anova(gml.riq.others, glm0_S, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~ + Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others, family = poisson)
summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.others, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.others, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.others, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.others, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.others, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.others, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm7,gml.riq.others, test="Chisq")
summary(gml.riq.others) # this is best model
#########
# graphic 
#########
ggplot(dados.others, aes(x = DAP, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
##################################################
# abundance analisis #
##################################################
head(dados.others)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.others, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+DAP+HF
gml.abd.others<-glm(Abundancia~1+Densidade_For+DAP+HF+PROF_COPA+DIAM_COPA, data = dados.others, family = "poisson")
summ(gml.abd.others)
summary(gml.abd.others)
############################
# null model #
glm0_Abd<-glm(Abundancia~1, data = dados.others, family = "poisson")
AICctab(glm0_Abd,gml.abd.others)
summ(glm0_Abd)
anova(gml.abd.others, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others, family = poisson)
summary(modelo.misto.abd.cheio)
summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.others, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.others, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.others, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.others, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.others, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.others, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm1,gml.abd.others, test="Chisq")
summary(glmm1) # this is best model 
summ(glmm1)
# with effect plot #

###################
# graphics #
###################
ggplot(dados.others, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = Inclinacao, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Slope", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = Cobertura_dossel, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = HF, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Trunk height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = PROF_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown depth", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = DAP, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Diameter at breast height", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others, aes(x = DIAM_COPA, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown diameter", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()
##########################################################################################################################################################################################################################################################################################################################################################################################
# Other ferns #
##########################################################################################################################################################################################################################################################################################################################################################################################
others.ferns<-read.table("Riqueza_outras samambaias_por_segmento.csv",header=T,sep=";", dec=",")
str(others.ferns)
#View(others.ferns)
# testing data normality #
# Richness #
shapiro.test(others.ferns$Riqueza)
#not normal
mS <- lm(Riqueza ~ 1, data = others.ferns)
summary(mS)
shapiro.test(resid(mS))
plotresid(mS, shapiro = T)
# BOX_COX<-dados %>%
# do(COX=boxcoxnc(.$S))
# abundance
shapiro.test(others.ferns$Abundancia)
# not normal
mAbd <- lm(Abundancia ~ 1, data = others.ferns)
summary(mAbd)
shapiro.test(resid(mAbd))
plotresid(mAbd, shapiro = T)
# not normal
# Richness analisis
dados.others.ferns<-data.frame(others.ferns, dados)
head(dados.others.ferns)
# GLM # 
Modelos.riqueza <- glmulti(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.others.ferns, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Riqueza~1+Cobertura_dossel+Inclinacao+Densidade_For+DAP+HF+PROF_COPA
gml.riq.others.ferns<-glm(Riqueza~1+Cobertura_dossel+Inclinacao+Densidade_For+DAP+HF+PROF_COPA, data = dados.others.ferns, family = "poisson")
summ(gml.riq.others.ferns)
summary(gml.riq.others.ferns) # 
# null model 
glm0_S<-glm(Riqueza~1, data = dados.others.ferns, family = "poisson")
AICctab(glm0_S,gml.riq.others.ferns)
summ(glm0_S)
anova(gml.riq.others.ferns, glm0_S, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.riq.cheio <- glmer(Riqueza ~ + Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others.ferns, family = poisson)
summary(modelo.misto.riq.cheio)
summ(modelo.misto.riq.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm2<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm3<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm4<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm5<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm6<-glmer(Riqueza ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm7<-glmer(Riqueza ~  Cobertura_dossel + (1 | Parcela), data = dados.others.ferns, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
summary(glmm5)
anova(glmm5,gml.riq.others.ferns, test="Chisq")
summary(gml.riq.others.ferns) # this is best model 
summ(gml.riq.others.ferns)
##########
# graphics # 
##########
ggplot(dados.others.ferns, aes(x = Cobertura_dossel, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others.ferns, aes(x = Inclinacao, y = Riqueza)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Slope", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others.ferns, aes(x = Densidade_For, y = Riqueza)) + geom_point(shape = 19, size = 2) + geom_smooth( method = glm, se=T) +
  labs( x = "Densidade_For", y = "Riqueza de outras samambaias") + theme_bw() + 
  scale_color_brewer(palette="Set1")
#####################################################
# abundance analisis #
#####################################################
head(dados.others.ferns)
# GLM # 
Modelos.abd <- glmulti(Abundancia  ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados.others.ferns, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Abundancia~1+Cobertura_dossel+Inclinacao+Densidade_For+DAP+HF+PROF_COPA
gml.abd.others.ferns<-glm(Abundancia~1+Cobertura_dossel+Inclinacao+Densidade_For+DAP+HF+PROF_COPA, data = dados.others.ferns, family = "poisson")
summ(gml.abd.others.ferns)
summary(gml.abd.others.ferns)
# null model
glm0_Abd<-glm(Abundancia~1, data = dados.others.ferns, family = "poisson")
AICctab(glm0_Abd,gml.abd.others.ferns)
summ(glm0_Abd)
anova(gml.abd.others.ferns, glm0_Abd, refit = T, test = "Chisq")
############ GLMM ##########
modelo.misto.abd.cheio <- glmer(Abundancia ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others.ferns, family = poisson)
summary(modelo.misto.abd.cheio)
summ(modelo.misto.abd.cheio,  confint = TRUE, digits = 3)
glmm1 <- glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm2<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm3<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm4<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm5<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + Densidade_For  + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm6<-glmer(Abundancia ~  Cobertura_dossel + Inclinacao  + (1 | Parcela), data = dados.others.ferns, family = "poisson")
glmm7<-glmer(Abundancia ~  Cobertura_dossel + (1 | Parcela), data = dados.others.ferns, family = "poisson")
AICctab(glmm1,glmm2, glmm3, glmm4, glmm5, glmm6, glmm7, base = T, weights = T)
anova(glmm5,gml.abd.others.ferns, test="Chisq")
summary(glmm5)
summary(gml.abd.others.ferns) # this is best model 
summ(gml.abd.others.ferns)
# not effect plot 
###########################################
##########
# graphics # 
##########
ggplot(dados.others.ferns, aes(x = Cobertura_dossel, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others.ferns, aes(x = Inclinacao, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Slope", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

ggplot(dados.others.ferns, aes(x = Densidade_For, y = Abundancia)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

