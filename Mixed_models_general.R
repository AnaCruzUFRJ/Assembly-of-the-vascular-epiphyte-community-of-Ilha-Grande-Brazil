
###############################################
#                                             #
# Date: 19-July-2022
# Name: Effect of environmental variables and tree morphometric in the epiphyte of the Brazilian Atlantic Forest 
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description: analysis of the effect of environmental and morphometric variations of trees on the richness and abundance of vascular epiphytes in 5 RAPELD plots of Ilha Grande, Brazil
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
dados1<-read.table("Dados_conjunto1.csv",header=T,sep=";", dec=",")
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

###################
# descriptive analysis
# Richness
Rest1<-dados %>% filter (Parcela == "1")
mean<-Rest1[,13]
describe(Rest1[,13])
Rest2<-dados %>% filter (Parcela == "2")
describe(Rest2[,13])
Baixas<-dados %>% filter (Parcela == "3")
describe(Baixas[,13])
Sub1<-dados %>% filter (Parcela == "4")
describe(Sub1[,13])
Sub2<-dados %>% filter (Parcela == "5")
describe(Sub2[,13])
# Let's compare if there is a difference between the wealth of the first plots #
ggplot(dados, aes(x=Parcela, y=S, fill=Parcela)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Plots",y="Richness") + 
  theme_bw() +
  theme(panel.grid = element_blank())

###################

# Mean and standard deviation graph 

Richness <- c(3.81,6.54,4.8,4.77,9.92)
sd <- c(2.5,2.03,1.81,1.88,3.95)
Plots<-c(1,2,3,4,5)

dat_plots<- data.frame(Plots, Richness, sd)
dat_plots$Plots<-as.factor(dat_plots$Plots)
str(dat_plots$Plots)

p_media_desvio <- ggplot(dat_plots, aes(factor(Plots), Richness)) + 
  geom_point(stat = "identity", size = 2) + 
  geom_errorbar(aes(ymin = Richness - sd, ymax =  Richness + sd), width = 0.2) +
  xlab("Forests")+ 
  theme_apa() 

p_media_desvio

###################
# difference na abundance now
###################
# descriptive analysis
Rest1<-dados %>% filter (Parcela == "1")
mean<-Rest1[,14]
describe(Rest1[,14])
Rest2<-dados %>% filter (Parcela == "2")
describe(Rest2[,14])
Baixas<-dados %>% filter (Parcela == "3")
describe(Baixas[,14])
Sub1<-dados %>% filter (Parcela == "4")
describe(Sub1[,14])
Sub2<-dados %>% filter (Parcela == "5")
describe(Sub2[,14])
###################
ggplot(dados, aes(x=Parcela, y=Abd, fill=Parcela)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Plots",y="Abundance") + 
  theme_bw() 
###################
# Mean and standard deviation graph 
Abundance <- c(19.62,34.77,19.8,14.08,36.69)
sd <- c(16.94,17.20,17.37,10.48,17.4)
Plots<-c(1,2,3,4,5)
dat_plots<- data.frame(Plots, Abundance, sd)
dat_plots$Plots<-as.factor(dat_plots$Plots)
str(dat_plots$Plots)
p_media_desvio <- ggplot(dat_plots, aes(factor(Plots), Abundance)) + 
  geom_point(stat = "identity", size = 2) + 
  geom_errorbar(aes(ymin = Abundance - sd, ymax =  Abundance + sd), width = 0.2) +
  xlab("Forests")+ 
  theme_apa() 
p_media_desvio
###################
glm.abd.parcela<-glm(Abd~Parcela, data=dados)
summ(glm.abd.parcela)
summary(glht (glm.abd.parcela, mcp (Parcela = "Tukey"))) #perform TukeyHSD

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
# Best model: S~1+Cobertura_dossel+Densidade_For+DAP
gml3<-glm(S ~ Cobertura_dossel  + Densidade_For  + DAP, data = dados, family = "poisson")
summary(gml3)
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
summary(glmm5)
summ(glmm5,  confint = TRUE, digits = 3)
# comparing the two best models: with random effect and without #
anova(glmm5,gml3, test="Chisq")
summary(glmm5)
summ(glmm5, scale=T) # put the response variable in units of standard deviation
summ(glmm5, confint = TRUE, digits = 2) # to report confidence intervals
# there is a difference between them, so we are left with the most complex, with random effect.
####### graphic for what was significant ######
###########################################
# Tree density and species richness # 
ggplot(dados, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa() + theme(axis.title = element_text(size = 8),
                        axis.text = element_text(size = 8))


ggplot(dados, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  facet_wrap(~ Parcela) +
  theme_base()

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
################################
#### GLM #### 
#### for the plots separately #### 
################################
head(dados)
# View(dados)
# restinga 1 #
rest1<- dados%>% filter(Parcela =="1")
# View(rest1)
glm_rest1 <- glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = rest1, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1+Cobertura_dossel+Densidade_For+DIAM_COPA
best_rest1<-glm (S ~ Cobertura_dossel  + Densidade_For  + DIAM_COPA, data = rest1, family = "poisson")
summary(best_rest1)
# null model
glm0_rest1<- glm (S ~ 1, data = rest1, family = "poisson")
summary(glm0_rest1)
AICctab(best_rest1,glm0_rest1)
anova(best_rest1, glm0_rest1, refit = T, test = "Chisq")

# graphic 
ggplot(rest1, aes(x = Cobertura_dossel, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover (%)", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(rest1, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(rest1, aes(x = DIAM_COPA, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Crown diameter (m)", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

# restinga 2 #
rest2<- dados%>% filter(Parcela =="2")
# View(rest2)
glm_rest2 <- glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = rest2, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Best model: S~1+Densidade_For
best_rest2<-glm (S ~ Densidade_For, data = rest2, family = "poisson")
summary(best_rest2)
summ(best_rest2)


# null model
glm0_rest2<- glm (S ~ 1, data = rest2, family = "poisson")
summary(glm0_rest2)
AICctab(best_rest2,glm0_rest2)
anova(best_rest2, glm0_rest2, refit = T, test = "Chisq")
# no difference from the null model
ggplot(rest2, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

# FlorestaTerrasBaixas #
terras.baixas<- dados%>% filter(Parcela =="3")
# View(terras.baixas)
glm_terras.baixas <- glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = terras.baixas, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1+Densidade_For
best_terras_baixas<-glm (S ~ Densidade_For, data = terras.baixas, family = "poisson")
summary(best_terras_baixas)
summ(best_terras_baixas)
# null model
glm0_terras_baixas<- glm (S ~ 1, data = terras.baixas, family = "poisson")
summary(glm0_terras_baixas)
AICctab(best_terras_baixas,glm0_terras_baixas)
anova(best_terras_baixas, glm0_terras_baixas, refit = T, test = "Chisq")
ggplot(terras.baixas, aes(x = Densidade_For, y = S)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Richness") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()
# FlorestaSubmontana1 #
submontana1<- dados%>% filter(Parcela =="4")
# View(submontana1)
glm_submontana1<-glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = submontana1, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1
glm_submontana1<-glm (S ~ 1, data = submontana1, family = "poisson")
summary(glm_submontana1)
 # none of the variables is significant
# null model
glm0_submontana1<- glm (S ~ 1, data = submontana1, family = "poisson")
summary(glm0_terras_baixas)

AICctab(glm_submontana1,glm0_submontana1)
anova(glm_submontana1, glm0_submontana1, refit = T, test = "Chisq")

# FlorestaSubmontana2 #
submontana2<- dados%>% filter(Parcela =="5")
# View(submontana2)
glm_submontana2<-glmulti(S ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = submontana2, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1
glm_submontana2<-glm (S ~ 1, data = submontana2, family = "poisson")
summary(glm_submontana2)
# null model
glm0_submontana2<- glm (S ~ 1, data = submontana2, family = "poisson")
summary(glm0_submontana2)
# none of the variables is significant

#############################################
#############################################
# 
# abundance analysis # 
#
#############################################
#############################################
# GLM  and GLMM

head(dados)
Modelos.abd <- glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = dados, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Abd~1+Cobertura_dossel+Densidade_For+DAP+DIAM_COPA
# model selection
best.model<-glm(Abd ~ Cobertura_dossel+Densidade_For+DAP+DIAM_COPA, data=dados, family = "poisson") 
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
summary(glmm5)
summ(glmm5)
# comparing the two best models: with random effect and without
anova(glmm5,best.model, test="Chisq")
summary(glmm5) # this is best model
AICctab(glmm5,best.model)
######################
##### graphics 
######################

ggplot(dados, aes(x = Densidade_For, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(dados, aes(x = Cobertura_dossel, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover (%)", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

################################################################################################################################
#### GLM #### plots separately #### 
################################################################################################################################
head(dados)
# View(dados)
# restinga 1 #
rest1<- dados%>% filter(Parcela =="1")
# View(rest1)
head(rest1)
Modelos.abd <- glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = rest1, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Abd~1+Cobertura_dossel+Densidade_For
best_rest1<-glm (Abd ~ Cobertura_dossel  + Densidade_For, data = rest1, family = "poisson")
summary(best_rest1)
summ(best_rest1)
# null model
glm0_abd<-glm(Abd~1, data = rest1, family = poisson)
summary(glm0_abd)
AICctab(glm0_abd,best_rest1)
anova(best_rest1, glm0_abd, refit = T, test = "Chisq")
# graphics #
ggplot(rest1, aes(x = Cobertura_dossel, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover (%)", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

ggplot(rest1, aes(x = Densidade_For, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()

# restinga 2 #
rest2<- dados%>% filter(Parcela =="2")
# View(rest2)
glm_rest2 <- glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = rest2, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: Abd~1+Inclinacao+Densidade_For+HF
best_rest2<-glm (S ~Inclinacao+Densidade_For+HF, data = rest2, family = "poisson")
summary(best_rest2)
summ(best_rest2)
# null model
glm0_abd<-glm(Abd~1, data = rest2, family = poisson)
summary(glm0_abd)
AICctab(glm0_abd,best_rest2)
anova(best_rest2, glm0_abd, refit = T, test = "Chisq")
# graphics #
ggplot(rest2, aes(x = Cobertura_dossel, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Canopy cover", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "black") +
  theme_apa()


# FlorestaTerrasBaixas #
terras.baixas<- dados%>% filter(Parcela =="3")
# View(terras.baixas)
glm_terras.baixas <- glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = terras.baixas, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1+Densidade_For
best_terras_baixas<-glm (Abd ~ Densidade_For, data = terras.baixas, family = poisson)
summary(best_terras_baixas)
summ(best_terras_baixas)
# null model 
glm0_abd<-glm(Abd~1, data = terras.baixas, family = poisson)
summary(glm0_abd)
AICctab(glm0_abd,best_terras_baixas)
anova(best_terras_baixas, glm0_abd, refit = T, test = "Chisq")

# graphic # 

ggplot(terras.baixas, aes(x = Densidade_For, y = Abd)) +    geom_point(shape = 21, size = 1.5, fill = "black") +
  labs(x = "Number of trees", y = "Abundance") +
  geom_smooth(method = "glm", 
              method.args = list(family = "poisson"), 
              se = F, size = 1, color = "grey45") +
  theme_apa()

# FlorestaSubmontana1 #
submontana1<- dados%>% filter(Parcela =="4")
# View(submontana1)
glm_submontana1<-glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = submontana1, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1
glm_submontana1<-glm (Abd ~ 1, data = submontana1, family = poisson)
summary(glm_submontana1, family = poisson)
summ(glm_submontana1)
# none of the variables is significant
glm0_abd<-glm(Abd~1, data = submontana1, family = poisson)
summary(glm0_abd)

# FlorestaSubmontana2 #
submontana2<- dados%>% filter(Parcela =="5")
# View(submontana2)
glm_submontana2<-glmulti(Abd ~ Cobertura_dossel + Inclinacao  + Densidade_For  + DAP + HF + PROF_COPA + DIAM_COPA, data = submontana2, level=1, fitfunction = 'glm', crit="aicc", confsetsize=100)
# Best model: S~1
glm_submontana2<-glm (Abd ~ 1, data = submontana2, family = "poisson")
summary(glm_submontana2)
summ(glm_submontana2)
# none of the variables is significant
glm0_abd<-glm(Abd~1, data = submontana2, family = poisson)
summary(glm0_abd)
summ(glm0_abd)
######################################################################################################################################################################################################################################################