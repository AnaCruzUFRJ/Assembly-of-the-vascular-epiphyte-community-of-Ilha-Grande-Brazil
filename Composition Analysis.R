###############################################
#                                             #
# Date: September-2021
# Name: Principal Coordinate Analysis, test Multivariate homogeneity of groups dispersions (variances), and other species composition analysis tests
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description: analyze the variation in the composition of epiphyte species in the RAPELD plots of Ilha Grande, Brazil
# 
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
#############################################################

dados<-read.table("MATRIZGERAL2.csv",header=T,sep=";")

#View(dados)

dados2<-dados[,-c(1:2)]

#View(dados2)
Plots<-dados[,2]

rowSums(dados2)
dados3<-decostand(dados2,"hellinger") #transform to solve asymmetry problems in abundance
#View(dados3)
dist.bray <- vegdist(dados3, method = "bray", binary = F)
dist.bray
#View(as.matrix(dist.bray))

cmdscale(dist.bray)

cmdscale(dist.bray,k=2)

pcoa<-cmdscale(dist.bray,k=2,eig=F) # to see the eigenvalues of each axis

pcoa
plot(pcoa)

### putting the names of the parcels in the pcoa ### 

plot(pcoa,main= "", type="n",xlab="PCoA1",ylab="PCoA2", ylim = c(-0.5, 0.3))
points(pcoa[1:26,1],pcoa[1:26,2],pch=21,bg="#f4f9f3", cex = 1)# plot 1
points(pcoa[27:39,1],pcoa[27:39,2],pch=21,bg="#9ec39a", cex = 1)# plot 2
points(pcoa[40:49,1],pcoa[40:49,2],pch=21,bg="#8af97f", cex = 1)# plot 3
points(pcoa[50:62,1],pcoa[50:62,2],pch=21,bg="#347c2c", cex = 1)# plot 4
points(pcoa[63:75,1],pcoa[63:75,2],pch=21,bg="#072c04", cex = 1) # plot 5


########################################################## testing whether the significant difference by the test Multivariate homogeneity of groups dispersions (variances) ##########################################################
#View(dados)
parc<-as.factor(dados$Plots)
parc
comparacao<- betadisper(d=dist.bray, group=parc)

results<-permutest(comparacao,pairwise =T, digits=2)
#View(results)
results
boxplot(comparacao)
plot(comparacao)
###################################################
# another way to make the PCOA graphic # 
###################################################
head(dados)
dados2 <- dados %>% select(-Segmento, -Plots)
#View(dados2)
Plots <- dados[,2]
Plots
dados3 <- decostand(dados2,"hellinger") 
dados3
dist.bray <- vegdist(dados3, method = "bray", binary = F)
dist.bray
#View(as.matrix(dist.bray))

pcoa <- cmdscale(dist.bray,k = 2,eig = F)
pcoa_auto <- cmdscale(dist.bray, k = 2, eig = T) # to see the eigenvalues of each axis
round((pcoa_auto$eig[1]/sum(pcoa_auto$eig))*100, 1) # eigenvalue of pcoa axis 1
round((pcoa_auto$eig[2]/sum(pcoa_auto$eig))*100, 1)  # # eigenvalue of pcoa axis 2

grafico <- data.frame( "Plots" = as.factor(dados$Plots), "PCoA1" = pcoa[, 1], "PCoA2" = pcoa[, 2])
#View(grafico)

colores <- c("#f4f9f3", "#9ec39a", "#8af97f", "#347c2c", "#072c04")
cores <- c("#f4f9f3", "#9ec39a", "#8af97f", "#347c2c", "#072c04")[grafico$Plots]
# it's in the order of the ggplot legend # 

grafico$cores <- cores

par(mar = c(5, 4, 4, 2) + 0.1)

ggplot(grafico, aes(x = PCoA1, y = PCoA2, fill = Plots)) +
  geom_point(shape = 21, size = 3) +
  #scale_fill_brewer(palette = "Greens") +
  scale_fill_manual(values = colores) + 
  labs(x = "PCoA1 (19.9%)", y = "PCoA2 (17.6%)") +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", size = .5, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", size = .5, color = "grey") +
  theme(panel.grid = element_blank())
#ggsave("grafico_pcoa_matriz_geral.jpeg", w = 6, h = 5, dpi = 300)

###################################################
############ similarity with bray curtis distance ###########
###################################################

clust<-vegdist (dados2, method = "bray")
clust 
 # dissimilarity values #
par(mfrow=c(1,1)) 
# plot(hclust(clust, method="average"), hang=-1, main = "Bray Curtis's Dendrogram", xlab = "Plots")
################################################################
################################################################################################################################### canonical analysis of principal coordinates#####################################################################################################################################################################################
#############
# Only CAP gives the flexibility to use any measure of distance, similarity or dissimilarity desired , such as the Bray-Curtis dissimilarity.# 
##############
setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/Manuscrito 2/Repositorio_Scripts")

# matriz de especies
especiescompleto<-read.table("MATRIZGERAL2.csv",header=T,sep=";")

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


###############################################################

####################################################################################################################################################################################################
#                                    #
#                                    #
# analyzing similarity using bray curtis distance #
#                                    #
#                                    #
######################################################################################################################################################################################################
dados<-read.table("MATRIZGERAL_braycurtis.csv",header=T,sep=";")
# View(dados)
dados<- dados %>% select(-Segmento, -Plots)

str(dados)
clust<-vegdist (dados, method = "bray")
clust # dissimilarity values
par(mfrow=c(1,1)) 
plot(hclust(clust, method="average"), hang=-2, main = "", xlab = "Zones") # we chose not to use the dendrogram in the manuscript
