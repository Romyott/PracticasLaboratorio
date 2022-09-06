require(ggplot2)
require(GGally) #este comando es igual al library
require(picante)
#if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
#BiocManager::install("phyloseq")
library(phyloseq)
library(ape)
library("phangorn")
#if (!requireNamespace("BiocManager", quietly = TRUE))
 #  install.packages("BiocManager")
#BiocManager::install("dada2")
library("dada2")
library("vegan")
library("factoextra")
library("labdsv")
library("dplyr")
library("corrplot")

setwd("D://Romy//FACU//doctoral//datos//52 lagunas//Correlograma")
#load("D://Romy//FACU//doctoral//datos//52 lagunas//Correlograma//correlograma-RIOSALADO.rdata")#*******

Otus<-read.csv("Core Otus.csv", header = T, row.names = 1)
Amb<-read.csv("LagunasPampeanasArroyos.csv", header=TRUE, row.names=1, sep = ";" , dec=",")#******
Tree <- read.nexus ("arbol-ARROYOGRANDE.nex")#("ArbolTotaltips.nex")#********************************



#Transformo las variables para trabajar#################################################
#necesito que la matriz ambiental sea solo numerica, sino no me sirve el CCA
rownames(Amb) <- Amb[ ,2 ] #asi coinciden los nombres de las lagunas para ambas matrices
AmbNum <- Amb[c(13:45)] #mismo recorte pero para los numericos:
#hago numericos los que eran otra cosa
AmbNum$Pico.Phyto = as.numeric(as.factor(AmbNum$Pico.Phyto))
AmbNum$X..Cropland.per.Department = as.numeric(as.factor(AmbNum$X..Cropland.per.Department))
AmbNum$DOC = as.numeric(as.factor(AmbNum$DOC))
AmbNum$Total.cattle.per.department.2015 = as.numeric(as.factor(AmbNum$Total.cattle.per.department.2015))
AmbNum$Altitude = as.numeric(as.factor(AmbNum$Altitude))
AmbNum$Department.surface..Km2. = as.numeric(as.factor(AmbNum$Department.surface..Km2.))
AmbNum$Altitude = as.numeric(as.factor(AmbNum$Altitude))
AmbNum$Water.Temperature = as.numeric(as.factor(AmbNum$Water.Temperature))
#eliminar las variables redundantes, saco las que saque en el analisis ambiental alla en su momento.
AmbNum2 <- AmbNum[-c(3,5,14,15,16,19,23,24,25)]
AmbEscalado <- scale(AmbNum2, center=TRUE , scale=TRUE)
AmbEscalado <- as.data.frame(AmbEscalado)

#sacar variables ambientales a ver que pasa con los VIF y el test en gral############################

#sacar las variables ambientales segun laguna ************************************************
#AmbNumX <- AmbNum2 [-c(12,19,20)] #RIOSALADO 
AmbNumX <- AmbNum2 #RIOARRECIFES

#sacar lagunas que no necesito para el arroyo*************************************************
#AmbNumX2 <- AmbNumX[-c(4,5,7,8,13,16,17,18,19,20,21,22,23,24,36,37,45), ] #RIO SALADO 
#AmbNumX2 <- AmbNumX[c(8,13), ] #RIOARRECIFES
AmbNumX2 <- AmbNumX[c(4,5,7), ] #ARROYOGRANDE

AmbEscaladoX <- AmbNumX2 

#me arruinan las variables por ser pocas observaciones????
#AmbEscaladoX <- scale(AmbNumX2, center=TRUE , scale=TRUE)
#AmbEscaladoX <- as.data.frame(AmbEscaladoX)



#analisis DCA##########################################################################
#AmbNum3 <-AmbNum2[-c(1,2)] #,13,20)] #Todas las lagunas****************************
#AmbNum3 <- AmbNumX2 [-c(1,2,12,19,20)] #RIOSALADO***********************************

DCA <- decorana(AmbNum3) # ,iweigh=1) # iweigh=1 turns on down-weighting of rare taxa
DCA #la long de 1,16 en el primer eje indica un RDA
summary(DCA)

ordiplot(DCA)#, xlim = c(-2,2), ylim = c(-3,2))
plot1 <- points (DCA, display=c("sites"), choices=1:2, pch=3, col="red", xlim=c(-2,2), ylim=c(-3,2))
plot2 <- text(DCA, display=c("species"), choices=1:2,cex=0.7, xlim=c(-2,2), ylim=c(-3,2))



#Defino las variables para trabajo#####################################################
TaxTree<-read.csv ("TaxaOTUarbol-ARROYOGRANDE.csv", header= T, row.names = 1) #("TaxaArbol.csv", header = T, row.names = 1)
TaxM<-as.matrix(TaxTree)
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(TaxM)
ps0 = phyloseq(OTU, Tree, TAX)
OtusHps <- hellinger(otu_table(ps0))
OtusHpsTras <- data.frame(t(OtusHps[])) 
OtusTras.ps <- na.omit(OtusHpsTras)   

#saco las lagunas que no sean de mi interes*******************************************
#OtusTras <- OtusTras.ps[-c(4,5,7,8,13,16,17,18,19,20,21,22,23,24,36,37,45), ] #RIOSALADO  
#OtusTras <- OtusTras.ps[c(8,13), ] #RIOARRECIFES
OtusTras <- OtusTras.ps[c(4,5,7), ] #ARROYOGRANDE

X= OtusTras
Y= AmbEscaladoX 



#ahora a analizar (CCA)################################################################
CCA<- cca(X, Y) 
CCA
summary(CCA) 
RsquareAdj(CCA)

#grafico los componentes principales del CCA
barplot<-screeplot(CCA, bstick = FALSE, type = c("barplot", "lines"), npcs = min(10, if (is.null(CCA$CCA) || CCA$CCA$rank == 0) CCA$CA$rank else CCA$CCA$rank), ptype = "o", bst.col = "red", bst.lty = "solid", xlab = "Component", ylab = "Inertia", main = deparse(substitute(CCA)), legend = bstick)
barplot2<- screeplot(CCA, bstick = FALSE, type = c("barplot", "lines"), npcs = 4, ptype = "o", bst.col = "red", bst.lty = "solid", xlab = "Component", ylab = "Inertia", main = deparse(substitute(CCA)))
#y despues, los ploteo en nube de puntos
CCA.plot1 <- plot(CCA,choices=c(1,2), display=c("sp","cn"), main = "Biplot CCA Otus") #xlim=c(-2,3), ylim=c(-5,5),
CCA.plot2 <- plot(CCA,choices=c(1,2), display=c("lc","cn"), main = "Biplot CCA Lagunas") #xlim=c(-2,3), ylim=c(-5,5),
CCA.plot3 <- plot(CCA,choices=c(1,2), display=c("sp","lc"), main = "Biplot CCA OtusxLagunas") #xlim=c(-2,3), ylim=c(-5,5),

CCA.plot4 <- plot(CCA,choices=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), display = c("sp","lc"),main="todos los ejes")

CCA.ANOVA <- anova.cca(CCA)
CCA.ANOVA



#Scores spxamb#########################################################################
## Default S3 method:
scores(CCA, display=c("lc","sp"))
scores(CCA.plot3$sp, display=c("sp"))
write.csv(CCA.plot3$species,file.path(getwd(),"OTUscoresHCCA-ps-ARROYOGRANDE.csv")) #*****************************
write.csv(CCA.plot4$species,file.path(getwd(),"CCA-ARROYOGRANDE.csv"))#********************************************



#VIF#######################################################################################

VIFCCA <- vif.cca(CCA)
VIFCCA
summary(VIFCCA)



#Distancia euclidiana ################################################################
OTUxSITE <- CCA.plot3$species

#calculo la distancia y lo pongo en formato matriz
DEucl <- stats::dist(OTUxSITE, method = "euclidean")
DEuclM <-as.matrix(DEucl)

#lleva h por la Hellinger
write.csv(DEuclM,file.path(getwd(),"DistEuclidh-ps-CCA-ARROYOGRANDE.csv"))#************************************



#ps##########################################################################################
ps1 = phyloseq(OTU, Tree, TAX)
tax_table(ps1)
otu_table(ps1)
phy_tree(ps1)
#saco la distancia cofenetica del arbolito
phydist1<-cophenetic(phy_tree(ps1))



#Correlograma########################################################################
#Euclidean<-read.csv("DistEuclidh.csv", header = T, row.names = 1)
#EuclidM <- as.matrix(Euclidean)
DEuclM <- as.matrix(DEucl)
phyM <-as.matrix(phydist1)

DEuclM
phyM

write.csv(DEuclM,file.path(getwd(),"DEuclM-ARROYOGRANDE.csv"))#******************************************
write.csv(phyM,file.path(getwd(),"phyM-ARROYOGRANDE.csv"))#**************************************************
#tengo que ver que sobra de phyM que falta en DEuclM....

Correlogram <-mantel.correlog(DEuclM, phyM, mult="bonferroni", progressive=TRUE)
Correlogram
plot(Correlogram)
write.csv(Correlogram$mantel.res,file.path(getwd(),"correlograma-CCA-ARROYOGRANDE.csv"))#****************************


save.image("D://Romy//FACU//doctoral//datos//52 lagunas//Correlograma//correlograma-ARROYOGRANDE.rdata")#******************




#preparo el ps para RC y bNTI ###########################################################

#eliminar singletons
psSinSigles <- filter_taxa(otu_table(ps1), function (x) {sum(x > 0) > 1}, prune=TRUE)
psSinSigles
#rarefaccionar y no, y ver cual va mejor
#rarefy_even_depth, paquete phyloseq......
set.seed(500); .Random.seed
Otusrarefaccionados <- rarefy_even_depth(otu_table(psSinSigles), rngseed = 500)
Otusrarefaccionados
#Distancia bNTI Stegen#########################################################################################

#intento Stegen 2013, bNTI
## read in the phylogeny
phylo <- phy_tree(ps1)
phylo; # a summary of the phylogeny
#plot.phylo(phylo,typ="fan"); # a quick plot

## make sure the names on the phylogeny are ordered the same as the names in otu table
match.phylo.otu = match.phylo.data(phylo, otu_table(psSinSigles));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'betaMNTD_weighted-psRARE.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
   
   rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
   
   print(c(date(),rep));
   
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
   for (rows in (columns+1):ncol(match.phylo.otu$data)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
      
   };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.table(weighted.bNTI,"weighted_bNTI-psRARE.txt",quote=TRUE, dec=".", sep=",");
pdf("weighted_bNTI_Histogram-psRARE.pdf")
hist(weighted.bNTI)
dev.off()

#selecciono aquellos con bNTI menor a 2

#Raup-Ckrick######################################################################################
#para calcular esto, sacar los de bNTI mayores a modulo de 2, luego correr el RC
#install.packages("iCAMP")
library("iCAMP")

#comm va a ser lo que sobro del bNTI con valores entre el modulo de 2, y ver que factores ambientales tienen asociados
#porque eso es el agente selectivo en la laguna para las bacterias

#comm es laguna por laguna lpm....#$"#%"$&#&%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
comm = beta.mntd.weighted
comm= t(otu_table(psSinSigles)) #si pongo la matriz traspuesta, toma las lagunas. si la pongo normal, toma por taxon
rand.time=500 # usually use 1000 for real data.
nworker=2 # parallel computing thread number
RC=RC.pc(comm=comm, rand = rand.time,
         nworker = nworker, weighted = TRUE,
         sig.index="RC")

RCM <- as.matrix(RC$index)
RCM

hist(RC$index)
write.table(RC$index,file="Raup-Crick-ps.txt")
write.csv(RC,file.path(getwd(),"Raup-Crick-ps.csv"))
write.csv(RCM,file.path(getwd(),"Raup-Crick-psMATRIZ.csv"))

#conteo#################################################################################
#conteos (hasta que arregle las matrices)
install.packages("tidyverse")
library("tidyverse")

RCtri <- RCM[lower.tri(RCM)]
RCtri
RCdata <- as.data.frame(RCtri)
conteo <- RCdata %>% tally(RCdata < -0.95)
hist(RCM)

BETAdata <- as.data.frame(weighted.bNTI)
conteo <- BETAdata %>% tally(BETAdata < 2)
hist(weighted.bNTI)

###################################################################################
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//Correlograma//sesionCorrelograma.rdata")
load("D://Romy//FACU//doctoral//datos//52 lagunas//Correlograma//sesionCorrelograma.rdata")

write.csv(phydist1,file.path(getwd(),"DistanciaFilogenetica.csv"))

#lo que descarto#################################################

#Evaluo la normalidad var amb #############################################################
plot(AmbNum2$Latitude)
shapiro.test(AmbNum2$Latitude)$p.value #muy al limite de 0,05 (da 0,0571)
plot(AmbNum2$Longitude)
shapiro.test(AmbNum2$Longitude)$p.value #no, 0,000345
plot(AmbNum2$Lake.Area.km2)
shapiro.test(AmbNum2$Lake.Area.km2)$p.value #no, 7,48 xe-10
plot(AmbNum2$C.D.)
shapiro.test(AmbNum2$C.D.)$p.value #no, 1,198 xe-05
plot(AmbNum2$Water.Temperature)
shapiro.test(AmbNum2$Water.Temperature)$p.value #normal, 0,1095
plot(AmbNum2$SD)
shapiro.test(AmbNum2$SD)$p.value #no, 1,861 e-11
plot(AmbNum2$DO)
shapiro.test(AmbNum2$DO)$p.value #no, 0,000823
plot(AmbNum2$pH)
shapiro.test(AmbNum2$pH)$p.value #normal, 0,71778
plot(AmbNum2$Conductivity)
shapiro.test(AmbNum2$Conductivity)$p.value #no, 3,394 e-15
plot(AmbNum2$Turbidity)
shapiro.test(AmbNum2$Turbidity)$p.value #no, 6,695 e-05
plot(AmbNum2$TON)
shapiro.test(AmbNum2$TON)$p.value #no, 0,000114
plot(AmbNum2$TP)
shapiro.test(AmbNum2$TP)$p.value #no, 1,489 e-10
plot(AmbNum2$TDP)
shapiro.test(AmbNum2$TDP)$p.value #no, 2,806 e-10
plot(AmbNum2$Chl.a)
shapiro.test(AmbNum2$Chl.a)$p.value #no, 5,089 e-12
plot(AmbNum2$DOC)
shapiro.test(AmbNum2$DOC)$p.value #casi, 0,04999
plot(AmbNum2$SMN.Anual.Precipitation)
shapiro.test(AmbNum2$SMN.Anual.Precipitation)$p.value #no, 2,309 e-06
plot(AmbNum2$Cultivo)
shapiro.test(AmbNum2$Cultivo)$p.value #no, 9,17 e-06
plot(AmbNum2$Doble.Cultivo)
shapiro.test(AmbNum2$Doble.Cultivo)$p.value #no, 6,399 e-08
plot(AmbNum2$Pasturas)
shapiro.test(AmbNum2$Pasturas)$p.value #no, 1,082 e-06
plot(AmbNum2$Cultivo.de.verano)
shapiro.test(AmbNum2$Cultivo.de.verano)$p.value #no, 1,918 e-06
plot(AmbNum2$Pastizalesnaturales)
shapiro.test(AmbNum2$Pastizalesnaturales)$p.value #no, 1,96 e-05
plot(AmbNum2$Wetland)
shapiro.test(AmbNum2$Wetland)$p.value #no, 1,766 e-06
plot(AmbNum2$Agua)
shapiro.test(AmbNum2$Agua)$p.value #no, 1,66 e-06
plot(AmbNum2$urbanizado) 
shapiro.test(AmbNum2$urbanizado)$p.value #no, 1,64 e-10

AmbEscalado <- scale(AmbNum2, center=TRUE , scale=TRUE)
AmbEscalado <- as.data.frame(AmbEscalado)

#LEER PARA VER CORREL: https://rpubs.com/pjmurphy/338798

plot(AmbEscalado$Latitude)
shapiro.test(AmbEscalado$Latitude)$p.value #muy al limite de 0,05 (da 0,0571)
plot(AmbEscalado$Longitude)
shapiro.test(AmbEscalado$Longitude)$p.value #no, 0,000345
plot(AmbEscalado$Lake.Area.km2)
shapiro.test(AmbEscalado$Lake.Area.km2)$p.value #no, 7,48 xe-10
plot(AmbEscalado$C.D.)
shapiro.test(AmbEscalado$C.D.)$p.value #no, 1,198 xe-05
plot(AmbEscalado$Water.Temperature)
shapiro.test(AmbEscalado$Water.Temperature)$p.value #normal, 0,1095
plot(AmbEscalado$SD)
shapiro.test(AmbEscalado$SD)$p.value #no, 1,861 e-11
plot(AmbEscalado$DO)
shapiro.test(AmbEscalado$DO)$p.value #no, 0,000823
plot(AmbEscalado$pH)
shapiro.test(AmbEscalado$pH)$p.value #normal, 0,71778
plot(AmbEscalado$Conductivity)
shapiro.test(AmbEscalado$Conductivity)$p.value #no, 3,394 e-15
plot(AmbEscalado$Turbidity)
shapiro.test(AmbEscalado$Turbidity)$p.value #no, 6,695 e-05
plot(AmbEscalado$TON)
shapiro.test(AmbEscalado$TON)$p.value #no, 0,000114
plot(AmbEscalado$TP)
shapiro.test(AmbEscalado$TP)$p.value #no, 1,489 e-10
plot(AmbEscalado$TDP)
shapiro.test(AmbEscalado$TDP)$p.value #no, 2,806 e-10
plot(AmbEscalado$Chl.a)
shapiro.test(AmbEscalado$Chl.a)$p.value #no, 5,089 e-12
plot(AmbEscalado$DOC)
shapiro.test(AmbEscalado$DOC)$p.value #casi, 0,04999
plot(AmbEscalado$SMN.Anual.Precipitation)
shapiro.test(AmbEscalado$SMN.Anual.Precipitation)$p.value #no, 2,309 e-06
plot(AmbEscalado$Cultivo)
#ver si estas variables se pueden poner como dummys o proporciones
#https://www.dummies.com/programming/r/how-to-calculate-data-proportions-and-find-the-center-in-r/
shapiro.test(AmbEscalado$Cultivo)$p.value #no, 9,17 e-06
plot(AmbEscalado$Doble.Cultivo)
shapiro.test(AmbEscalado$Doble.Cultivo)$p.value #no, 6,399 e-08
plot(AmbEscalado$Pasturas)
shapiro.test(AmbEscalado$Pasturas)$p.value #no, 1,082 e-06
plot(AmbEscalado$Cultivo.de.verano)
shapiro.test(AmbEscalado$Cultivo.de.verano)$p.value #no, 1,918 e-06
plot(AmbEscalado$Pastizalesnaturales)
shapiro.test(AmbEscalado$Pastizalesnaturales)$p.value #no, 1,96 e-05
plot(AmbEscalado$Wetland)
shapiro.test(AmbEscalado$Wetland)$p.value #no, 1,766 e-06
plot(AmbEscalado$Agua)
shapiro.test(AmbEscalado$Agua)$p.value #no, 1,66 e-06
plot(AmbEscalado$urbanizado) 
shapiro.test(AmbEscalado$urbanizado)$p.value #no, 1,64 e-10


#transformacion de matriz Otu inicial################################################
#realizar la transformacion de Hellinger a Otus como data.frame
OtusH <- hellinger(Otus)
# Transpone todas las columnas copiando la primera de molde
OtusTras <- data.frame(t(OtusH[])) #era Otus 
# Añadimos los nombres de las columnas
colnames(OtusTras) <- Otus[ ,1]
colnames(OtusTras) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52")
OtusTras
#como le borro las filas NA???
OtusTras <- na.omit(OtusTras)
#por esto la traspuse! Sino tengo demasiados NA y no puedo correr el CCA
rowSums(OtusH)
colSums(OtusH)


#Ahora pruebo con RDA #######################################################################
RDA<- rda(X, Y, scale=FALSE)
RDA
summary(RDA) #sacar los scores para las distancias euclideanas de cada OTU
# R^2 Ajustado y estandard
RsquareAdj(RDA)

#grafico los componentes principales del CCA
barplotr<-screeplot(RDA, bstick = FALSE, type = c("barplot", "lines"), npcs = min(10, if (is.null(RDA$RDA) || RDA$RDA$rank == 0) RDA$CA$rank else RDA$RDA$rank), ptype = "o", bst.col = "red", bst.lty = "solid", xlab = "Component", ylab = "Inertia", main = deparse(substitute(RDA)), legend = bstick)
barplotr2<- screeplot(RDA, bstick = FALSE, type = c("barplot", "lines"), npcs = 4, ptype = "o", bst.col = "red", bst.lty = "solid", xlab = "Component", ylab = "Inertia", main = deparse(substitute(RDA)))
#y despues, los ploteo en nube de puntos
RDA.plotr1 <- plot(RDA,choices=c(1,2), xlim=c(-1,1), ylim=c(-1,1), display=c("sp","cn"), main = "Biplot RDA Otus") 
RDA.plotr2 <- plot(RDA,choices=c(1,2), xlim=c(-1,1), ylim=c(-1,1), display=c("lc","cn"), main = "Biplot RDA Lagunas")
RDA.plotr3 <- plot(RDA,choices=c(1,2), xlim=c(-1,1), ylim=c(-1,1), display=c("sp","lc"), main = "Biplot RDA OtusxLagunas") 

#prueba de hipotesis y p-valor
RDA.ANOVA <- anova (RDA, permutations = how(nperm = 999))
RDA.ANOVA

scores(RDA, display=c("lc","sp"))
scores(RDA.plot3$sp, display=c("sp"))
write.csv(RCA.plot3$species,file.path(getwd(),"OTUscoresHRDA.csv"))


VIFRDA <- vif.cca(RDA)
VIFRDA
summary(VIFRDA)


#RC con rarefaccionados ###############################################
#ahora lo mismo pero rarefaccionados:
Oturaretras <- as.data.frame(t(Otusrarefaccionados)) #mismo planteo, es por OTU o por sitio?
spXsiteRARE <- Oturaretras

comm=spXsiteRARE
rand.time=100 # usually use 1000 for real data.
nworker=2 # parallel computing thread number
RC.RARE=RC.pc(comm=comm, rand = rand.time,
              nworker = nworker, weighted = TRUE,
              sig.index="RC")

RC.RARE
write.table(RC.RARE,file="Raup-Crick-RARE.txt")
write.csv(RC.RARE,file.path(getwd(),"Raup-Crick-RARE.csv"))


pdf("Raup-Crick index.pdf")
hist(RC$index)