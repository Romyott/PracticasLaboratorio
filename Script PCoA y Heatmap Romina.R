setwd("D://Romy//FACU//doctoral//datos//52 lagunas//Otros analisis//Va_script_para_PLS_y_yapa")
#partial least square analysis
#metodo para analisis de multicolinealidad
#install.packages("pls")
#install.packages("labdsv")
library("pls")
library("labdsv")
library(phyloseq)
require(ggplot2)
require(GGally)
require(picante)
library(ape)
library("phangorn")
#install.packages("spls")
library("spls")
#BiocManager::install("mixOmics")
library(mixOmics)
#install.packages("rlang")
library("rlang")

Otus<-read.csv("Core Otus.csv", header = T, row.names = 1)
Amb<-read.csv("LagunasPampeanasArroyos.csv", header=TRUE, row.names=1, sep = ";" , dec=",")
Tree <- read.nexus("arbol-ARROYODELOSPADRES.nex")
TaxTree<-read.csv ("TaxaOTUarbol-ARROYODELOSPADRES.csv", header = T, row.names = 1)


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
#transformo los Otus asi me da algo
TaxM<-as.matrix(TaxTree)
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(TaxM)
ps0 = phyloseq(OTU, Tree, TAX)
OtusHps <- hellinger(otu_table(ps0))
OtusHpsTras <- data.frame(t(OtusHps[])) 
OtusTras.ps <- na.omit(OtusHpsTras)

#tutorial: https://www.statology.org/partial-least-squares-in-r/

#X= OtusTras.ps
#Y= AmbEscalado 
#me gusta mas como queda con este orden
X<-AmbEscalado
Y<-OtusTras.ps
#Z = 0.5
dim(X); dim(Y)
MyResult.spls <- spls(X,Y) #,1,0.5, )
#plotIndiv(MyResult.spls)     ## sample plot
#plotVar(MyResult.spls)       ## variable plot
X11()
cim(MyResult.spls, comp = 1)


