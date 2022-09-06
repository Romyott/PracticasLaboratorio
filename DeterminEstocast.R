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
#install.packages("iCAMP")
library("iCAMP")

setwd("D://Romy//FACU//doctoral//datos//52 lagunas//Deterministico Estocastico")
#load("D://Romy//FACU//doctoral//datos//52 lagunas//Deterministico Estocastico//sesion-ENCADENADASOESTE.rdata")

Otus<-read.csv("Core Otus.csv", header = T, row.names = 1)
Tree <- read.nexus ("arbol-RIOSAUCEGRANDE.nex")#("ArbolTotaltips.nex")#********************************
TaxTree<-read.csv ("TaxaOTUarbol-RIOSAUCEGRANDE.csv", header= T, row.names = 1) #("TaxaArbol.csv", header = T, row.names = 1)

TaxM<-as.matrix(TaxTree)
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(TaxM)


#ps##########################################################################################
ps1 = phyloseq(OTU, Tree, TAX)
tax_table(ps1)
otu_table(ps1)
phy_tree(ps1)
#saco la distancia cofenetica del arbolito
#phydist1<-cophenetic(phy_tree(ps1))
#write.csv(phydist1,file.path(getwd(),"DistanciaFilogenetica-RIOSALADO.csv"))#*****************

#preparo el ps para RC y bNTI ###########################################################
#eliminar singletons
psSinSigles <- filter_taxa(otu_table(ps1), function (x) {sum(x > 0) > 1}, prune=TRUE)
psSinSigles

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
#write.csv(beta.mntd.weighted,'betaMNTD_weighted-psRARE.csv',quote=F);

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
hist(weighted.bNTI)

#write.csv(weighted.bNTI,file.path(getwd(),"bNTI-RIOSALADO.csv"))


#Raup-Ckrick######################################################################################
comm= t(otu_table(psSinSigles)) #si pongo la matriz traspuesta, toma las lagunas. si la pongo normal, toma por taxon
rand.time=500 # usually use 1000 for real data.
nworker=2 # parallel computing thread number
RC=RC.pc(comm=comm, rand = rand.time,
         nworker = nworker, weighted = TRUE,
         sig.index="RC")

RCM <- as.matrix(RC$index)
RCM

hist(RC$index)
#write.csv(RC,file.path(getwd(),"Raup-Crick-RIOSALADO.csv")) #**************************
#write.csv(RCM,file.path(getwd(),"Raup-Crick-RIOSALADO-MATRIZ.csv")) #******************
###################################################################################
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//Deterministico Estocastico//sesion-RIOSAUCEGRANDE.rdata")

#conteo#################################################################################
#conteos (hasta que arregle las matrices)
install.packages("tidyverse")
library("tidyverse")

RCM[lower.tri(RCM, diag = TRUE)] <- NA
RCMm <- as.data.frame(RCM)

conteo0 <- RCMm %>% tally(RCMm > -10)
conteo1 <- RCMm %>% tally(RCMm > 0.95)
conteo2 <- RCMm %>% tally(RCMm < -0.95)


BETAdata <- as.data.frame(weighted.bNTI)
conteoA <- BETAdata %>% tally(BETAdata > -10)
conteoB <- BETAdata %>% tally(BETAdata > 2)
conteoC <- BETAdata %>% tally(BETAdata < -2)




