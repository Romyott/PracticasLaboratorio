rm(list=ls())
setwd("D://Romy//FACU//doctoral//datos//52 lagunas")
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(data.table)
library("gplots")
library(stats)
Amb<-read.csv("LagunasPampeanasArroyos.csv", header=TRUE, row.names=1, sep = ";" , dec=",")
Otus<-read.csv("Core Otus.csv", header = T, row.names = 1)
Tax<-read.csv("Core Taxonomy.csv", header = T, row.names = 1)

#Como phyloseq hace lio con tabla taxoncomica, la transformo en matriZ
TaxM<-as.matrix(Tax)
#puedo subir el arból y luego también agregarlo como objeto phyloseq
#CREACION DE OBJETOS PHYLOSEQ DE FORMA MANUAL:
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(TaxM)
AMB=sample_data(as.data.frame(Amb)) #sin esto funciona el ps :/
rownames(AMB) <- c("SE","KH","SY","SG","HR","MA2","LP","BV","CH","LI","TR","SM","PI","FS","FN","CS","SU","TI","LU","LJ","SP","ST","EC","LB","VI","TA","BS","CC","BU","AD","CO","MO","VE","EP","PA","CT","PU","AL","BG","LO","RC","MA1","GO","CA","AZ","BR","ZO","QU","HG","PE","IM","MT")

common.ids <- intersect(rownames(TAX), colnames(OTU))

ps1 <- phyloseq(OTU, TAX, AMB)
ps1

heatmap(otu_table(ps1))
plot_heatmap(ps1)
plot_heatmap(ps1, species.label = "Family")
plot_heatmap(ps1, method = "DCA", distance = "bray",
             sample.label = NULL, taxa.label = "Phylum", low = "deepskyblue2",
             high = "brown2", na.value = "black", max.label = 250, title = "Heatmap Taxa-Laguna", taxa.order = "Phylum")

#https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
#https://joey711.github.io/phyloseq/plot_heatmap-examples.html
                                                  