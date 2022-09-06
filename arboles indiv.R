library("treeio")
library("ggplot2")
library("phyloseq")
library("dplyr")
library("dada2")
library("phangorn")
library("corrplot")
library("DECIPHER")
library(dada2); packageVersion("dada2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("dada2")

#remove.packages("dada2")


setwd("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo")
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//errores.rdata")

pathFX <- file.path(getwd(), "pathFX")
pathRX <- file.path(getwd(), "pathRX")

filtpathF <- file.path(pathFX, "filtered") 
filtpathR <- file.path(pathRX, "filtered") 

fastqFs <- sort(list.files(pathFX, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathRX, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#trimeo segun cada laguna, el truncLen y el trimLeft van a variar **************** 
#out <- filterAndTrim(fwd=file.path(pathFX, fastqFs), filt=file.path(filtpathF, fastqFs),
#                     rev=file.path(pathRX, fastqRs), filt.rev=file.path(filtpathR, fastqRs),trimLeft = 25,
#                     truncLen=c(300,190), maxEE=c(2,5), maxN = 0, truncQ = 3, rm.phix=TRUE,
#                     compress=TRUE, multithread=FALSE)

#head(out)


###################################################################################
filtpathF <- file.path(getwd(), "pathFX","filtered")
filtpathR <- file.path(getwd(), "pathRX","filtered")

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R1"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_R2"), `[`, 1) 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)


################################################################################
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
dadaFs <- vector("list", length(sample.names))
names(dadaFs) <- sample.names
dadaRs <- vector("list", length(sample.names))
names(dadaRs) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  dadaFs[[sam]]<-ddF
  dadaFs[[sam]]
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  dadaRs[[sam]]<-ddR
  dadaFs[[sam]]
  merger <- mergePairs(ddF, derepF, ddR, derepR, maxMismatch=1)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)


######################################################################################
#Tabla de secuencias y les sacamos las quimeras
seqtab <- makeSequenceTable(mergers)
#saveRDS(seqtab, file.path(getwd(),"seqtab-ARROYOCHASICO.rds")) #*************************************
str(seqtab)
#Sacamos quimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, minFoldParentOverAbundance = 8)
dim(seqtab_nochim)
paste("% of non removed reads = ", sum(seqtab_nochim)/sum(seqtab)*100, sep="")
#pdf(file.path(getwd(),"seqlen-RIOSALADO.pdf")) #*************************************************
plot(table(nchar(getSequences(seqtab_nochim))))
dev.off()

#saveRDS(seqtab_nochim, file.path(getwd(),"seqtab_nochim-RIOSALADO.rds")) #***********************


######################################################################
#agregar la taxonomia obtenida
#tax <- assignTaxonomy(seqtab_nochim, "D:/Romy/FACU/doctoral/datos/52 lagunas/estudio evolutivo/entrenamiento v132/silva_nr_v132_train_set.fa", multithread=TRUE)
#tax_add <- addSpecies(tax, "D:/Romy/FACU/doctoral/datos/52 lagunas/Estudio evolutivo/entrenamiento v132/silva_species_assignment_v132.fa")

#saveRDS(tax, file.path(getwd(),"Tax-RIOSALADO.rds")) #lo salvamos porlas***************************************
#write.csv(tax,file.path(getwd(),"Tax-RIOSALADO.csv")) #lo salvamos porlas*************************************


######################################################################
seqs <- getSequences(seqtab_nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#saveRDS(alignment, file.path(getwd(),"alignment-RIOSALADO.rds")) #***************************


#########################################################################
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session-RIOSALADOalineacion.rdata")
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session-RIOSALADOalineacion.rdata")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "ratchet", control = pml.control(trace = 0)) #probar otros modelos evol y otros rearreglos, a ver que obtengo...


plot_tree(fitGTR$tree, method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = NULL, text.size = 3, sizebase = 5,
          base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
          title = NULL, treetheme = NULL, justify = "jagged") #REVISAR ACA PORLAS EL PLOTEO


#ahora guardo el arbol, preferentemente en formato .nexus o .tree
#Cambio la X por el numero de laguna que corresponda *********************************
#ape::write.nexus(fitGTR$tree, file='arbolGTR-ARROYOQUEQUEN.nex')
#save(fitGTR, file='fitGTR-ARROYOQUEQUEN.pml')
#saveRDS(fitGTR, file.path(getwd(),"fitGTR-ARROYOQUEQUEN.rds"))

#guardo la sesion tal cual la deje###############################################
#cambiar la X por quien corresponda **************************************************
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session-RIOSALADO.rdata")









###############################################################################
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session18.rdata")
#!!!!!!!!!!!!!!!!REPETIR ESTA PARTE CON RIO ARRECIFES PORQUE LO SOBREESCRIBI CON LO DE SALADO
#ponerle raiz al arbol
set.seed(800)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

#armo el objeto phyloseq, pero los tips son las secuencuas
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), #
               tax_table(tax),phy_tree(fitGTR$tree)) 

#asi que le tengo que poner nombre a las puntas
#unir taxa y otu en una tabla :O
#https://stackoverflow.com/questions/53032504/combine-otu-and-tax-table-and-replace-actual-sequences-with-otu-ids-phyloseq-da
#https://github.com/joey711/phyloseq/issues/213
#https://github.com/joey711/phyloseq/issues/833
#https://benjjneb.github.io/dada2/assign.html

taxa_names(ps)
n_seqs <- seq(ntaxa(ps)) #
len_n_seqs <- nchar(max(n_seqs)) 
taxa_names(ps) <- paste("Seq", formatC(n_seqs, 
                                       width = len_n_seqs, 
                                       flag = "0"), sep = "_")
taxa_names(ps)

#generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 
                             sep = "__"))  # to distinguish from "_" within tax ranks

#turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_table(ps))
tmp <- names(otu_export)

#paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

#overwrite old names
names(otu_export) <- names(tmp)
head(otu_export)[5]

taxa_names(ps) <- paste0("OTU_", seq(ntaxa(ps))) #y con este comando le pongo el nombre artificial que yo quiera a los labels del arbol

plot_tree(phy_tree(ps), method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = "taxa_names", text.size = 3, sizebase = 5, #perra taxonomia...
          base.spacing = 0.2, ladderize = FALSE, plot.margin = 0.5,
          title = NULL, treetheme = NULL, justify = "jagged") 



tax_table(ps)

saveRDS(tax_table(ps), file.path(getwd(),"TaxaOTUarbol-ARROYOQUEQUENSALADO.rds")) #***************************************
write.csv(tax_table(ps),file.path(getwd(),"TaxaOTUarbol-ARROYOQUEQUENSALADO.csv")) #*************************************
saveRDS(phy_tree(ps), file.path(getwd(),"Arbol-ARROYOQUEQUENSALADO.rds")) #***************************************
ape::write.nexus(phy_tree(ps), file='arbol-ARROYOQUEQUENSALADO.nex') #**************************************
saveRDS(ps, file.path(getwd(),"psArbol-ARROYOQUEQUENSALADO.rds")) #***************************************
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session-ARROYOQUEQUENSALADOb.rdata")






#save.image(phy_tree(ps), file="graficoarbolXX.jpg")
#guardo la sesion tal cual la deje###############################################
#cambiar la X por quien corresponda **************************************************








load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session10.rdata")
########################################################################################
#armo el pbjeto phyloseq
Otus<-read.csv("OTU_Table_R.csv", header = T, row.names = 1)
Tax<- tax     #read.csv("Tax02.csv", header = T, row.names = 1) #**********************************
#Amb<-read.csv("LagunasPampeanas.csv", header=TRUE, row.names=1, sep = ";" , dec=",") #aca no tiene sentido porque es una unica laguna
Tree<-fitGTR$tree

#TaxM<-as.matrix(Tax) #no hace falta, tax ya es matriz
OtuM <-as.matrix(Otus)

#armo el objeto phyloseq con otu table, las secuencias sin quimeras, la matrix tutoreada de taxones y el arbol optimizado
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(Tax)
AMB=sample_data(as.data.frame(Amb)) #sin esto funciona el ps :/


ps <-phyloseq(OTU,TAX,phy_tree(fitGTR$tree)) #probanding




ps



saveRDS(ps, file.path(getwd(),"ps10.rds")) #**********************************************

plot_tree(phy_tree(ps), method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = "ps@tax_table", text.size = 3, sizebase = 5, #perra taxonomia...
          base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
          title = NULL, treetheme = NULL, justify = "jagged") #REV
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//psdelorto.rdata")



#########################################################

saveRDS(otu_export, file.path(getwd(),"otu asignado.rds")) #lo salvamos porlas***************************************
write.csv(otu_export,file.path(getwd(),"otu asignado.csv")) #lo salvamos porlas*************************************


#funciona!!!
###################################################################################
#ahora hacerlo funcionar en ps lpm!!!!





tax_table(ps)

OtuM <-as.matrix(otu_export)

tax_table <- cbind(tax_table, Strain=otu_export) #esto funca pero... me queda medio rari... cada otu tiene una seq y un taxon asociados

ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), #aca iba FALSE
               tax_table(tax),phy_tree(fitGTR$tree)) #asi queda pero los tips son las sec, no los taxones

taxa_names(ps) <- paste0("OTU", seq(ntaxa(ps)))

plot_tree(phy_tree(ps), method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = "taxa_names", text.size = 3, sizebase = 5, #perra taxonomia...
          base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
          title = NULL, treetheme = NULL, justify = "jagged") #REV

#################################
#otra manera
short_otu = head(otu_table(ps))
short_tax = tax_table(ps)[rownames(short_otu), ]
substr
dput

###################################################################################
#de aca en adelante puedo correr el programa de phyloseq sin problema (EstudioOTUstaller.R)



#guardo la sesion tal cual la deje###############################################
#cambiar la X por quien corresponda **************************************************
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_session01.rdata")




# Close R, Re-open R
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//my_sessionX.rdata")


