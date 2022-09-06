library("treeio")
library("ggplot2")
library("phyloseq")
library("dplyr")
library("dada2")
library("phangorn")
library("corrplot")
library("DECIPHER")

#remove.packages("vctrs") 
#install.packages("vctrs")
########################################################
setwd("D://Romy//FACU//doctoral//datos//52 lagunas//estudio secuencias//secuencias//Analisis Maria")
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//errores.rdata")
load("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//arbolcompleto.rdata")


pathF <- file.path(getwd(), "pathF")
pathR <- file.path(getwd(), "pathR")

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#################################################################################
out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                      rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),trimLeft = 30,
                      truncLen=c(330,210), maxEE=c(2,5), maxN = 0, truncQ = 3, rm.phix=TRUE,
                      compress=TRUE, multithread=FALSE)
head(out)

## seteo la carpeta con las sec
filtpathF <- file.path(getwd(), "pathF","filtered")
filtpathR <- file.path(getwd(), "pathR","filtered")

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

## Important: adjust characters that need to act as delimiters

sample.names <- sapply(strsplit(basename(filtFs), "_R1"), `[`, 1) # Fix according to filename
sample.namesR <- sapply(strsplit(basename(filtRs), "_R2"), `[`, 1) # Fix according to filename
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, multithread=TRUE) #nbases defoult 1e+8
# Learn reverse error rates
errR <- learnErrors(filtRs, multithread=TRUE) #nbases defoult 1e+8


save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//errores.rdata")

# Sample inference and merger of paired-end reads
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

###################################################################################
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(getwd(),"seqtab.rds"))
str(seqtab)
# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, minFoldParentOverAbundance = 8)
dim(seqtab_nochim)
paste("% of non removed reads = ", sum(seqtab_nochim)/sum(seqtab)*100, sep="")
pdf(file.path(getwd(),"seqlen.pdf"))
plot(table(nchar(getSequences(seqtab_nochim))))
dev.off()

save(seqtab_nochim, file.path(getwd(),"seqtab_nochim.csv")) #no funca...
saveRDS(seqtab_nochim, file.path(getwd(),"seqtab_nochim.rds"))


##################################################################################
#agregar la taxonomia obtenida
tax <- assignTaxonomy(seqtab_nochim, "D:/Romy/FACU/doctoral/datos/52 lagunas/Estudio secuencias/secuencias/Analisis Maria/entrenamiento v132/silva_nr_v132_train_set.fa", multithread=TRUE)
tax_add <- addSpecies(tax, "D:/Romy/FACU/doctoral/datos/52 lagunas/Estudio secuencias/secuencias/Analisis Maria/entrenamiento v132/silva_species_assignment_v132.fa")

##################################################################################
#alineacion de secuencias:
seqs <- getSequences(tax) #va seqtab_nochim pero voy a probar esto a ver si funciona
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

saveRDS(alignment, file.path(getwd(),"alignment.rds")) #lo salvamos porlas


###################################################################################
#luego se arma un arbol de ML para comenzar a trabajar
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#phy_tree(fitGTR$tree) #formas de mostrar el arbolito

plot_tree(fitGTR$tree, method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = NULL, text.size = NULL, sizebase = 5,
          base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
          title = NULL, treetheme = NULL, justify = "jagged")

#write.tree(treeNJ, file = "", append = FALSE, #para lograr un arbol en formato parentetico, siempre que el original sea formato phylo
           #digits = 4, tree.names = FALSE)

#print(treeNJ) #a ver como quedo el arbolito; hay que enraizarlo







#ahora guardo el arbol, preferentemente en formato .nexus o .tree

ape::write.nexus(treeNJ, file='arbolNJ.nex')
ape::write.nexus(fitGTR, file='fitGTR.nex')


write.nexus(..., file = "", translate = TRUE) #para escribir un arbol de formato phylo a nexus

saveRDS(treeNJ, file.path(getwd(),"arbolito.rds")) #salvemos todo al final, siempre
saveRDS(fitGTR, file.path(getwd(), "arbolmodelado.rds"))





####################################################################
#guardo la sesion tal cual la deje, para seguir despues
save.image("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//arbolcompleto.rdata")


#esto es para objetos RDS
saveRDS(seqtab, "path/to/seqtab.rds")
# Close and re-open R
seqtab <- readRDS("path/to/seqtab.rds")












######################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")


if (!requireNamespace("Bioconductor", quietly = TRUE))
  install.packages("Bioconductor")

BiocManager::install("decipher")

install.packages("decipher")
s
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
a

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("treeio", "ggtree"))
a
install.packages("Biostrings")
install.packages("BiocGenerics")
install.packages("parallel")
install.packages("dada2")
install.packages("phyloseq")
install.packages("ellipsis")

                                              
                                              