library("phangorn")
library("dada2")
library("phyloseq")

error

setwd("D://Romy//FACU//doctoral//datos//52 lagunas//estudio evolutivo//arboles obtenidos")
treeopt <- readRDS("fitGTR12.rds")
treeraw <- read.nexus("arbolNJ52.nex")

plot_tree(treeraw, method = "sampledodge", nodelabf = NULL,
          color = NULL, shape = NULL, size = NULL, min.abundance = Inf,
          label.tips = NULL, text.size = NULL, sizebase = 5,
          base.spacing = 0.02, ladderize = FALSE, plot.margin = 0.2,
          title = NULL, treetheme = NULL, justify = "jagged")

#limpio de quimeras###################################################################
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, minFoldParentOverAbundance = 8)
dim(seqtab_nochim)

#los clasifico taxonomicamente########################################################
tax <- assignTaxonomy(seqtab_nochim, "D:/Romy/FACU/doctoral/datos/52 lagunas/Estudio secuencias/secuencias/Analisis Maria/entrenamiento v132/silva_nr_v132_train_set.fa", multithread=TRUE)
tax_add <- addSpecies(tax, "D:/Romy/FACU/doctoral/datos/52 lagunas/Estudio secuencias/secuencias/Analisis Maria/entrenamiento v132/silva_species_assignment_v132.fa")

#If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. If using DECIPHER for taxonomy, try IdTaxa (..., strand="both").

#salvamos ************************************************************************
saveRDS(seqtab_nochim, file.path(getwd(),"seqtab_final.rds"))
saveRDS(tax, file.path(getwd(),"tax_final.rds"))
saveRDS(tax_add, file.path(getwd(),"tax_add_final.rds"))

######################################################################################
#otra forma de clasificar:
#usando el paquete DECIPHER
library(DECIPHER); packageVersion("DECIPHER")
#bajar su comprimido de entrenamiento de la pag " http://DECIPHER.codes/Downloads.html. Download the SILVA SSU r132 (modified)"

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#######################################################################################
#armo el objeto phyloseq final

Otus<-read.csv("OTU_Table_R.csv", header = T, row.names = 1)
Tax<-read.csv("Taxonomy.csv", header = T, row.names = 1)
TaxM<-as.matrix(Tax)
OTU=otu_table(Otus, taxa_are_rows = TRUE)
TAX=tax_table(TaxM)
common.ids <- intersect(rownames(TAX), colnames(OTU))

#ps = phyloseq (OTU, TAX) #este es el que funciono

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample



add.tips(tree, tips, where, edge.length = NULL)
#Arguments
#tree an object of class "phylo".
#tips a character vector containing the names of the tips.
#where an integer or character vector of the same length as tips giving the number of the
#node or tip of the tree x where the tree y is binded.
#edge.length optional numeric vector with edge length



