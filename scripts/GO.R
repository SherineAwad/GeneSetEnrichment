library("edgeR")
library(biomaRt)
library(dplyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(data.table)
library("gplots")
library(magrittr)
library(pathview)
library("org.Mm.eg.db")
library(gage)
library(gageData)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
dge = args[1]
ORGANISM  =args[2] 
compare_type = args[3]
outname = args[4]
 
res <- read.csv(file = dge, header = TRUE)
summary(res) 
foldchanges = res$logFC
names(foldchanges) = res$entrez
outname = paste(outname, compare_type, sep ="-")
outGO = paste(outname, "GO.csv", sep ="-")

if (ORGANISM == "HUMAN")
{
data(kegg.sets.hs)
data(go.sets.hs)
data(carta.hs)
data(sigmet.idx.hs)
data(go.subs.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
#GO
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=kegg.sets.hs, same.dir =TRUE, compare =compare_type,cutoff=0.05)
lapply(gobpres, head)
write.csv(gobpres, outGO)

} else if (ORGANISM == "MOUSE"){
data(kegg.sets.mm)
data(go.sets.mm)
data(carta.mm)
data(sigmet.idx.mm)
data(go.subs.mm)


kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm,3)
data(go.sets.mm)
data(go.subs.mm)
gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=kegg.sets.mm, same.dir =TRUE, compare =compare_type,cutoff=0.05)
lapply(gobpres, head)
write.csv(gobpres, outGO)
} 

