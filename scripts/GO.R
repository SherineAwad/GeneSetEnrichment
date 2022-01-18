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
outname = args[3]
 
res <- read.csv(file = dge, header = TRUE)
summary(res) 
#---------------------------
#KEGG Analysis 
#---------------------------
if (ORGANISM == "HUMAN")
{
data(kegg.sets.hs)
data(go.sets.hs)
data(carta.hs)
data(sigmet.idx.hs)
data(go.subs.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
foldchanges = res$logFC
names(foldchanges) = res$entrez
#GO
outGO = paste(outname, "GO.csv", sep ="-")
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=kegg.sets.hs, same.dir =TRUE, compare ="unpaired",cutoff=0.05)
lapply(gobpres, head)
write.csv(gobpres, outGO)

} else if (ORGANISM == "MOUSE"){
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

data(kegg.sets.mm)
data(go.sets.mm)
data(carta.mm)
data(sigmet.idx.mm)
data(go.subs.mm)


kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm,3)
foldchanges = res$logFC
names(foldchanges) = res$entrez

#GO
outGO = paste(outname, "GO.csv", sep ="-")
data(go.sets.mm)
data(go.subs.mm)
gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=kegg.sets.mm, same.dir =TRUE, compare ="unpaired",cutoff=0.05)
lapply(gobpres, head)
write.csv(gobpres, outGO)
} 

