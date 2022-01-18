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
outname 
 
#GO
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
#---------------------------------------------------Use Kegg and gage to get upregulated and downregulated pathways
out = paste(outname, "KEGG.csv", sep="-")
outUP = paste(outname, "KEGG_UP.csv", sep ="-")
outDOWN = paste(outname, "KEGG_DOWN.csv", sep ="-")

data(kegg.gs)
keggres = gage(foldchanges, gsets =kegg.sets.hs, same.dir = TRUE, compare="paired",make.plot = TRUE)
lapply(keggres, head)
write.csv(keggres,out)
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
write.csv(keggrespathwaysup, outUP)
write.csv(keggrespathwaysdn, outDOWN)
#-------------------------------------------------------------------------------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
#gobpres = gage(foldchanges, gsets=kegg.sets.hs, same.dir =FALSE, compare ="paired",cutoff=0.05)
#lapply(gobpres, head)
#---------------------------------------------------------Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human", new.signature=FALSE)
#---------------------------------------------------------plot multiple pathways ups and downs 
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human"))
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
#---------------------------------------------------Use Kegg and gage to get upregulated and downregulated pathways
data(kegg.gs)
keggres = gage(foldchanges, gsets =kegg.sets.mm, same.dir = TRUE, compare="paired",make.plot = TRUE)
lapply(keggres, head)
write.csv(keggres,out)
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
write.csv(keggrespathwaysup, outUP)
write.csv(keggrespathwaysdn, outDOWN)
#-------------------------------------------------------------------------------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
#gobpres = gage(foldchanges, gsets=kegg.sets.mm, same.dir =FALSE, compare ="paired",cutoff=0.05)
#lapply(gobpres, head)
#---------------------------------------------------------Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse", new.signature=FALSE)
#---------------------------------------------------------plot multiple pathways ups and downs 
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse"))
} 

