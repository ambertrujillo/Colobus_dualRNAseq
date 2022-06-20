#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Find Hepatocystis genes with expression correlated with parasitemia
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

# To install BioConductor packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("GO.db")

library(limma)
library(edgeR)

library(ggplot2)
library(ggrepel)

library(plyr)
library(xtable)

library(gprofiler2)    

library("GO.db")

# --- Load feature count data

load("edgeR_results/hepato.fc.Rdata")

fc = hepato.fc

# Show number of genes and individuals
dim(fc$counts)
# 45274    29

# --- Load parasitemia info

parasitemia = read.csv("parasitemia_proxy.csv", header=TRUE)

# Tally things up since numbers are per-lane

parasitemia$Colobus_Reads_Mapped=
    as.numeric(gsub(",", "", parasitemia$Colobus_Reads_gene))
parasitemia$Hepatocystis_Reads_Mapped =
    as.numeric(gsub(",", "", parasitemia$Hepatocystis_Reads_gene))

parasitemia.info = aggregate(parasitemia[,c("Sample_name", "Sex")],
    by=list(parasitemia$Sample_name), FUN=unique)

parasitemia.sum = aggregate(parasitemia[,c("Colobus_Reads_Mapped", "Hepatocystis_Reads_Mapped")],
    by=list(parasitemia$Sample_name), FUN=sum)

parasitemia = merge(parasitemia.info, parasitemia.sum)
names(parasitemia)[1] = "Sample_name"

parasitemia$parasitemia.proxy =
    parasitemia$Hepatocystis_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts=fc$counts, genes=fc$annotation)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge) > 5) >= 2

# Make reference table list
ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))
nrow(ref.list)
#4799

write.table(ref.list, file="edgeR_results/reference_list.hepato.txt",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Get CPM just normalized by library composition
cpm = cpm(fc.dge.norm)

# --- Estimate dispersion the complicated way (using CR method)

# Ensure samples are in correct order
# (i.e. make sure fc.dge.norm$samples and parasitemia$Individual have same name)
table(gsub(".hepato.bam", "_001", row.names(fc.dge.norm$samples)) == parasitemia$Sample_name)
#  TRUE
#   29

design = model.matrix(~ parasitemia.proxy, data=parasitemia)

# Removed: parasitemia$Sex    

#rownames(design) = rownames(fc.dge.norm$samples)

# colnames(design)

disp = estimateDisp(fc.dge.norm, design, robust = TRUE)

# CPM
cpm.disp = cpm(disp)

# Write list of significant up- and down-regulated genes
sig.genes.up = subset(tt$table, FDR < 0.2 & logFC > 0, select=c(logFC, FDR))
nrow(sig.genes.up)
#3

sig.genes.dn = subset(tt$table, FDR < 0.2 & logFC < 0, select=c(logFC, FDR))
nrow(sig.genes.dn)
#38

# Write background list of genes
bg.genes = tt$table$GeneID
length(bg.genes)
# 4799

# Write results to file
write.table(sig.genes.up, file="edgeR_results/hepato_upreg_sig_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(sig.genes.dn, file="edgeR_results/hepato_dnreg_sig_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

write.table(bg.genes, file="edgeR_results/hepato_bg_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
