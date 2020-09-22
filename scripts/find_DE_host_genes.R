#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Find baboon genes with expression correlated with parasitemia
# ----------------------------------------------------------------------------------------

# cd /scratch/aet359/colobus_hep

module load r/intel/3.6.0
R

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

#library(gprofiler2)    

library("GO.db")

# --- Load feature count data

load("hepato.ind.fc.Rdata")

fc = hepato.ind.fc

# Show number of genes and individuals
dim(fc$counts)

# --- Load parasitemia info

parasitemia = read.table("parasitemia_proxy.txt", header=TRUE)

# Tally things up since numbers are per-lane

parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping =
    as.numeric(gsub(",", "", parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping))
parasitemia$Hepatocystis_Reads_Mapped =
    as.numeric(gsub(",", "", parasitemia$Hepatocystis_Reads_Mapped))

parasitemia.info = aggregate(parasitemia[,c("BioSample", "Sex")],
    by=list(parasitemia$Sample_name), FUN=unique)

parasitemia.sum = aggregate(parasitemia[,c("Colobus_Reads_Mapped_Hepatocystis_mapping", "Hepatocystis_Reads_Mapped")],
    by=list(parasitemia$Sample_name), FUN=sum)

parasitemia = merge(parasitemia.info, parasitemia.sum)
names(parasitemia)[1] = "Sample_name"

parasitemia$parasitemia.proxy =
    parasitemia$Hepatocystis_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts=fc$counts, genes=fc$annotation)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge) > 5) >= 2

# Make reference table list
ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))

write.table(ref.list, file="results/reference_list.colobus.txt",
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
table(gsub(".hepato.bam", "", row.names(fc.dge.norm$samples)) == parasitemia$Sample)
#  TRUE
#   29

design = model.matrix(~ parasitemia.proxy, data=parasitemia)

# Removed: parasitemia$Sex    

#rownames(design) = rownames(fc.dge.norm$samples)

# colnames(design)

disp = estimateDisp(fc.dge.norm, design, robust = TRUE)

# CPM
cpm.disp = cpm(disp)

# LRT method ---------------------------------------------------------------
# fit = glmFit(disp, design, robust = TRUE)
# lrt = glmLRT(fit, coef = "parasitemia.proxy")
# tt  = topTags(lrt, n=Inf, adjust.method = "BH", p.value = 1)

# Or to add a log fold change cutoff ---------------------------------------------------------------   
# fit = glmFit(disp, design, robust = TRUE)
# tr = glmTreat(fit, coef = "parasitemia.proxy", lfc=1)
# tt = topTags(tr, n=Inf, adjust.method = "BH", p.value = 1)

# Or QLFTest method ---------------------------------------------------------------   
fit = glmQLFit(disp, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = "parasitemia.proxy")
tt  = topTags(qlf, n=Inf, adjust.method = "BH", p.value = 1)

    

# Write list of significant up- and down-regulated genes
sig.genes.up = tt$table$GeneID[tt$table$FDR < 0.1 & tt$table$logFC > 0]
length(sig.genes.up)

sig.genes.dn = tt$table$GeneID[tt$table$FDR < 0.1 & tt$table$logFC < 0]
length(sig.genes.dn)

# Write background list of genes
bg.genes = tt$table$GeneID
length(bg.genes)

# Write results to file
write.table(sig.genes.up, file="reports/edgeR_colobus_upreg_sig_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(sig.genes.dn, file="reports/edgeR_colobus_dnreg_sig_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(bg.genes, file="reports/edgeR_colobus_bg_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
