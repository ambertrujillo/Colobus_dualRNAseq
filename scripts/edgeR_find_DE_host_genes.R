##############3 Find colobus genes with expression correlated with parasitemia ##########

# cd /scratch/aet359/colobus_hep

options(stringsAsFactors=FALSE)
options(scipen = 100)

BiocManager::install("GO.db")

library(limma)
library(edgeR)
library(statmod)

library(ggplot2)
library(ggrepel)

library(plyr)
library(xtable)

library(gprofiler2)    

library("GO.db")

################## EdgeR analysis (Colobus) ################
# --- Load feature count data

load("edgeR_results/colobus.fc.Rdata")

Colobusfc = colobus.fc

# Show number of genes and individuals
dim(Colobusfc$counts)
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
parasitemia = subset(parasitemia, select=-c(1))

parasitemia$parasitemia.proxy =
    parasitemia$Hepatocystis_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped

# --- Load PCA data (the contribution of each variable to the PC based on cell type composition. i.e., res.var$contrib)
PCA = read.table("Results_gene/PCA_results")
PCA$Sample_name = parasitemia$Sample_name
PCA_1and2 = subset(PCA, select=c(1:2, 10))
names(PCA_1and2)[1] = "PCA_1"
names(PCA_1and2)[2] = "PCA_2"

parasitemia = merge(parasitemia, PCA_1and2, by="Sample_name")

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts=Colobusfc$counts, genes=Colobusfc$annotation)


# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Get CPM just normalized by library composition
cpm = cpm(fc.dge.norm)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge.norm) > 5) >= 2

# Make reference table list
ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))

write.table(ref.list, file="Results_gene/reference_list.colobus.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

fc.dge.norm = fc.dge.norm[keep, , keep.lib.sizes=FALSE]

# --- Estimate dispersion the complicated way (using CR method)

# Ensure samples are in correct order
# (i.e. make sure fc.dge.norm$samples and parasitemia$Individual have same name)
table(gsub(".colobus.bam", "_001", row.names(fc.dge.norm$samples)) == parasitemia$Sample_name)
#  TRUE
#   29

#Model Designs
design1 = model.matrix(~ parasitemia.proxy, data=parasitemia)
design2 =  model.matrix(~ parasitemia.proxy + Sex, data=parasitemia)
design3 = model.matrix(~ parasitemia.proxy + PCA_1 + PCA_2, data=parasitemia)
#parasitemia is not capturing the immune cell composition (PCs) so they should be included in the model. Maybe add the plots of PC ~ parasitemia.proxy

#rownames(design) = rownames(fc.dge.norm$samples)

# colnames(design)

disp1 = estimateDisp(fc.dge.norm, design1, robust = TRUE)
disp2 = estimateDisp(fc.dge.norm, design2, robust = TRUE)
disp3 = estimateDisp(fc.dge.norm, design3, robust = TRUE)

# LRT method ---------------------------------------------------------------
# fit = glmFit(disp, design, robust = TRUE)
# lrt = glmLRT(fit, coef = "parasitemia.proxy")
# tt  = topTags(lrt, n=Inf, adjust.method = "BH", p.value = 1)

# Or to add a log fold change cutoff ---------------------------------------------------------------   
# fit = glmFit(disp, design, robust = TRUE)
# tr = glmTreat(fit, coef = "parasitemia.proxy", lfc=1)
# tt = topTags(tr, n=Inf, adjust.method = "BH", p.value = 1)

# Or QLFTest method ---------------------------------------------------------------   
fit1 = glmQLFit(disp1, design1, robust = TRUE)
qlf1 = glmQLFTest(fit1, coef = "parasitemia.proxy") # same as coef=2
tt1  = topTags(qlf1, n=Inf, adjust.method = "BH", p.value = 1)

fit2 = glmQLFit(disp2, design2, robust = TRUE)
qlf2 = glmQLFTest(fit2, coef = "parasitemia.proxy")
tt2  = topTags(qlf2, n=Inf, adjust.method = "BH", p.value = 1)

qlf2.sex = glmQLFTest(fit2, coef = "Sexmale")
tt2.sex  = topTags(qlf2.sex, n=Inf, adjust.method = "BH", p.value = 1)

fit3 = glmQLFit(disp3, design3, robust = TRUE)
head(coef(fit3))
qlf3 = glmQLFTest(fit3, coef = "parasitemia.proxy")
tt3 = topTags(qlf3, n=Inf, adjust.method = "BH", p.value = 1)

qlf3.PC1 = glmQLFTest(fit3, coef = "PCA_1")
tt3.PC1 = topTags(qlf3.PC1, n=Inf, adjust.method = "BH", p.value = 1)

qlf3.PC2 = glmQLFTest(fit3, coef = "PCA_2")
tt3.PC2 = topTags(qlf3.PC2, n=Inf, adjust.method = "BH", p.value = 1)


# Fix 3rd model

#try = tt3[["table"]]
#try = data.frame(tt3[["table"]])
#try$adj = p.adjust(try$PValue, method = "fdr")
#hist(tt3$table$PValue) # There is variation in the P values but after FDR correction, none of the variables are significant. Alex "too much uncertainty when interaction term is added"
# Doesnt think an interaction term is necessary

# compare sex included in model
hist(tt1$table$FDR)
hist(tt2$table$FDR)
hist(tt2.sex$table$FDR)
hist(tt3$table$FDR)
hist(tt3.PC1$table$FDR)
hist(tt3.PC2$table$FDR)

hist(qlf1$coefficients)
hist(qlf2$coefficients)
hist(qlf2.sex$coefficients)
hist(qlf3$coefficients)
hist(qlf3.PC1$coefficients)
hist(qlf3.PC2$coefficients)

#------- Write list of significant up- and down-regulated genes ------
#Up
sig.genes.up1 = tt1$table$GeneID[tt1$table$FDR < 0.2 & tt1$table$logFC > 0]
length(sig.genes.up1)
#641

hist(tt2$table$FDR, xlim = c(0,0.05))
sig.genes.up2 = tt2$table$GeneID[tt2$table$FDR < 0.2 & tt2$table$logFC > 0]
length(sig.genes.up2)
#518

sig.genes.up2.sex = tt2.sex$table$GeneID[tt2.sex$table$FDR < 0.2 & tt2.sex$table$logFC > 0]
length(sig.genes.up2.sex)
#0

hist(tt3$table$FDR, xlim = c(0,0.05))
sig.genes.up3 = tt3$table$GeneID[tt3$table$FDR < 0.2 & tt3$table$logFC > 0]
length(sig.genes.up3)
#394

sig.genes.up3.PC1 = tt3.PC1$table$GeneID[tt3.PC1$table$FDR < 0.2 & tt3.PC1$table$logFC > 0]
length(sig.genes.up3.PC1)
#2786

sig.genes.up3.PC2 = tt3.PC2$table$GeneID[tt3.PC2$table$FDR < 0.2 & tt3.PC2$table$logFC > 0]
length(sig.genes.up3.PC2)
#0

#Down
sig.genes.dn1 = tt1$table$GeneID[tt1$table$FDR < 0.2 & tt1$table$logFC < 0]
length(sig.genes.dn1)
#1497

sig.genes.dn2 = tt2$table$GeneID[tt2$table$FDR < 0.2 & tt2$table$logFC < 0]
length(sig.genes.dn2)
#1289

sig.genes.dn2.sex = tt2.sex$table$GeneID[tt2.sex$table$FDR < 0.2 & tt2.sex$table$logFC < 0]
length(sig.genes.dn2.sex)
#2

sig.genes.dn3 = tt3$table$GeneID[tt3$table$FDR < 0.2 & tt3$table$logFC < 0]
length(sig.genes.dn3)
#756

sig.genes.dn3.PC1 = tt3.PC1$table$GeneID[tt3.PC1$table$FDR < 0.2 & tt3.PC1$table$logFC < 0]
length(sig.genes.dn3.PC1)
#1918

sig.genes.dn3.PC2 = tt3.PC2$table$GeneID[tt3.PC2$table$FDR < 0.2 & tt3.PC2$table$logFC < 0]
length(sig.genes.dn3.PC2)
#0

# Write background list of genes
bg.genes = tt1$table$GeneID
length(bg.genes)
#12306

# Write results to file (for model 3 pulling out those affected by PC1)
sig.genes.up = subset(tt1$table, FDR < 0.2 & logFC > 0, select=c(logFC, FDR))
sig.genes.up2 = subset(tt2$table, FDR < 0.2 & logFC > 0, select=c(logFC, FDR))
sig.genes.up3 = subset(tt3$table, FDR < 0.2 & logFC > 0, select=c(logFC, FDR))
sig.genes.up3.PC1 = subset(tt3.PC1$table, FDR < 0.2 & logFC > 0, select=c(logFC, FDR))

sig.genes.dn = subset(tt1$table, FDR < 0.2 & logFC < 0, select=c(logFC, FDR))
sig.genes.dn2 = subset(tt2$table, FDR < 0.2 & logFC < 0, select=c(logFC, FDR))
sig.genes.dn3 = subset(tt3$table, FDR < 0.2 & logFC < 0, select=c(logFC, FDR))
sig.genes.dn3.PC1 = subset(tt3.PC1$table, FDR < 0.2 & logFC < 0, select=c(logFC, FDR))

all.sig.genes = subset(tt1$table, FDR < 0.2, select=c(logFC, FDR))
all.sig.genes2 = subset(tt2$table, FDR < 0.2, select=c(logFC, FDR))
all.sig.genes3 = subset(tt3$table, FDR < 0.2, select=c(logFC, FDR))
all.sig.genes3.PC1 = subset(tt3.PC1$table, FDR < 0.2, select=c(logFC, FDR))

bg.genes = subset(tt1$table, select=c(logFC, FDR, PValue))

write.table(sig.genes.up, file="Results_gene/colobus_upreg_sig_genes_mod1.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up2, file="Results_gene/colobus_upreg_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up3, file="Results_gene/colobus_upreg_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up3.PC1, file="Results_gene/colobus_upreg_sig_genes_mod3_PC1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)

write.table(sig.genes.dn, file="Results_gene/colobus_dnreg_sig_genes_mod1.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn2, file="Results_gene/colobus_dnreg_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn3, file="Results_gene/colobus_dnreg_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn3.PC1, file="Results_gene/colobus_dnreg_sig_genes_mod3_PC1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)

write.table(bg.genes, file="Results_gene/colobus_bg_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
    
write.table(all.sig.genes, file="Results_gene/colobus_all_sig_genes_mod1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes2, file="Results_gene/colobus_all_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes3, file="Results_gene/colobus_all_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes3.PC1, file="Results_gene/colobus_all_sig_genes_mod3_PC1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)

# --- Do functional enrichment test on GMT files - Convert genes to human ------

sig.genes.up_list=rownames(sig.genes.up)
sig.genes.up_list2=rownames(sig.genes.up2)
sig.genes.up_list3=rownames(sig.genes.up3)
sig.genes.up_list3_PC1=rownames(sig.genes.up3.PC1)

sig.genes.dn_list=rownames(sig.genes.dn)
sig.genes.dn_list2=rownames(sig.genes.dn2)
sig.genes.dn_list3=rownames(sig.genes.dn3)
sig.genes.dn_list3_PC1=rownames(sig.genes.dn3.PC1)

bg.genes_list=rownames(bg.genes)

set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15") # Sometimes it works, sometimes it doesnt, just keep trying it until it works!

# Up-regulated Model 1-3

genes.up.human = gorth(
    sig.genes.up_list,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

genes.up.human2 = gorth(
  sig.genes.up_list2,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

genes.up.human3 = gorth(
  sig.genes.up_list3,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

genes.up.human3.PC1 = gorth(
  sig.genes.up_list3_PC1,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)


genes.up.hs = genes.up.human$ortholog_ensg[genes.up.human$ortholog_ensg != "N/A"]
genes.up.hs2 = genes.up.human2$ortholog_ensg[genes.up.human2$ortholog_ensg != "N/A"]
genes.up.hs3 = genes.up.human3$ortholog_ensg[genes.up.human3$ortholog_ensg != "N/A"]
genes.up.hs3.PC1 = genes.up.human3.PC1$ortholog_ensg[genes.up.human3.PC1$ortholog_ensg != "N/A"]

# Down-regulated Model 1-3

genes.dn.human = gorth(
    sig.genes.dn_list,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

genes.dn.human2 = gorth(
  sig.genes.dn_list2,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15") 
genes.dn.human3 = gorth(
  sig.genes.dn_list3,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15") 
genes.dn.human3.PC1 = gorth(
  sig.genes.dn_list3_PC1,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

genes.dn.hs = genes.dn.human$ortholog_ensg[genes.dn.human$ortholog_ensg != "N/A"]
genes.dn.hs2 = genes.dn.human2$ortholog_ensg[genes.dn.human2$ortholog_ensg != "N/A"]
genes.dn.hs3 = genes.dn.human3$ortholog_ensg[genes.dn.human3$ortholog_ensg != "N/A"]
genes.dn.hs3.PC1 = genes.dn.human3.PC1$ortholog_ensg[genes.dn.human3.PC1$ortholog_ensg != "N/A"]

# Background
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15") 
bg.genes.human = gorth(
    bg.genes_list,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

bg.genes.hs = bg.genes.human$ortholog_ensg[bg.genes.human$ortholog_ensg != "N/A"]

# --- Test for significant overlap between GMTs and up- or down-regulated genes ------
# --- using hypergeometric test

gmt.files = c(
    "gene_sets/CTD_gene_disease_malaria.gmt",
    "gene_sets/DISEASE_textmining_malaria.gmt",
    "gene_sets/GAD_gene_disease_malaria.gmt",
    "gene_sets/Ebel_interacting_proteins.PPIPs.gmt",
    "gene_sets/Ebel_interacting_proteins.PPIPs-human.gmt",
    "gene_sets/HALLMARK_HEME_METABOLISM.gmt",
    "gene_sets/REACTOME_ERYTHROCYTES_TAKE_UP_OXYGEN_AND_RELEASE_CARBON_DIOXIDE.gmt",
    "gene_sets/REACTOME_ERYTHROCYTES_TAKE_UP_CARBON_DIOXIDE_AND_RELEASE_OXYGEN.gmt",
    "gene_sets/STEINER_ERYTHROCYTE_MEMBRANE_GENES.gmt"
)

hyper.p.up = do.call(rbind, lapply(gmt.files, function(gmt.file) {
    
    gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
    gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
    
    sig.in.gmt = sum(genes.up.hs %in% gmt)
    sig        = length(genes.up.hs)
    bg         = length(bg.genes.hs)
    gmt        = sum(gmt %in% bg.genes.hs)
    
    return(c(
        gmt.name,
        phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
               lower.tail = FALSE, log.p = FALSE)
    ))
}))

hyper.p.up2 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.up.hs2 %in% gmt)
  sig        = length(genes.up.hs2)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))

hyper.p.up3 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.up.hs3 %in% gmt)
  sig        = length(genes.up.hs3)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))

hyper.p.up3.PC1 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.up.hs3.PC1 %in% gmt)
  sig        = length(genes.up.hs3.PC1)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))


hyper.p.dn = do.call(rbind, lapply(gmt.files, function(gmt.file) {
    
    gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
    gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
    
    sig.in.gmt = sum(genes.dn.hs %in% gmt)
    sig        = length(genes.dn.hs)
    bg         = length(bg.genes.hs)
    gmt        = sum(gmt %in% bg.genes.hs)
    
    return(c(
        gmt.name,
        phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
               lower.tail = FALSE, log.p = FALSE)
    ))
}))

hyper.p.dn2 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.dn.hs2 %in% gmt)
  sig        = length(genes.dn.hs2)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))

hyper.p.dn3 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.dn.hs3 %in% gmt)
  sig        = length(genes.dn.hs3)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))

hyper.p.dn3.PC1 = do.call(rbind, lapply(gmt.files, function(gmt.file) {
  
  gmt.name = read.table(gmt.file, sep="\t", header=FALSE)[1]
  gmt = read.table(gmt.file, sep="\t", header=FALSE)[-c(1:2)]
  
  sig.in.gmt = sum(genes.dn.hs3.PC1 %in% gmt)
  sig        = length(genes.dn.hs3.PC1)
  bg         = length(bg.genes.hs)
  gmt        = sum(gmt %in% bg.genes.hs)
  
  return(c(
    gmt.name,
    phyper(sig.in.gmt - 1, sig, bg   - sig, gmt,
           lower.tail = FALSE, log.p = FALSE)
  ))
}))

# Write tables of a priori gene set enrichment results

write.table(hyper.p.up,
            file=paste0("Results_gene/malaria_GMT_overrep.upreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up2,
            file=paste0("Results_gene/malaria_GMT_overrep.upreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up3,
            file=paste0("Results_gene/malaria_GMT_overrep.upreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up3.PC1,
            file=paste0("Results_gene/malaria_GMT_overrep.upreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)


write.table(hyper.p.dn,
            file=paste0("Results_gene/malaria_GMT_overrep.downreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn2,
            file=paste0("Results_gene/malaria_GMT_overrep.downreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn3,
            file=paste0("Results_gene/malaria_GMT_overrep.downreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn3.PC1,
            file=paste0("Results_gene/malaria_GMT_overrep.downreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)


# --- Test whether parasitemia-linked genes are enriched for selected genes-------

# Van der Lee et al. 2017; doi: 10.1093/nar/gkx704, Table S4 for Models 1-3
sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel = sum(genes.up.hs %in% sel)
sig        = length(genes.up.hs)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(up.reg.sel.test = phyper(sig.in.sel - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel2 = sum(genes.up.hs2 %in% sel)
sig        = length(genes.up.hs2)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(up.reg.sel.test2 = phyper(sig.in.sel2 - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel3 = sum(genes.up.hs3 %in% sel)
sig        = length(genes.up.hs3)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(up.reg.sel.test3 = phyper(sig.in.sel3 - 1, sig, bg   - sig, sel,
                           lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel3.PC1 = sum(genes.up.hs3.PC1 %in% sel)
sig        = length(genes.up.hs3.PC1)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(up.reg.sel.test3.PC1 = phyper(sig.in.sel3.PC1 - 1, sig, bg   - sig, sel,
                           lower.tail = FALSE, log.p = FALSE))


#Down
sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel.dn = sum(genes.dn.hs %in% sel)
sig        = length(genes.dn.hs)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(dn.reg.sel.test = phyper(sig.in.sel.dn - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel.dn2 = sum(genes.dn.hs2 %in% sel)
sig        = length(genes.dn.hs2)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(dn.reg.sel.test2 = phyper(sig.in.sel.dn2 - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel.dn3 = sum(genes.dn.hs3 %in% sel)
sig        = length(genes.dn.hs3)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(dn.reg.sel.test3 = phyper(sig.in.sel.dn3 - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel.dn3.PC1 = sum(genes.dn.hs3.PC1 %in% sel)
sig        = length(genes.dn.hs3.PC1)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(dn.reg.sel.test3.PC1 = phyper(sig.in.sel.dn3.PC1 - 1, sig, bg   - sig, sel,
                           lower.tail = FALSE, log.p = FALSE))

# Write tables of selected gene set enrichment results for Models 1-3

write.table(up.reg.sel.test,
            file=paste0("Results_gene/dNdS_selected_overrep.upreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test2,
            file=paste0("Results_gene/dNdS_selected_overrep.upreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test3,
            file=paste0("Results_gene/dNdS_selected_overrep.upreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test3.PC1,
            file=paste0("Results_gene/dNdS_selected_overrep.upreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(dn.reg.sel.test,
            file=paste0("Results_gene/dNdS_selected_overrep.downreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test2,
            file=paste0("Results_gene/dNdS_selected_overrep.downreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test3,
            file=paste0("Results_gene/dNdS_selected_overrep.downreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test3.PC1,
            file=paste0("Results_gene/dNdS_selected_overrep.downreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Make plots for overrep analysis ----------
# Up-regulated (Nothing significant in down-regulated genes except for PC1 of Model 3)
#Model1
gmt.up1=read.table("Results_gene/malaria_GMT_overrep.upreg_mod1.txt")
colnames(gmt.up1)=c("GMT", "P_Val")
gmt.up1$log_p_val=-log10(gmt.up1$P_Val)
gmt.up1 <- gmt.up1[order(gmt.up1$log_p_val), ]
gmt.up1$GMT <- factor(gmt.up1$GMT, levels = gmt.up1$GMT[order(gmt.up1$log_p_val)])

ggplot(gmt.up1, aes(GMT, log_p_val, fill=log_p_val)) +
    scale_fill_gradient2(name="Log10(P-Val)", low="black", mid="grey",
                         high="firebrick") +
    geom_bar(stat="identity") +
    ylab("-Log10(P-Value)") +
    coord_flip() +
    theme_classic() +
    ggtitle("Model 1") +
    geom_hline(yintercept = 1.3, linetype="dotted", 
               color = "red")

#Model2
gmt.up2=read.table("Results_gene/malaria_GMT_overrep.upreg_mod2.txt")
colnames(gmt.up2)=c("GMT", "P_Val")
gmt.up2$log_p_val=-log10(gmt.up2$P_Val)
gmt.up2 <- gmt.up2[order(gmt.up2$log_p_val), ]
gmt.up2$GMT <- factor(gmt.up2$GMT, levels = gmt.up2$GMT[order(gmt.up2$log_p_val)])

ggplot(gmt.up2, aes(GMT, log_p_val, fill=log_p_val)) +
  scale_fill_gradient2(name="Log10(P-Val)", low="black", mid="grey",
                       high="firebrick") +
  geom_bar(stat="identity") +
  ylab("-Log10(P-Value)") +
  coord_flip() +
  theme_classic() +
  ggtitle("Model 2") +
  geom_hline(yintercept = 1.3, linetype="dotted", 
             color = "red")

#Model3
gmt.up3=read.table("Results_gene/malaria_GMT_overrep.upreg_mod3.txt")
colnames(gmt.up3)=c("GMT", "P_Val")
gmt.up3$log_p_val=-log10(gmt.up3$P_Val)
gmt.up3 <- gmt.up3[order(gmt.up3$log_p_val), ]
gmt.up3$GMT <- factor(gmt.up3$GMT, levels = gmt.up3$GMT[order(gmt.up3$log_p_val)])

ggplot(gmt.up3, aes(GMT, log_p_val, fill=log_p_val)) +
  scale_fill_gradient2(name="Log10(P-Val)", low="black", mid="grey",
                       high="firebrick") +
  geom_bar(stat="identity") +
  ylab("-Log10(P-Value)") +
  coord_flip() +
  theme_classic() +
  ggtitle("Model 3") +
  geom_hline(yintercept = 1.3, linetype="dotted", 
             color = "red")


gmt.up3.PC1=read.table("Results_gene/malaria_GMT_overrep.upreg_mod3_PC1.txt")
colnames(gmt.up3.PC1)=c("GMT", "P_Val")
gmt.up3.PC1$log_p_val=-log10(gmt.up3.PC1$P_Val)
gmt.up3.PC1 <- gmt.up3.PC1[order(gmt.up3.PC1$log_p_val), ]
gmt.up3.PC1$GMT <- factor(gmt.up3.PC1$GMT, levels = gmt.up3.PC1$GMT[order(gmt.up3.PC1$log_p_val)])

ggplot(gmt.up3.PC1, aes(GMT, log_p_val, fill=log_p_val)) +
  scale_fill_gradient2(name="Log10(P-Val)", low="black", mid="grey",
                       high="firebrick") +
  geom_bar(stat="identity") +
  ylab("-Log10(P-Value)") +
  coord_flip() +
  theme_classic() +
  ggtitle("Model 3 (coef=PC1 Up-regulated genes)") +
  geom_hline(yintercept = 1.3, linetype="dotted", 
             color = "red")

gmt.dn3.PC1=read.table("Results_gene/malaria_GMT_overrep.downreg_mod3_PC1.txt")
colnames(gmt.dn3.PC1)=c("GMT", "P_Val")
gmt.dn3.PC1$log_p_val=-log10(gmt.dn3.PC1$P_Val)
gmt.dn3.PC1 <- gmt.dn3.PC1[order(gmt.dn3.PC1$log_p_val), ]
gmt.dn3.PC1$GMT <- factor(gmt.dn3.PC1$GMT, levels = gmt.dn3.PC1$GMT[order(gmt.dn3.PC1$log_p_val)])

ggplot(gmt.dn3.PC1, aes(GMT, log_p_val, fill=log_p_val)) +
  scale_fill_gradient2(name="Log10(P-Val)", low="black", mid="grey",
                       high="firebrick") +
  geom_bar(stat="identity") +
  ylab("-Log10(P-Value)") +
  coord_flip() +
  theme_classic() +
  ggtitle("Model 3 (coef=PC1 Down-regulated genes)") +
  geom_hline(yintercept = 1.3, linetype="dotted", 
             color = "red")


### --> Overrep analysis using piN/piS

# Look at distribution of piN/piS across all genes

piN_piS = read.csv("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Second_analysis/snp_results/snp_genie/piN_piS.csv")

piN_piS$piN_piS = as.numeric(piN_piS$piN_piS)
piN_piS = na.omit(piN_piS)

hist(piN_piS$piN_piS, xlim=c(0,1.5))

# EMMREML overrep (FDR=0.2)

sel = subset(piN_piS, piN_piS > 1, select = c(product, piN_piS))

sig.in.sel = sum(sig.genes.up_list %in% sel)
sig        = length(sig.genes.up_list)
bg         = length(bg.genes_list)
sel        = sum(sel %in% bg.genes_list)

(up.reg.sel.test = phyper(sig.in.sel - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sig.in.sel = sum(sig.genes.dn_list %in% sel)
sig        = length(sig.genes.dn_list)
bg         = length(bg.genes_list)
sel        = sum(sel %in% bg.genes_list)

(dn.reg.sel.test = phyper(sig.in.sel - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

# Write results tables

write.table(up.reg.sel.test,
            file=paste0("piNpiS_selected_overrep.upreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test,
            file=paste0("piNpiS_selected_overrep.downreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)


#### Investigate gene vs. exon nonsense ####
library(ggplot2)
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR")

proxy=read.csv("parasitemia_proxy.csv")

# Correlation between "exon" and "gene" proxies
ggplot(proxy, aes(Parasitemia_proxy_edgeR, Parasitemia_proxy_gene)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  xlab("Proxy by 'exon'") +
  ylab("Proxy by 'gene'")

# Correlation between "exon" and "gene" read counts - Colobus
ggplot(proxy, aes(Colobus_Reads_Mapped, Colobus_Reads_gene)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  xlab("Read count by 'exon'") +
  ylab("Read count by 'gene'") +
  labs(title="Colobus")

# Correlation between "exon" and "gene read counts - Hepatocystis
ggplot(proxy, aes(Hepatocystis_Reads_Mapped, Hepatocystis_Reads_gene)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme_minimal() +
  xlab("Read count by 'exon'") +
  ylab("Read count by 'gene'") +
  labs(title="Hepatocystis")

# Look at total read counts per individual

proxy2 = data.frame(proxy$Sample_name)
colnames(proxy2)[1] = "Individuals"
proxy2$Total_reads = proxy$Total_Reads_Hepatocystis_mapping
proxy2$Mapping_method = "exon"

proxy3 = data.frame(proxy$Sample_name)
colnames(proxy3)[1] = "Individuals"
proxy3$Total_reads = proxy$Total_Reads_gene
proxy3$Mapping_method = "gene"

proxy4 = rbind(proxy2, proxy3)

ggplot(proxy4, aes(fill=Mapping_method, y=Total_reads, x=Individuals)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Plot genes of interest
plot.gene.expr = function(gene.name) {
  
  cpm.gene_of_interest = cpm[row.names(cpm) == gene.name,]
  ct = data.frame(counts      = cpm.gene_of_interest,
                  parasitemia = parasitemia$parasitemia.proxy)
  
  p = ggplot(ct, aes(parasitemia, counts)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    xlab("Inferred parasitemia proxy") +
    ylab("Gene Expression\n(Normalized count per million reads)") +
    theme_bw()
  ggsave(p, file=paste0("Results_gene/Figures/cpm_by_parasitemia.", gene.name, ".", suffix= "pdf"))
}

#up
lapply(c("PAPPA2", "ABCB11", "FBXW11", "TMEM132D"), plot.gene.expr)

# Individual genes
cpm.gene_of_interest.PAPPA2 = data.frame(cpm[row.names(cpm) == "PAPPA2",])
colnames(cpm.gene_of_interest.PAPPA2)[1] = "CPM"
cpm.gene_of_interest.PAPPA2$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.PAPPA2$gene = "PAPPA2"

cpm.gene_of_interest.ABCB11 = data.frame(cpm[row.names(cpm) == "ABCB11",])
colnames(cpm.gene_of_interest.ABCB11)[1] = "CPM"
cpm.gene_of_interest.ABCB11$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.ABCB11$gene = "ABCB11"

cpm.gene_of_interest.FBXW11 = data.frame(cpm[row.names(cpm) == "FBXW11",])
colnames(cpm.gene_of_interest.FBXW11)[1] = "CPM"
cpm.gene_of_interest.FBXW11$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.FBXW11$gene = "FBXW11"

cpm.gene_of_interest.TMEM132D = data.frame(cpm[row.names(cpm) == "TMEM132D",])
colnames(cpm.gene_of_interest.TMEM132D)[1] = "CPM"
cpm.gene_of_interest.TMEM132D$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.TMEM132D$gene = "TMEM132D"

cpm.gene_of_interest = rbind(cpm.gene_of_interest.PAPPA2, cpm.gene_of_interest.ABCB11, cpm.gene_of_interest.FBXW11, cpm.gene_of_interest.TMEM132D)

ggplot(cpm.gene_of_interest, aes(parasitemia, CPM)) +
  geom_point(color='firebrick') +
  geom_smooth(method="lm", se=FALSE, color="black") +
  facet_wrap(~ gene, scales="free") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw()

#Down
cpm.gene_of_interest.TBC1D9B = data.frame(cpm[row.names(cpm) == "TBC1D9B",])
colnames(cpm.gene_of_interest.TBC1D9B)[1] = "CPM"
cpm.gene_of_interest.TBC1D9B$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.TBC1D9B$gene = "TBC1D9B"

cpm.gene_of_interest.HRH2 = data.frame(cpm[row.names(cpm) == "HRH2",])
colnames(cpm.gene_of_interest.HRH2)[1] = "CPM"
cpm.gene_of_interest.HRH2$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HRH2$gene = "HRH2"

cpm.gene_of_interest.NAGA = data.frame(cpm[row.names(cpm) == "NAGA",])
colnames(cpm.gene_of_interest.NAGA)[1] = "CPM"
cpm.gene_of_interest.NAGA$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.NAGA$gene = "NAGA"

cpm.gene_of_interest_down = rbind(cpm.gene_of_interest.TBC1D9B, cpm.gene_of_interest.HRH2, cpm.gene_of_interest.NAGA)

ggplot(cpm.gene_of_interest_down, aes(parasitemia, CPM)) +
  geom_point(color='firebrick') +
  geom_smooth(method="lm", se=FALSE, color="black") +
  facet_wrap(~ gene, scales="free") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw()

#ACKR1
cpm.gene_of_interest.ACKR1 = data.frame(cpm[row.names(cpm) == "ACKR1",])
colnames(cpm.gene_of_interest.ACKR1)[1] = "CPM"
cpm.gene_of_interest.ACKR1$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.ACKR1$gene = "ACKR1"
cpm.gene_of_interest.ACKR1$Individuals = rownames(cpm.gene_of_interest.ACKR1)
cpm.gene_of_interest.ACKR1$Individuals=str_sub(cpm.gene_of_interest.ACKR1$Individuals,1,5)

#With outlier
ggplot(cpm.gene_of_interest.ACKR1, aes(parasitemia, CPM)) +
  geom_point(color='firebrick') +
  geom_smooth(method="lm", se=FALSE, color="black") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw() +
  geom_text(data=subset(cpm.gene_of_interest.ACKR1, CPM > 15),
            aes(parasitemia,CPM,label=Individuals), vjust=2)

#withoutoutlier
ggplot(cpm.gene_of_interest.ACKR1, aes(parasitemia, CPM)) +
  geom_point(color='firebrick') +
  geom_smooth(method="lm", se=FALSE, color="black") +
  ylim(0, 7.5) +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  ggtitle("ACKR1") +
  theme_bw() 

# Make plots comparing Model 1 and Model 3
Mod1_Mod3_mod1 = data.frame(rownames(all.sig.genes))
colnames(Mod1_Mod3_mod1) <- c('GeneID')
Mod1_Mod3_mod1$logFC_mod1 = all.sig.genes$logFC
Mod1_Mod3_mod1$FDR_mod1 = all.sig.genes$FDR

Mod1_Mod3_mod3 = data.frame(rownames(all.sig.genes3))
colnames(Mod1_Mod3_mod3) <- c('GeneID')
Mod1_Mod3_mod3$logFC_mod3 = all.sig.genes3$logFC
Mod1_Mod3_mod3$FDR_mod3 = all.sig.genes3$FDR

Mod1_Mod3 = merge(Mod1_Mod3_mod1, Mod1_Mod3_mod3, by='GeneID', all=TRUE)

Mod1_Mod3$Expression = "blank"

Mod1_Mod3$Expression[is.na(Mod1_Mod3$logFC_mod1)==TRUE] <- "Model 3"
Mod1_Mod3$Expression[is.na(Mod1_Mod3$logFC_mod3)==TRUE] <- "Model 1"

Bobowik_falciparum$Study[is.na(Bobowik_falciparum$Bobowik_falciparumIndoLogFC)==FALSE & is.na(Bobowik_falciparum$falciparumLogFC)==FALSE] <- "Both"
Bobowik_falciparum$Study=as.factor(Bobowik_falciparum$Study)


if (Mod1_Mod3_mod1$logFC_mod1 > 0
