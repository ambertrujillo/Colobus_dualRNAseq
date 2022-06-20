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
PCA = read.table("edgeR_esults/PCA_results")
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

write.table(ref.list, file="edgeR_results/reference_list.colobus.txt",
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
qlf3 = glmQLFTest(fit3, coef = "parasitemia.proxy")
tt3 = topTags(qlf3, n=Inf, adjust.method = "BH", p.value = 1)

qlf3.PC1 = glmQLFTest(fit3, coef = "PCA_1")
tt3.PC1 = topTags(qlf3.PC1, n=Inf, adjust.method = "BH", p.value = 1)

qlf3.PC2 = glmQLFTest(fit3, coef = "PCA_2")
tt3.PC2 = topTags(qlf3.PC2, n=Inf, adjust.method = "BH", p.value = 1)

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

write.table(sig.genes.up, file="edgeR_results/colobus_upreg_sig_genes_mod1.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up2, file="edgeR_results/colobus_upreg_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up3, file="edgeR_results/colobus_upreg_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.up3.PC1, file="edgeR_results/colobus_upreg_sig_genes_mod3_PC1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)

write.table(sig.genes.dn, file="edgeR_results/colobus_dnreg_sig_genes_mod1.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn2, file="edgeR_results/colobus_dnreg_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn3, file="edgeR_results/colobus_dnreg_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(sig.genes.dn3.PC1, file="edgeR_results/colobus_dnreg_sig_genes_mod3_PC1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)

write.table(bg.genes, file="edgeR_results/colobus_bg_genes.txt",
    sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
    
write.table(all.sig.genes, file="edgeR_results/colobus_all_sig_genes_mod1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes2, file="edgeR_results/colobus_all_sig_genes_mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes3, file="edgeR_results/colobus_all_sig_genes_mod3.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE)
write.table(all.sig.genes3.PC1, file="edgeR_results/colobus_all_sig_genes_mod3_PC1.txt",
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
            file=paste0("edgeR_results/malaria_GMT_overrep.upreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up2,
            file=paste0("edgeR_results/malaria_GMT_overrep.upreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up3,
            file=paste0("edgeR_results/malaria_GMT_overrep.upreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.up3.PC1,
            file=paste0("edgeR_results/malaria_GMT_overrep.upreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)


write.table(hyper.p.dn,
            file=paste0("edgeR_results/malaria_GMT_overrep.downreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn2,
            file=paste0("edgeR_results/malaria_GMT_overrep.downreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn3,
            file=paste0("edgeR_results/malaria_GMT_overrep.downreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn3.PC1,
            file=paste0("edgeR_results/malaria_GMT_overrep.downreg_mod3_PC1", suffix=".txt"),
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
            file=paste0("edgeR_results/dNdS_selected_overrep.upreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test2,
            file=paste0("edgeR_results/dNdS_selected_overrep.upreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test3,
            file=paste0("edgeR_results/dNdS_selected_overrep.upreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(up.reg.sel.test3.PC1,
            file=paste0("edgeR_results/dNdS_selected_overrep.upreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(dn.reg.sel.test,
            file=paste0("edgeR_results/dNdS_selected_overrep.downreg_mod1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test2,
            file=paste0("edgeR_results/dNdS_selected_overrep.downreg_mod2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test3,
            file=paste0("edgeR_results/dNdS_selected_overrep.downreg_mod3", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test3.PC1,
            file=paste0("edgeR_results/dNdS_selected_overrep.downreg_mod3_PC1", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
