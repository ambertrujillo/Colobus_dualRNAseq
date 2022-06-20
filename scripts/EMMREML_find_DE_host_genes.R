#################################
####### Voom --> EMMREML  #######
#################################

######### --> NORMALIZED EXPRESSION MATRIX
# save expression matrix (log2CPM normalized) from voom
#library("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)
load("EMMREML_results/colobus.EMMREML.fc.Rdata")

fc = colobus.EMMREML.fc
head(fc$counts)
dim(fc$counts)
# 45274 29

# Make DGEList of counts
d0 <- DGEList(fc$counts)

# Calculate normalization factors for voom 
d0 <- calcNormFactors(d0)
dim(d0)
#[1] 45274    29

# Apply a cutoff
cutoff <- 1 # 1 cpm 
drop <- which(apply(cpm(d0), 1, max) < cutoff) # in at least 1 sample
d <- d0[-drop,] 
dim(d) # number of genes left
#[1] 16221    29

################## Voom transformation and calculation of variance weights ######################
parasitemia=read.csv("EMMREMLparasitemia_proxy.csv")
parasitemia$parasitemia.proxy =
  parasitemia$Hepatocystis_reads_gene /
  parasitemia$Colobus_reads_gene

# Create Model Matrix
mm <- model.matrix(~parasitemia.proxy + Sex, data = parasitemia)

# Voom
y <- voom(d, mm, plot = T)

# Save expression matrix from voom
colnames(y$E)=gsub(".colobus.bam", "_001", colnames(y$E))
e.keep = y$E

########## --> RELATEDNESS MATRIX (k matrix)

# load metadata (parasitemia table), check if in correct order
table(gsub("_001", "", parasitemia$Sample_name==colnames(e.keep)))

# Put together k and z matrices
#animals = sort(unique(meta$Individual))

# Calculate k matrix (pairwise relatedness) and create an identity matrix for the second model

library(cultevo)
k1 <- read.csv("Results_gene/colobus.plink.relatedness.csv")
rownames(k1) = k1$X
k1 <- k1[,-1]
k1 <- as.matrix(k1, ncol=29, nrow=29, byrow=TRUE, bycol=TRUE) 
diag(k1) = 1
class(k1)

# Identity matrix for second (no relatedness model)
k2 = matrix(0, 29, 29)
rownames(k2) = rownames(k1)
colnames(k2) = colnames(k1)
diag(k2) = 1

########## --> MATRIX LINKING SAMPLES TO ANIMALS (z matrix) 

z1 = matrix(0,nrow=nrow(k1),ncol=ncol(k1))
rownames(z1) = rownames(k1)
colnames(z1) = colnames(k1)
diag(z1) = 1

z2 = matrix(0,nrow=nrow(k2),ncol=ncol(k2))
rownames(z2) = rownames(k2)
colnames(z2) = colnames(k2)
diag(z2) = 1

########## --> MODEL DESIGN

library(parallel)
library(doParallel)

# Initialize output
out = list()

n.cores = detectCores() - 2

# Design model covariates

model.covariates = c('parasitemia.proxy', 'Sex')
design = model.matrix(as.formula(paste('~',paste(model.covariates,collapse=' + '))), data=parasitemia)

# set up cluster for parallel processing model 1_relate

clus = parallel::makeCluster(n.cores, setup_timeout = 0.5)
registerDoParallel(cores=n.cores)  
clusterExport(clus,varlist=c('e.keep','k1','z1','design'),envir=environment())

# run model

out1 = t(parApply(clus, e.keep,1,function(y) {
  require(EMMREML)
  
  emma=emmreml(y = y,X = design,Z = z1,K = k1,varbetahat = T,varuhat = T,PEVuhat = T,test = T)
  
  p = emma$pvalbeta
  varb = emma$varbetahat
  b = emma$betahat
  
  c(b,varb,p[,"none"])
}))

stopCluster(clus)

colnames(out1)[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
colnames(out1)[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
colnames(out1)[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

out1 = data.frame(out1)

# Benjamini & Hochberg correction
out1$parasitema.fdr = p.adjust(out1$pval.parasitemia.proxy, method = 'BH')
out1$sex.fdr = p.adjust(out1$pval.Sexmale, method = 'BH')

hist(out1$parasitema.fdr, breaks = 500, xlim = c(0,0.1))
sum(out1$parasitema.fdr < 0.2)
#3556 genes are related to parasitemia

# Look at genes significant for/impacted by sex
hist(out1$sex.fdr, breaks = 500, xlim = c(0,0.1))
sum(out1$sex.fdr < 0.2)
#372 genes are related to sex
sig.genes.up.0.2.sex.Mod1 = subset(out1, beta.Sexmale > 0 & sex.fdr < 0.2, select=c(beta.Sexmale, sex.fdr), quote=FALSE)
#194
sig.genes.dn.0.2.sex.Mod1 = subset(out1, beta.Sexmale < 0 & sex.fdr < 0.2, select=c(beta.Sexmale, sex.fdr))
#178

# Get up-regulated, down-regulated, and affected significant genes
## up
sig.genes.up.0.2.Mod1 = subset(out1, beta.parasitemia.proxy > 0 & parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr), quote=FALSE)
# 1297

#for overrep analysis
sig.up.0.2_list.Mod1=rownames(sig.genes.up.0.2.Mod1)

## down
sig.genes.dn.0.2.Mod1 = subset(out1, beta.parasitemia.proxy < 0 & parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr))
# 2259

#for overrep analysis
sig.dn.0.2_list.Mod1=rownames(sig.genes.dn.0.2.Mod1)

bg_list_Mod1=rownames(out1) # background list

## all
sig.genes.all.0.2.Mod1 = subset(out1, parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr))
# 3556

# Write tables of results
#up
write.table(sig.genes.up.0.2.Mod1, file="EMMREML_results/emmreml_sig_up_02_Mod1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
#down
write.table(sig.genes.dn.0.2.Mod1, file="EMMREML_results/emmreml_sig_dn_02_Mod1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
#all
write.table(sig.genes.all.0.2.Mod1, file="EMMREML_results/emmreml_sig_all_02_Mod1.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

# set up cluster for parallel processing model 2_relate

clus = parallel::makeCluster(n.cores, setup_timeout = 0.5)
registerDoParallel(cores=n.cores)  
clusterExport(clus,varlist=c('e.keep','k2','z2','design'),envir=environment())

# run model

out2 = t(parApply(clus, e.keep,1,function(y) {
  require(EMMREML)
  
  emma=emmreml(y = y,X = design,Z = z2,K = k2,varbetahat = T,varuhat = T,PEVuhat = T,test = T)
  
  p = emma$pvalbeta
  varb = emma$varbetahat
  b = emma$betahat
  
  c(b,varb,p[,"none"])
}))

stopCluster(clus)

colnames(out2)[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
colnames(out2)[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
colnames(out2)[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

out2 = data.frame(out2)

# Benjamini & Hochberg correction
out2$parasitema.fdr = p.adjust(out2$pval.parasitemia.proxy, method = 'BH')
out2$sex.fdr = p.adjust(out2$pval.Sexmale, method = 'BH')

hist(out2$parasitema.fdr, breaks = 500, xlim = c(0,0.1))
sum(out2$parasitema.fdr < 0.2)
#2940 genes are related to parasitemia

# Look at genes significant for/impacted by sex
hist(out2$sex.fdr, breaks = 500, xlim = c(0,0.1))
sum(out2$sex.fdr < 0.2)
#36 genes are related to sex
sig.genes.up.0.2.sex.Mod2 = subset(out2, beta.Sexmale > 0 & sex.fdr < 0.2, select=c(beta.Sexmale, sex.fdr), quote=FALSE)
#23
sig.genes.dn.0.2.sex.Mod2 = subset(out2, beta.Sexmale < 0 & sex.fdr < 0.2, select=c(beta.Sexmale, sex.fdr))
#13

# Get up-regulated, down-regulated, and affected significant genes
## up
sig.genes.up.0.2.Mod2 = subset(out2, beta.parasitemia.proxy > 0 & parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr), quote=FALSE)
# 1061

#for overrep analysis
sig.up.0.2_list.Mod2=rownames(sig.genes.up.0.2.Mod2)

## down
sig.genes.dn.0.2.Mod2 = subset(out2, beta.parasitemia.proxy < 0 & parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr))
# 1879

#for overrep analysis
sig.dn.0.2_list.Mod2=rownames(sig.genes.dn.0.2.Mod2)

bg_list_Mod2=rownames(out2) # background list

## all
sig.genes.all.0.2.Mod2 = subset(out2, parasitema.fdr < 0.2, select=c(beta.parasitemia.proxy, parasitema.fdr))
# 2940

# Write tables of results
#up
write.table(sig.genes.up.0.2.Mod2, file="EMMREML_results/emmreml_sig_up_02_Mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
#down
write.table(sig.genes.dn.0.2.Mod2, file="EMMREML_results/emmreml_sig_dn_02_Mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
#all
write.table(sig.genes.all.0.2.Mod2, file="EMMREML_results/emmreml_sig_all_02_Mod2.txt",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

########## --> See if beta (voom-EMMREML) and LogFC (edgeR) are correlated

Logfc = read.csv("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR/Results_gene/colobus_bg_genes.csv")
beta1 = out1

beta1$GeneID = rownames(out1)

# Merge the two data frames to match Gene IDs
DE_cor = merge(Logfc, beta1, by="GeneID", all = T)
DE_cor = na.omit(DE_cor) # Take out rows with NAs

library(ggplot2)
plot(DE_cor$logFC, DE_cor$beta.parasitemia.proxy)

d = lm(beta.parasitemia.proxy ~ logFC, data=DE_cor)
summary(d)
# pval = < 0.05, R2 = 0.2656 

ggplot(DE_cor, aes(logFC, beta.parasitemia.proxy)) +
  geom_point() +
  xlim(-15, 20) +
  ylim(-20, 20) +
  xlab("LogFC") +
  ylab("Effect Size (Beta; Inferred Parasitemia)") +
  geom_abline(intercept = 0.0217, slope = 0.9002) +
  theme_classic()

cor(DE_cor$logFC, DE_cor$beta.parasitemia.proxy)

# See if Betas correlate between Model 1 and Model 2
cor_betas1 = data.frame(rownames(out1))
colnames(cor_betas1)[1]="GeneID"
cor_betas1$beta.parasitemia.proxy1 = out1$beta.parasitemia.proxy
cor_betas1$beta.Sexmale1 = out1$beta.Sexmale
cor_betas1$parasitemia.fdr1 = out1$parasitema.fdr
cor_betas1$sex.fdr1 = out1$sex.fdr

cor_betas2 = data.frame(rownames(out2))
colnames(cor_betas2)[1]="GeneID"
cor_betas2$beta.parasitemia.proxy2 = out2$beta.parasitemia.proxy
cor_betas2$beta.Sexmale2 = out2$beta.Sexmale
cor_betas2$parasitemia.fdr2 = out2$parasitema.fdr
cor_betas2$sex.fdr2 = out2$sex.fdr

cor_betas = merge(cor_betas1, cor_betas2, by="GeneID", all = T)
cor_betas = na.omit(cor_betas) # Take out rows with NAs

#betas - parasitemia proxy
ggplot(cor_betas, aes(beta.parasitemia.proxy1, beta.parasitemia.proxy2)) +
  geom_point() +
  theme_bw() +
  ggtitle("Betas - Parasitemia Proxy") +
  xlab("Model 1 - Relatedness") +
  ylab("Model 2 - No Relatedness")

betas_parasitemia = lm(beta.parasitemia.proxy1 ~ beta.parasitemia.proxy2, data=cor_betas)
summary(betas_parasitemia)

#betas - sex
ggplot(cor_betas, aes(beta.Sexmale1, beta.Sexmale2)) +
  geom_point() +
  theme_bw() +
  ggtitle("Betas - Sex") +
  xlab("Model 1 - Relatedness") +
  ylab("Model 2 - No Relatedness")
  

betas_sex = lm(beta.Sexmale1 ~ beta.Sexmale2, data=cor_betas)
summary(betas_sex)

#FDR - parasitemia proxy
ggplot(cor_betas, aes(parasitemia.fdr1, parasitemia.fdr2)) +
  geom_point() +
  theme_bw() +
  ggtitle("FDR - Parasitemia Proxy") +
  xlab("Model 1 - Relatedness") +
  ylab("Model 2 - No Relatedness")

#FDR - sex
ggplot(cor_betas, aes(sex.fdr1, sex.fdr2)) +
  geom_point() + 
  theme_bw() +
  ggtitle("FDR - Sex") +
  xlab("Model 1 - Relatedness") +
  ylab("Model 2 - No Relatedness")

# Write a table with specific genes in mind

rownames(DE_cor) = DE_cor$GeneID
export <- DE_cor[c("ACKR1", "NAGA", "UBE2K", "TMEM167A", "LSM14A", "UBE2B", "TTLL12", "TBC1D9B", "AGFG1", "APOBEC2", "ZNF682"),]
genes_int = subset(export, select = c(beta.parasitemia.proxy, parasitema.fdr, logFC, FDR))
#PP2D1 was not in both data sets

# Change the number of decimal places
is.num <- sapply(genes_int, is.numeric)
genes_int[is.num] <- lapply(genes_int[is.num], round, 8)

write.table(genes_int, file="results/emmreml_edgeR_genes.csv",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

genes_int$GeneID = rownames(genes_int)

cor(genes_int$beta.parasitemia.proxy, genes_int$logFC)
lm(logFC ~ beta.parasitemia.proxy, data=genes_int)
ggplot(genes_int, aes(beta.parasitemia.proxy, logFC, color=GeneID)) +
  geom_point() +
  geom_abline(intercept = 0.04406, slope = 1.17309) +
  theme_classic()

# Only include columns of interest that are significant
genes_int_sig = subset(genes_int, parasitema.fdr < 0.1 & FDR < 0.1)

genes_int_sig$GeneID = rownames(genes_int_sig)

lm(logFC ~ beta.parasitemia.proxy, data=genes_int_sig)
ggplot(genes_int_sig, aes(beta.parasitemia.proxy, logFC, color=GeneID)) +
  geom_point() +
  geom_abline(intercept = -0.05493, slope = 1.03125) +
  theme_classic()

### --> Look at blood genes that were significant in original analysis

export2 <- DE_cor[c("RHAG", "SPTA1", "KLF1", "SLC4A1", "ABCB6", "SLC2A1", "STEAP3", "JAK2"),]
genes_int_blood = subset(export2, select = c(bvar.parasitemia.proxy,beta.parasitemia.proxy, parasitema.fdr, logFC, FDR))

# plot them
is.num <- sapply(genes_int_blood, is.numeric)
genes_int_blood[is.num] <- lapply(genes_int_blood[is.num], round, 8)

write.table(genes_int_blood, file="results/emmreml_edgeR_blood_genes.csv",
            sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

genes_int_blood$GeneID = rownames(genes_int_blood)

cor(genes_int_blood$beta.parasitemia.proxy, genes_int_blood$logFC)
lm(logFC ~ beta.parasitemia.proxy, data=genes_int_blood)
ggplot(genes_int_blood, aes(beta.parasitemia.proxy, logFC, color=GeneID)) +
  geom_point() +
  geom_abline(intercept = 0.4027, slope = 1.4566) +
  theme_classic()

# --- Do functional enrichment test on GMT files - Convert genes to human
library(gprofiler2)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15") 

# Up-regulated FDR = 0.2

genes.up.human = gorth(
  sig.up.0.2_list,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

genes.up.hs = genes.up.human$ortholog_ensg[genes.up.human$ortholog_ensg != "N/A"]

# Down-regulated

genes.dn.human = gorth(
  sig.dn.0.2_list,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

genes.dn.hs = genes.dn.human$ortholog_ensg[genes.dn.human$ortholog_ensg != "N/A"]

# Background

bg.genes.human = gorth(
  bg_list,
  source_organism = "ptephrosceles",
  target_organism = "hsapiens",
  numeric_ns = "ENTREZGENE_ACC",
  mthreshold = Inf,
  filter_na = TRUE
)

bg.genes.hs = bg.genes.human$ortholog_ensg[bg.genes.human$ortholog_ensg != "N/A"]

# --- Test for significant overlap between GMTs and up- or down-regulated genes
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

# Write tables of a priori gene set enrichment results

write.table(hyper.p.up,
            file=paste0("Results_gene/malaria_GMT_overrep.upreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn,
            file=paste0("Results_gene/malaria_GMT_overrep.downreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- Test whether parasitemia-linked genes are enriched for selected genes

# Van der Lee et al. 2017; doi: 10.1093/nar/gkx704, Table S4
sel = read.table("gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
                 sep="\t", header=TRUE)$Ensembl.Gene.ID

sig.in.sel = sum(genes.up.hs %in% sel)
sig        = length(genes.up.hs)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(up.reg.sel.test = phyper(sig.in.sel - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

sig.in.sel = sum(genes.dn.hs %in% sel)
sig        = length(genes.dn.hs)
bg         = length(bg.genes.hs)
sel        = sum(sel %in% bg.genes.hs)

(dn.reg.sel.test = phyper(sig.in.sel - 1, sig, bg   - sig, sel,
                          lower.tail = FALSE, log.p = FALSE))

# Write tables of selected gene set enrichment results

write.table(up.reg.sel.test,
            file=paste0("Results_gene/dNdS_selected_overrep.upreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test,
            file=paste0("Results_gene/dNdS_selected_overrep.downreg.0.2", suffix=".txt"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
