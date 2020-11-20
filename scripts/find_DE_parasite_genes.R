#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# To install BioConductor packages:
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install("edgeR")
#  BiocManager::install("GO.db")

library(limma)
library(edgeR)

library(ggplot2)
library(ggrepel)

library(plyr)
library(xtable)

library(gprofiler2)

library("GO.db")

# --- Load variables from simple DE analysis

load("after_DE_analysis.Aunin_med.Rdata")

# --- Load feature count data

aunin = data.frame(t(read.table("data/Aunin_etal_2020_read_counts.csv",
    header=TRUE, row.names=1)))

ind.info = read.table("data/ind_info.txt", header=TRUE)

aunin = merge(aunin, ind.info[,c(1,3)], by.x="row.names", by.y="Sample")

aunin.ind = aggregate(aunin[2:(ncol(aunin) - 1)],
    by=list(aunin$Sample_name), FUN=sum)

names(aunin.ind)[1] = "Sample"

# Prep data for use in edgeR
aunin.counts = data.frame(t(aunin.ind[,-c(1)]))
names(aunin.counts) = paste0(aunin.ind[,1], ".hepato.bam")
aunin.anno = data.frame(GeneID=names(aunin.ind)[-c(1)])

# Remove RC127 outlier
#aunin.counts = aunin.counts[,!grepl("RC127", names(aunin.counts))]       

# --- Load parasitemia info

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

# Remove RC127 outlier
#parasitemia = parasitemia[parasitemia$Sample_name != "RC127",]         

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts=aunin.counts, genes=aunin.anno)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge) > 10) >= 10

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

design = model.matrix(~ parasitemia.proxy, data=parasitemia)

# colnames(design)

disp = estimateDisp(fc.dge.norm, design, robust = TRUE)

# CPM
cpm.disp = cpm(disp)

fit = glmQLFit(disp, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = "parasitemia.proxy")
tt  = topTags(qlf, n=Inf, adjust.method = "BH", p.value = 1)

# Write list of significant up- and down-regulated genes
# Note based on raw, unadjusted p-value
sig.genes.up = tt$table$GeneID[tt$table$PValue < 0.01 & tt$table$logFC > 0]
length(sig.genes.up)

sig.genes.dn = tt$table$GeneID[tt$table$PValue < 0.01 & tt$table$logFC < 0]
length(sig.genes.dn)

# Write background list of genes
bg.genes = tt$table$GeneID
length(bg.genes)

# Write results to file
write.table(sig.genes.up, file="reports/edgeR_hepato_upreg_sig_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(sig.genes.dn, file="reports/edgeR_hepato_dnreg_sig_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(bg.genes, file="reports/edgeR_hepato_bg_genes.txt",
    sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# --- Plot genes of interest

plot.gene.expr = function(gene.name) {

    cpm.disp.gene_of_interest = cpm.disp[row.names(cpm.disp) == gene.name,]
    ct = data.frame(counts      = cpm.disp.gene_of_interest,
                    parasitemia = parasitemia$parasitemia.proxy)

    p = ggplot(ct, aes(parasitemia, counts)) +
        geom_point() +
        geom_smooth(method="lm") +
        xlab("Inferred parasitemia") +
        ylab("Normalized count per million reads") +
        theme_bw()
    ggsave(p, file=paste0("reports/cpm_by_parasitemia.", gene.name, ".pdf"))
}

lapply(c("HEP_00087700", "HEP_00401500", "HEP_00477000",
         "HEP_00319200", "HEP_00421400", "HEP_00390600",
         "HEP_00097800"), plot.gene.expr)

# --- Make volcano plot

keygenes = c("HEP_00087700", "HEP_00401500", "HEP_00477000",
             "HEP_00319200", "HEP_00421400", "HEP_00390600",
             "HEP_00097800")

p = ggplot(tt$table, aes(logFC, -log10(PValue), col=abs(logFC))) +
    geom_point(size=0.5) +
    geom_text_repel(data=tt$table[tt$table$GeneID %in% keygenes & tt$table$logFC > 0,],
        aes(label = GeneID), size = 2, color="black",
        force = 10, min.segment.length = 0,
        segment.size=0.1,
        nudge_x=300, nudge_y=0.25) +
    geom_text_repel(data=tt$table[tt$table$GeneID %in% keygenes & tt$table$logFC < 0,],
        aes(label = GeneID), size = 2, color="black",
        force = 10, min.segment.length = 0,
        segment.size=0.1,
        nudge_x=-300, nudge_y=0.25) +
    geom_point(data=tt$table[tt$table$PValue < 0.05,],
        aes(logFC, -log10(PValue)), pch=21, fill="#CC0033", size=1.5) +
    xlab(expression(paste(log[2], "(Fold Change)"))) +
    ylab(expression(paste(-log[10], "(unadj. p-value)"))) +
    scale_color_gradient(
        low = "grey25",
        high = "grey40",
    ) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "#4D0000"),
          axis.text.x = element_text(color = "grey45"),
          axis.text.y = element_text(color = "grey45"))

ggsave(p, file="reports/volcano_Hepato_genes.pdf",
    height=3, width=4)

save.image("after_DE_analysis_parasite.Rdata")
