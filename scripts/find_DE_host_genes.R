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

# --- Load feature count data

load("data/Hepato_mapping.colobus.ind.fc.Rdata")

fc = hepato.ind.fc

# Show number of genes and individuals
# dim(fc$counts)

# --- Load parasitemia info

# From Google Sheet. Download > CSV for sheet "Merged file read counts"
parasitemia = read.csv("data/parasitemia_proxy.csv", header=TRUE)

parasitemia$parasitemia.proxy =
    parasitemia$Hepatocystis_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts=fc$counts, genes=fc$annotation)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge) > 5) >= 4

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
table(gsub(".colobus.bam", "", row.names(fc.dge.norm$samples)) == parasitemia$Sample)

design = model.matrix(~ parasitemia.proxy, data=parasitemia)

# colnames(design)

disp = estimateDisp(fc.dge.norm, design, robust = TRUE)

# CPM
cpm.disp = cpm(disp)

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

# --- Make Table for paper of sig up-reg genes

sig.up.reg.tbl = tt$table[tt$table$FDR < 0.05 & tt$table$logFC > 0, c(1,7,10:11)]

write.table(sig.up.reg.tbl, "reports/sig_up_reg_genes.txt",
    quote=FALSE, row.names=FALSE)

# --- Make Table for suppl of up- and down-reg genes

sig.reg.tbl = rbind(tt$table[tt$table$FDR < 0.1 & tt$table$logFC > 0, c(1,7,10:11)],
                    tt$table[tt$table$FDR < 0.1 & tt$table$logFC < 0, c(1,7,10:11)])

write.table(sig.reg.tbl, "reports/sig_reg_genes.txt",
    quote=FALSE, row.names=FALSE)

save.image("tmp_after_analysis.Rdata")

# --- Look up DARC gene, using new name

sink(file="reports/DARC_info.txt", split=TRUE)
    tt$table[tt$table$GeneID == "ACKR1",c(1,6:11)]
sink()

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

lapply(c("ACKR1", "LSM14A", "APOBEC2", "UBE2K", "CPNE3"), plot.gene.expr)

# --- Do functional overrep test using gProfiler - GO and HP, etc. annotations

# Up-regulated

# Get temporary link to results
(sig.up.link = gost(sig.genes.up,
                organism = "ptephrosceles", ordered_query = FALSE,
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.1, correction_method = "fdr",
                domain_scope = "custom",
                custom_bg = bg.genes,
                numeric_ns = "ENTREZGENE_ACC",
                sources = c("GO:MF", "GO:CC", "GO:BP", "HP"),
                as_short_link = TRUE))

sig.up.gp = gost(sig.genes.up,
                organism = "ptephrosceles", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.1, correction_method = "fdr",
                domain_scope = "custom",
                custom_bg = bg.genes,
                numeric_ns = "ENTREZGENE_ACC",
                sources = c("GO:MF", "GO:CC", "GO:BP", "HP"),
                as_short_link = FALSE)

sig.up.gp.simp = sig.up.gp$result

sig.up.gp.simp = sig.up.gp.simp[order(sig.up.gp.simp$p_value),
                    c("source", "term_id", "term_name",
                      "term_size", "query_size", "intersection_size", "p_value")]

# Include only terms with more than 2 genes
sig.up.gp.simp.flt = sig.up.gp.simp[sig.up.gp.simp$intersection_size > 2,]

# Down-regulated

sig.dn.gp = gost(sig.genes.dn,
                organism = "ptephrosceles", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.1, correction_method = "fdr",
                domain_scope = "custom",
                custom_bg = bg.genes,
                numeric_ns = "ENTREZGENE_ACC",
                sources = c("GO:MF", "GO:CC", "GO:BP", "HP"),
                as_short_link = FALSE)

sig.dn.gp.simp = sig.dn.gp$result

sig.dn.gp.simp = sig.dn.gp.simp[order(sig.dn.gp.simp$p_value),
                    c("source", "term_id", "term_name",
                      "term_size", "query_size", "intersection_size", "p_value")]

# Include only terms with more than 2 genes
sig.dn.gp.simp.flt = sig.dn.gp.simp[sig.dn.gp.simp$intersection_size > 2,]

# Write tables

write.table(sig.up.gp.simp.flt, file="reports/GO_HP_overrep.upreg.txt",
    sep="\t", row.names=FALSE)
write.table(sig.dn.gp.simp.flt, file="reports/GO_HP_overrep.dnreg.txt",
    sep="\t", row.names=FALSE)

# --- Get DE p-values for blood disorder genes

sig.blood.genes = c("RHAG", "SPTA1", "KLF1", "ABCB6", "SLC4A1")

sink(file="reports/blood_gene_info.txt", split=TRUE)
    tt$table[row.names(tt$table) %in% sig.blood.genes, c(1,7,10:11)]
sink()

# --- Do functional enrichment test on GMT files - Convert genes to human

# Up-regulated

genes.up.human = gorth(
    sig.genes.up,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

genes.up.hs = genes.up.human$ortholog_ensg[genes.up.human$ortholog_ensg != "N/A"]

# Down-regulated

genes.dn.human = gorth(
    sig.genes.dn,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

genes.dn.hs = genes.dn.human$ortholog_ensg[genes.dn.human$ortholog_ensg != "N/A"]

# Background

bg.genes.human = gorth(
    bg.genes,
    source_organism = "ptephrosceles",
    target_organism = "hsapiens",
    numeric_ns = "ENTREZGENE_ACC",
    mthreshold = Inf,
    filter_na = TRUE
)

bg.genes.hs = bg.genes.human$ortholog_ensg[bg.genes.human$ortholog_ensg != "N/A"]

# --- Test for significant overlap between GMTs and up- or down-regulated genes
# --- using hypergeometric test

gmt.files = c("data/gene_sets/CTD_gene_disease_malaria.gmt",
              "data/gene_sets/DISEASE_experimental_malaria.gmt",
              "data/gene_sets/DISEASE_textmining_malaria.gmt",
              "data/gene_sets/GAD_gene_disease_malaria.gmt",
              "data/gene_sets/Ebel_interacting_proteins.PPIPs.gmt",
              "data/gene_sets/Ebel_interacting_proteins.PPIPs-human.gmt",
              "data/gene_sets/HALLMARK_HEME_METABOLISM.gmt",
              "data/gene_sets/REACTOME_ERYTHROCYTES_TAKE_UP_OXYGEN_AND_RELEASE_CARBON_DIOXIDE.gmt",
              "data/gene_sets/REACTOME_ERYTHROCYTES_TAKE_UP_CARBON_DIOXIDE_AND_RELEASE_OXYGEN.gmt",
              "data/gene_sets/STEINER_ERYTHROCYTE_MEMBRANE_GENES.gmt"
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

write.table(hyper.p.up, file="reports/malaria_GMT_overrep.upreg.txt",
    row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(hyper.p.dn, file="reports/malaria_GMT_overrep.downreg.txt",
    row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- Test whether parasitemia-linked genes are enriched for selected genes

# Van der Lee et al. 2017; doi: 10.1093/nar/gkx704, Table S4
sel = read.table("data/gene_sets/Van_der_Lee_et_al__TableS4__331_PSG__info_statistics.txt",
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

write.table(up.reg.sel.test, file="reports/dNdS_selected_overrep.upreg.txt",
    row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(dn.reg.sel.test, file="reports/dNdS_selected_overrep.downreg.txt",
    row.names=FALSE, col.names=FALSE, quote=FALSE)

save.image("after_DE_analysis.Rdata")
