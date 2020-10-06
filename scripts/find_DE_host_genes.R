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

# Get median and range of reads mapping to colobus
cat(paste("Median per-individual colobus reads:", median(colSums(fc$counts))), "\n")
cat("Range  per-individual colobus reads:", "\n")
cat(range(colSums(fc$counts)), "\n")

# --- Calculate parasitemia proxy

# From Google Sheet. Download > CSV for sheet "Merged file read counts"
parasitemia = read.csv("data/parasitemia_proxy.csv", header=TRUE)

# Remove RC127
#parasitemia = parasitemia[parasitemia$Sample_name != "RC127",]

parasitemia$parasitemia.proxy.ours =
    parasitemia$Hepatocystis_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping

# --- Compare our read counts to published ones: Hepato reads per individual

aunin = data.frame(t(read.table("data/Aunin_etal_2020_read_counts.csv",
    header=TRUE, row.names=1)))

ind.info = read.table("data/ind_info.txt", header=TRUE)

aunin = merge(aunin, ind.info[,c(1,3)], by.x="row.names", by.y="Sample")

aunin.ind = aggregate(aunin[2:(ncol(aunin) - 1)],
    by=list(aunin$Sample_name), FUN=sum)

names(aunin.ind)[1] = "Sample"

# Also compute sum of middle 50% of genes (toss top and bottom 25% of genes by count)
aunin.gene.sums = colSums(aunin.ind[,-c(1)])
aunin.mid50.bounds = quantile(aunin.gene.sums, probs=c(0.25, 0.75))
to.count = aunin.gene.sums >= aunin.mid50.bounds[1] &
           aunin.gene.sums <= aunin.mid50.bounds[2]
aunin.sums.mid50 = rowSums(aunin.ind[c(FALSE, to.count)])

aunin.hep.ct = data.frame("Sample_name"   = aunin.ind[,1],
                          "Aunin_hep_ct"  = rowSums(aunin.ind[,-c(1)]),
                          "Aunin_hep_med" = aunin.sums.mid50)

# --- Compute parasitemia proxy and compare

parasitemia = merge(parasitemia, aunin.hep.ct, all=TRUE)

parasitemia$parasitemia.proxy.aunin =
    parasitemia$Aunin_hep_ct /
    parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping

parasitemia$parasitemia.proxy.aunin_med =
    parasitemia$Aunin_hep_med /
    parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping

p = ggplot(parasitemia, aes(parasitemia.proxy.aunin, parasitemia.proxy.ours)) +
    geom_point() +
    geom_abline(intercept=0, slope=1)
ggsave(p, file="reports/hepato_reads_per_ind_ours_vs_Aunin_etal.pdf")

parasitemia$parasitemia.proxy.plasmo =
    parasitemia$Plasmodium_Reads_Mapped /
    parasitemia$Colobus_Reads_Mapped_Plasmodium_mapping

p = ggplot(parasitemia, aes(parasitemia.proxy.plasmo, parasitemia.proxy.ours)) +
    geom_point() +
    geom_abline(intercept=0, slope=1)
ggsave(p, file="reports/hepato_reads_per_ind_hepato_vs_plasmo.pdf")

# Plot parasite counts from the Hepato vs. Plasmo mapping
p = ggplot(parasitemia, aes(Plasmodium_Reads_Mapped, Hepatocystis_Reads_Mapped)) +
    geom_point() +
    geom_abline(intercept=0, slope=1)
ggsave(p, file="reports/hepato_reads_raw_count_per_ind_hepato_vs_plasmo.pdf")

# Plot colobus counts from the Hepato vs. Plasmo mapping
p = ggplot(parasitemia, aes(Colobus_Reads_Mapped_Plasmodium_mapping, Colobus_Reads_Mapped_Hepatocystis_mapping)) +
    geom_point() +
    geom_abline(intercept=0, slope=1)
ggsave(p, file="reports/colobus_reads_raw_count_per_ind_hepato_vs_plasmo.pdf")

# --- Compare our read counts to published ones: total reads per Hepato gene

load("data/Hepato_mapping.Hepato.gene.fc.Rdata")

fc.parasite = genehepato.ind.fc

aunin.gene.ct = data.frame("Gene"         = names(aunin.ind)[-c(1)],
                           "Aunin_hep_ct" = colSums(aunin.ind[-c(1)]))
our.gene.ct   = data.frame("Gene"         = row.names(fc.parasite$counts),
                           "our_hep_ct"   = rowSums(fc.parasite$counts))

both.gene.ct = merge(aunin.gene.ct, our.gene.ct, all=TRUE)

both.gene.ct = both.gene.ct[grepl("^HEP_", both.gene.ct$Gene),]

p = ggplot(both.gene.ct, aes(Aunin_hep_ct, our_hep_ct)) +
    geom_point() +
    xlim(c(0,2000000)) + ylim(c(0,2000000)) +
    geom_abline(intercept=0, slope=1)
ggsave(p, file="reports/hepato_reads_per_gene_ours_vs_Aunin_etal.pdf")

# Genes that are only in our dataset (high ones are entirely rRNA)
only.us = both.gene.ct[is.na(both.gene.ct$Aunin_hep_ct),]
only.us = only.us[only.us$our_hep_ct != 0,]

# --- Function to do full DE analysis

do.host.de.analysis = function(in.counts, in.annotation, parasitemia.colname, suffix) {

    parasitemia$parasitemia.proxy = parasitemia[[parasitemia.colname]]

    # --- Normalize and remove low count genes

    # Can add lib.sizes=c() to this to not recalculate
    fc.dge = DGEList(counts = in.counts, genes = in.annotation)

    # Filter for genes with low counts across conditions
    keep = (rowSums(cpm(fc.dge) > 5) >= 4)

    # Make reference table list
    ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))

    fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

    # Normalize for library compositional bias
    fc.dge.norm  = calcNormFactors(fc.dge)
    fc.dge.norm$samples

    # Get CPM just normalized by library composition
    cpm = cpm(fc.dge.norm)

    # --- Estimate dispersion the complicated way (using CR method)
    # --- and ID DE genes

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
    sig.genes.up = tt$table$GeneID[tt$table$FDR < 0.2 & tt$table$logFC > 0]
    length(sig.genes.up)

    sig.genes.dn = tt$table$GeneID[tt$table$FDR < 0.2 & tt$table$logFC < 0]
    length(sig.genes.dn)

    # Write background list of genes
    bg.genes = tt$table$GeneID
    length(bg.genes)

    # Write results to file
    write.table(sig.genes.up,
        file=paste0("reports/edgeR_colobus_upreg_sig_genes.", suffix, ".txt"),
        sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

    write.table(sig.genes.dn,
        file=paste0("reports/edgeR_colobus_dnreg_sig_genes.", suffix, ".txt"),
        sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

    write.table(bg.genes,
        file=paste0("reports/edgeR_colobus_bg_genes.", suffix, ".txt"),
        sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

    # --- Make table for paper of sig up-reg genes

    sig.up.reg.tbl = tt$table[tt$table$FDR < 0.05 & tt$table$logFC > 0, c(1,7,10:11)]

    write.table(sig.up.reg.tbl, paste0("reports/sig_up_reg_genes.", suffix, ".txt"),
        quote=FALSE, row.names=FALSE)

    # --- Make Table for suppl of up- and down-reg genes

    sig.reg.tbl = rbind(tt$table[tt$table$FDR < 0.2 & tt$table$logFC > 0, c(1,7,10:11)],
                        tt$table[tt$table$FDR < 0.2 & tt$table$logFC < 0, c(1,7,10:11)])

    write.table(sig.reg.tbl, paste0("reports/sig_reg_genes.", suffix, ".txt"),
        quote=FALSE, row.names=FALSE)

    # --- Make table of all genes

    full.reg.tbl = tt$table[, c(1,7,10:11)]

    write.table(full.reg.tbl, paste0("reports/full_reg_genes.", suffix, ".txt"),
        quote=FALSE, row.names=FALSE)

    # --- Look up DARC gene, using new name

    darc.info = tt$table[tt$table$GeneID == "ACKR1",c(1,6:11)]
    write.table(darc.info, file=paste0("reports/DARC_info.", suffix, ".txt"))

    # --- Plot genes of interest

    plot.gene.expr = function(gene.name) {

        cpm.disp.gene_of_interest = cpm.disp[row.names(cpm.disp) == gene.name,]
        ct = data.frame(counts      = cpm.disp.gene_of_interest,
                        parasitemia = parasitemia$parasitemia.proxy)

        p = ggplot(ct, aes(parasitemia, counts)) +
            geom_point() +
            geom_smooth(method="lm") +
            xlab("Inferred parasitemia proxy") +
            ylab("Normalized count per million reads") +
            theme_bw()
        ggsave(p, file=paste0("reports/cpm_by_parasitemia.", gene.name, ".", suffix, ".pdf"))
    }

    lapply(c("ACKR1", "UBE2K", "LSM14A", "APOBEC2", "PP2D1"), plot.gene.expr)

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

    write.table(sig.up.gp.simp.flt,
        file=paste0("reports/GO_HP_overrep.upreg.", suffix, ".txt"),
        sep="\t", row.names=FALSE)
    write.table(sig.dn.gp.simp.flt,
        file=paste0("reports/GO_HP_overrep.dnreg.", suffix, ".txt"),
        sep="\t", row.names=FALSE)

    # --- Get DE p-values for blood disorder genes

    # All 3: Stomatocytosis, Poikilocytosis, Reticulocytosis
    sig.blood.genes = c("RHAG", "SPTA1", "SLC4A1", "ABCB6", "SLC2A1")
    # Only Reticulocytosis
    sig.blood.genes = c(sig.blood.genes, "KLF1")
    # Only Poikilocytosis
    sig.blood.genes = c(sig.blood.genes, "STEAP3", "JAK2")

    blood.gene.info = tt$table[row.names(tt$table) %in% sig.blood.genes, c(1,7,10:11)]
    write.table(blood.gene.info, file=paste0("reports/blood_gene_info.", suffix, ".txt"))

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

    gmt.files = c(
        "data/gene_sets/CTD_gene_disease_malaria.gmt",
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

    write.table(hyper.p.up,
        file=paste0("reports/malaria_GMT_overrep.upreg.", suffix, ".txt"),
        row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(hyper.p.dn,
        file=paste0("reports/malaria_GMT_overrep.downreg.", suffix, ".txt"),
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

    write.table(up.reg.sel.test,
        file=paste0("reports/dNdS_selected_overrep.upreg.", suffix, ".txt"),
        row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(dn.reg.sel.test,
        file=paste0("reports/dNdS_selected_overrep.downreg.", suffix, ".txt"),
        row.names=FALSE, col.names=FALSE, quote=FALSE)

    save.image(paste0("after_DE_analysis.", suffix, ".Rdata"))
}

# Our original Hepato mapping
do.host.de.analysis(fc$counts, fc$annotation,
    "parasitemia.proxy.ours",      "hepato_mapping")

# Aunin et al. read counts for computing parasitemia
do.host.de.analysis(fc$counts, fc$annotation,
    "parasitemia.proxy.aunin",     "Aunin_mapping")

# Aunin et al. read counts for computing parasitemia,
# but using middle 50% of genes by expression
do.host.de.analysis(fc$counts, fc$annotation,
    "parasitemia.proxy.aunin_med", "Aunin_med")
