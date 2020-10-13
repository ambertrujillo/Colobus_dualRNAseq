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

# --- Bring in immune cell composition

immune.comp = read.table("data/Simons_etal_2019_immune_components.txt",
    header=TRUE)

immune.comp$Sample[immune.comp$Sample == "RC106"] = "RC106R"

immune.comp = immune.comp[order(immune.comp$Sample),]

# immune.comp$totalCD4 = immune.comp$naiveCD4plus +
#     immune.comp$MemrestCD4plus + immune.comp$MemactCD4plus

# immune.comp$Mono.Neu.ratio = immune.comp$Monocytes / immune.comp$Neutrophils

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

parasitemia.immune = merge(parasitemia, immune.comp,
    by.x = "Sample_name", by.y = "Sample")

# --- Test for correlation of cell types with parasitemia

cell.types = c("MemB", "CD8", "naiveCD4plus",
               "MemrestCD4plus", "MemactCD4plus",
               "restNK", "Monocytes",
               "Macrophages", "Neutrophils")

cor.test.res = lapply(cell.types, function (cell.type) {
    cor.test(parasitemia.immune$parasitemia.proxy, parasitemia.immune[[cell.type]])
})
names(cor.test.res) = cell.types

cor.test.res.simp = data.frame(do.call(rbind, lapply(cell.types, function (cell.type) {
    x = cor.test.res[[cell.type]]
    cbind(cell.type, x$estimate, x$p.value)
})))

write.table(cor.test.res.simp, file="reports/immune_cell_component_parasitemia_cor.txt",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# cor.test.res.simp[cor.test.res.simp$V3 < 0.05,]

# Plot cell type against parasitemia

tmp = lapply(cell.types, function (cell.type) {

    p = ggplot(parasitemia.immune, aes_string("parasitemia.proxy", cell.type)) +
            geom_point() +
            geom_smooth(method="lm") +
            xlab("Inferred parasitemia") +
            ylab(paste("Immune cell type proportion: ", cell.type)) +
            theme_bw()
    ggsave(p, file=paste0("reports/cell_type_prop_by_parasitemia.", cell.type, ".pdf"))

})

# --- Normalize and remove low count genes

# Can add lib.sizes=c() to this to not recalculate
fc.dge = DGEList(counts = fc$counts, genes = fc$annotation)

# Filter for genes with low counts across conditions

keep = rowSums(cpm(fc.dge) > 5) >= 4

# Make reference table list
ref.list = data.frame(Colobus.id = names(keep[keep == TRUE]))

fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# --- All types, no interactions: Estimate dispersion and ID DE genes

# Ensure samples are in correct order
table(gsub(".colobus.bam", "", row.names(fc.dge.norm$samples)) ==
    parasitemia.immune$Sample)

design.c = model.matrix(~ parasitemia.proxy + MemB + CD8 + naiveCD4plus +
                        MemrestCD4plus + MemactCD4plus +
                        restNK + Monocytes +
                        Macrophages + Neutrophils, data=parasitemia.immune)

# colnames(design)

disp.c = estimateDisp(fc.dge.norm, design.c, robust = TRUE)

# CPM
cpm.disp.c = cpm(disp.c)

fit.c = glmQLFit(disp.c, design.c, robust = TRUE)

tt.cc = lapply(colnames(design.c)[-c(1)], function (cell.component) {

    qlf.c = glmQLFTest(fit.c, coef = cell.component)
    tt.c  = topTags(qlf.c, n=Inf, adjust.method = "BH", p.value = 1)
    cbind(cell.component, tt.c$table[, c(1,7,10:11)])
})

names(tt.cc) = c("parasitemia.proxy", cell.types)

# Get count of sig DE genes for each category
lapply(tt.cc, function (x) { nrow(x[x$FDR < 0.1,]) })

# Write parasitemia DE gene results to file
w_cell_comp.parasitemia_DE = tt.cc[[1]][tt.cc[[1]]$FDR < 0.1,]
write.table(w_cell_comp.parasitemia_DE,
    file="results/sig_DE_genes.w_cell_comp.parasitemia.txt",
    sep="\t", row.names=FALSE)

# Examine genes DE for parasitemia in simple model without cell comp
simple.model.sig.genes = read.table("reports/sig_reg_genes.Aunin_med.txt",
    header=TRUE)

simp.and.complex.model.res = merge(simple.model.sig.genes, tt.cc[[1]],
    by="GeneID", suffixes=c(".simplemodel", ".wcellcomp"), sort=FALSE)

top.DE.genes = simp.and.complex.model.res[
    simp.and.complex.model.res$PValue.simplemodel < 0.0001,
]
table(top.DE.genes$PValue.wcellcomp < 0.0001)
table(top.DE.genes$PValue.wcellcomp < 0.001)
table(top.DE.genes$PValue.wcellcomp < 0.01)

# --- Test for interaction between cell type proportion and parasitemia
# --- limited to just cell types with significant correlation to parasitemia

design.inter = model.matrix(~ parasitemia.proxy +
                        Monocytes + Neutrophils +
                        parasitemia.proxy:Monocytes +
                        parasitemia.proxy:Neutrophils, data=parasitemia.immune)

# colnames(design.inter)

disp.inter = estimateDisp(fc.dge.norm, design.inter, robust = TRUE)

# CPM
cpm.disp.inter = cpm(disp.inter)

fit.inter = glmQLFit(disp.inter, design.inter, robust = TRUE)

tt.inter = lapply(colnames(design.inter)[-c(1)], function (cell.component) {

    qlf.inter = glmQLFTest(fit.inter, coef = cell.component)
    tt.inter  = topTags(qlf.inter, n=Inf, adjust.method = "BH", p.value = 1)
    cbind(cell.component, tt.inter$table[, c(1,7,10:11)])
})

names(tt.inter) = colnames(design.inter)[-c(1)]

tt.inter.m = do.call(rbind, tt.inter)
table(tt.inter.m[tt.inter.m$FDR < 0.1,]$cell.component)

save.image("after_cell_comp_DE.Rdata")
