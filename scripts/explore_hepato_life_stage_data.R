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

# --- Bring in Hepatocystis life stage proportion data

hepato.stage = read.csv("data/Aunin_allstagecell_props.csv",
    header=TRUE)
names(hepato.stage)[1] = "stage"

hepato.stage = data.frame(names(hepato.stage), t(hepato.stage))
names(hepato.stage) = hepato.stage[1,]
hepato.stage = hepato.stage[-c(1),]
names(hepato.stage)[1] = "Sample"

for (life.stage in names(hepato.stage)[-c(1)]) {
    hepato.stage[[life.stage]] = as.numeric(hepato.stage[[life.stage]])
}

# Combine miscellaneous columsn to match Aunin et al.
other.cols = c("bbSpz", "sgSpz", "EEF", "oocyst", "ook", "ookoo")
hepato.stage$Other = rowSums(hepato.stage[,other.cols])
hepato.stage = hepato.stage[,!(names(hepato.stage) %in% other.cols)]

hepato.stage = hepato.stage[order(hepato.stage$Sample),]

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

parasitemia.immune = merge(parasitemia, hepato.stage,
    by.x = "Sample_name", by.y = "Sample")

life.stages = names(hepato.stage)[-c(1)]

# --- Test for correlation of life stages with parasitemia

para.imm.high = parasitemia.immune[
    parasitemia.immune$parasitemia.proxy >
        quantile(parasitemia.immune$parasitemia.proxy, probs = c(0.05)),]

cor.test.res = lapply(life.stages, function (life.stage) {
    cor.test(para.imm.high$parasitemia.proxy,
             para.imm.high[[life.stage]])
})
names(cor.test.res) = life.stages

cor.test.res.simp = data.frame(do.call(rbind, lapply(life.stages, function (life.stage) {
    x = cor.test.res[[life.stage]]
    cbind(life.stage, x$estimate, x$p.value)
})))
names(cor.test.res.simp) = c("life.stage", "cor", "cor.p")

write.table(cor.test.res.simp, file="reports/hepatocystis_life_stage_parasitemia_cor.txt",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# cor.test.res.simp[cor.test.res.simp$V3 < 0.05,]

# Plot cell type against parasitemia

tmp = lapply(life.stages, function (life.stage) {

    p = ggplot(para.imm.high, aes_string("parasitemia.proxy", life.stage)) +
            geom_point() +
            geom_smooth(method="lm") +
            xlab("Inferred parasitemia") +
            ylab(paste("Hepatocystis life stage proportion: ", life.stage)) +
            theme_bw()
    ggsave(p, file=paste0("reports/life_stage_prop_by_parasitemia.", life.stage, ".pdf"))

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

# --- All stages, all interactions: Estimate dispersion and ID DE genes

# Ensure samples are in correct order
table(gsub(".colobus.bam", "", row.names(fc.dge.norm$samples)) ==
    parasitemia.immune$Sample)

design.c = model.matrix(~ parasitemia.proxy + Female + Male +
                        Merozoite + Ring + Schizont + Trophozoite + Other +
                        parasitemia.proxy:Female + parasitemia.proxy:Male +
                        parasitemia.proxy:Merozoite + parasitemia.proxy:Ring +
                        parasitemia.proxy:Trophozoite,
                        data=parasitemia.immune)

# colnames(design)

disp.c = estimateDisp(fc.dge.norm, design.c, robust = TRUE)

# CPM
cpm.disp.c = cpm(disp.c)

fit.c = glmQLFit(disp.c, design.c, robust = TRUE)

tt.ls = lapply(colnames(design.c)[-c(1)], function (life.stage) {

    qlf.c = glmQLFTest(fit.c, coef = life.stage)
    tt.c  = topTags(qlf.c, n=Inf, adjust.method = "BH", p.value = 1)
    cbind(life.stage, tt.c$table[, c(1,7,10:11)])
})

names(tt.ls) = colnames(design.c)[-c(1)]

# Get count of sig DE genes for each category
lapply(tt.ls, function (x) { nrow(x[x$FDR < 0.2,]) })

save.image("after_Hepato_life_stage_DE.Rdata")
