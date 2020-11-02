#!/usr/bin/env Rscript

library(ggplot2)

# install.packages("ggpubr")
library(ggpubr)

library(xtable)

#--> Scatter plot of parasitemia

# Load Rdata object
load("after_DE_analysis.Aunin_med.Rdata")

hepato = data.frame(Sample_name   = parasitemia$Sample_name,
                    hepato_reads  = parasitemia$Aunin_hep_ct,
                    colobus_reads = parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping)
hepato$col_hunthou = hepato$colobus_reads / 100000
hepato$hep_per_hunthoucol = hepato$hepato_reads / hepato$col_hunthou

p = ggplot(hepato, aes(x=Sample_name, y=hep_per_hunthoucol)) +
        geom_point() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        xlab("Sample Name") +
        ylab("Hepatocystis Reads per 100,000 Colobus Reads")

hepato$total_reads = hepato$hepato_reads + hepato$colobus_reads

hepato.l = rbind(data.frame(Sample_name = as.character(hepato$Sample_name),
                            Reads = as.numeric(hepato$hepato_reads),
                            read_type="Hepato"),
                 data.frame(Sample_name = as.character(hepato$Sample_name),
                            Reads = as.numeric(hepato$colobus_reads),
                            read_type="Colobus"))

hepato.l$Sample_name = factor(hepato.l$Sample_name,
    levels = hepato$Sample_name[order(-1 * hepato$total_reads)])

hepato.l$read_type = factor(hepato.l$read_type, levels=c("Colobus", "Hepato"))

p = ggplot(hepato.l, aes(fill=read_type, y=Reads / 10^6, x=Sample_name)) +
    geom_bar(position="stack", stat="identity") +
    xlab("Blood Sample") +
    ylab("Uniquely Mapped Sequencing Read Count (10^6)") +
    ylab(expression(paste("Uniquely Mapped Sequencing Read Count (", 10^6, ")"))) +
    scale_fill_manual(values = c("#CC0033", "black"),  # "#5F6A72"
                      breaks = c("Colobus", "Hepato"),
                      labels = c("Red Colobus Monkey", expression(italic("Hepatocystis")))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position = c(0.78, 0.93),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.text.align = 0)

ggsave(p, file="reports/parasitemia_barplot.pdf",
    height=4, width=4)

#--> Cell type proportions against parasitemia

rm(list=ls())
load("after_cell_comp_DE.Rdata")

# Monocytes

p = ggplot(parasitemia.immune, aes(parasitemia.proxy, Monocytes)) +
        geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
        geom_point(col="#CC0033") +
        xlab("Parasitemia Proxy") +
        ylab("Monocyte Proportion") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "#4D0000"),
              axis.text.x = element_text(color = "grey45"),
              axis.text.y = element_text(color = "grey45"))

ggsave(p, file="reports/fancy_cell_type_prop_by_parasitemia.Monocytes.pdf",
    height=2.5, width=2.5)

# Neutrophils

p = ggplot(parasitemia.immune, aes(parasitemia.proxy, Neutrophils)) +
        geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
        geom_point(col="#CC0033") +
        xlab("Parasitemia Proxy") +
        ylab("Neutrophil Proportion") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "#4D0000"),
              axis.text.x = element_text(color = "grey45"),
              axis.text.y = element_text(color = "grey45"))

ggsave(p, file="reports/fancy_cell_type_prop_by_parasitemia.Neutrophils.pdf",
    height=2.5, width=2.5)

#--> Expression of genes of interest by parasitemia

rm(list=ls())
load("after_DE_analysis.Aunin_med.Rdata")

library(edgeR)

fc.dge = DGEList(counts = fc$counts, genes = fc$annotation)
keep = (rowSums(cpm(fc.dge) > 5) >= 4)
fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]
fc.dge.norm  = calcNormFactors(fc.dge)

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

design = model.matrix(~ parasitemia.proxy, data=parasitemia)
disp = estimateDisp(fc.dge.norm, design, robust = TRUE)
cpm.disp = cpm(disp)

plot.gene.expr = function(gene.name) {
    cpm.disp.gene_of_interest = cpm.disp[row.names(cpm.disp) == gene.name,]
    ct = data.frame(counts      = cpm.disp.gene_of_interest,
                    parasitemia = parasitemia$parasitemia.proxy.aunin_med)

    # Remove outlier in ACKR1
    if (gene.name == "ACKR1") {
        ct = ct[ct$counts < 30,]
    }

    p = ggplot(ct, aes(parasitemia, counts)) +
            geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
            geom_point(pch=21, fill="#CC0033", col="black") +
            ggtitle(bquote(paste(italic(.(gene.name))))) +
            xlab("Parasitemia Proxy") +
            ylab("Gene expression\nNormalized count per million reads") +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(color = "#4D0000"),
                  axis.text.x = element_text(color = "grey45"),
                  axis.text.y = element_text(color = "grey45"))

    ggsave(p, file=paste0("reports/cpm_by_parasitemia.", gene.name, ".Aunin_med.pdf"),
        height=2.5, width=2.5)

    p
}

# Plot ACKR1
expr.plots.ackr1 = lapply(c("ACKR1"), plot.gene.expr)

# Plot up-regulated genes
up.genes = c("UBE2K", "PP2D1", "TMEM167A", "LSM14A",
             "UBE2B", "TTLL12", "AGFG1", "APOBEC2")

expr.plots.up = lapply(up.genes, plot.gene.expr)

names(expr.plots.up) = up.genes

expr.plots.up = lapply(expr.plots.up, function (x) { x = x + rremove("xy.title") })

p.all = ggarrange(plotlist = expr.plots.up,
    ncol = 4, nrow = 2)

p.all = annotate_figure(p.all,
                bottom = text_grob("Parasitemia Proxy"),
                left = text_grob("Gene expression\n(Normalized count per million reads)",
                    rot = 90)
                )

ggsave(p.all, file="reports/cpm_by_parasitemia.up.Aunin_med.pdf",
    height=4, width=8)

# Plot down-regulated genes
dn.genes = c("NAGA", "TBC1D9B", "ZNF682")

expr.plots.dn = lapply(dn.genes, plot.gene.expr)

names(expr.plots.dn) = dn.genes

expr.plots.dn = lapply(expr.plots.dn, function (x) { x = x + rremove("xy.title") })

p.all = ggarrange(plotlist = expr.plots.dn,
    ncol = 3, nrow = 1)

p.all = annotate_figure(p.all,
                bottom = text_grob("Parasitemia Proxy"),
                left = text_grob("Gene expression\n(Normalized count\nper million reads)",
                    rot = 90)
                )

ggsave(p.all, file="reports/cpm_by_parasitemia.down.Aunin_med.pdf",
    height=2.1, width=6.4)

#--> Volcano Plot

rm(list=ls())
load("after_DE_analysis.Aunin_med.Rdata")

library(edgeR)
library(ggrepel)

fc.dge = DGEList(counts = fc$counts, genes = fc$annotation)
keep = (rowSums(cpm(fc.dge) > 5) >= 4)
fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]
fc.dge.norm  = calcNormFactors(fc.dge)

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

design = model.matrix(~ parasitemia.proxy, data=parasitemia)
disp = estimateDisp(fc.dge.norm, design, robust = TRUE)
cpm.disp = cpm(disp)

fit = glmQLFit(disp, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = "parasitemia.proxy")
tt  = topTags(qlf, n=Inf, adjust.method = "BH", p.value = 1)

# Remove outlier
tt$table = tt$table[abs(tt$table$logFC) < 20,]

# Remove LOC genes
tt$table = tt$table[!grepl("^LOC", tt$table$GeneID),]

bloodgenes = c("RHAG", "SPTA1", "KLF1", "ABCB6",
               "SLC4A1", "SLC2A1", "STEAP3", "JAK2")

p = ggplot(tt$table, aes(logFC, -log10(FDR), col=abs(logFC))) +
    geom_point(size=0.5) +
    geom_text_repel(data=tt$table[tt$table$FDR < 0.05 & tt$table$logFC > 0,],
        aes(label = GeneID), size = 2, color="black",
        force = 10, min.segment.length = 0,
        segment.size=0.1,
        nudge_x=3, nudge_y=0.25) +
    geom_text_repel(data=tt$table[tt$table$GeneID == "ACKR1",],
        aes(label = GeneID), size = 2, color="black",
        force = 10, min.segment.length = 0,
        segment.size=0.1,
        nudge_x=1, nudge_y=0.25) +
    geom_text_repel(data=tt$table[tt$table$FDR < 0.05 & tt$table$logFC < 0,],
        aes(label = GeneID), size = 2, color="black",
        force = 10, min.segment.length = 0,
        segment.size=0.1,
        nudge_x=-3, nudge_y=0.25) +
    geom_point(data=tt$table[tt$table$FDR < 0.05,],
        aes(logFC, -log10(FDR)), pch=21, fill="#CC0033", size=1.5) +
    geom_point(data=tt$table[tt$table$GeneID %in% bloodgenes,],
        aes(logFC, -log10(FDR)), pch=23, col="#4D0000", fill="white", size=1.5) +
    xlab(expression(paste(log[2], "(Fold Change)"))) +
    ylab(expression(paste(-log[10], "(FDR)"))) +
    xlim(c(-20,20)) +
    ylim(c(0, 2.25)) +
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

ggsave(p, file="reports/volcano_Aunin_med.pdf",
    height=3, width=4)

#--> Nice table of top DE genes

rm(list=ls())
load("after_DE_analysis.Aunin_med.Rdata")

library(edgeR)
library(ggrepel)

fc.dge = DGEList(counts = fc$counts, genes = fc$annotation)
keep = (rowSums(cpm(fc.dge) > 5) >= 4)
fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]
fc.dge.norm  = calcNormFactors(fc.dge)

parasitemia$parasitemia.proxy = parasitemia$parasitemia.proxy.aunin_med

design = model.matrix(~ parasitemia.proxy, data=parasitemia)
disp = estimateDisp(fc.dge.norm, design, robust = TRUE)
cpm.disp = cpm(disp)

fit = glmQLFit(disp, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = "parasitemia.proxy")
tt  = topTags(qlf, n=Inf, adjust.method = "BH", p.value = 1)

# Remove outlier
tt$table = tt$table[abs(tt$table$logFC) < 20,]

# Remove LOC genes
tt$table = tt$table[!grepl("^LOC", tt$table$GeneID),]

# Write LaTeX table

to.print = tt$table[tt$table$FDR <= 0.05, c(1,7,10:11)]

to.print = to.print[order(to.print$logFC < 1, to.print$FDR),]
names(to.print) = c("Gene", "$\\log_{2}$FC", "$p$-value", "adj. $p$-value")

to.print$Gene = paste0("\\emph{", to.print$Gene, "}")
to.print[,3] = sanitize.numbers(format(to.print[,3], scientific = TRUE, digits=3),
                 type = "latex", math.style.exponents = TRUE)

xt = xtable(to.print, digits=c(0, 0, 2, 0, 4))

out.file = "reports/top_DE_genes_table.tex"

sink(out.file)

cat("\\documentclass[a4paper,landscape]{article}",
    "\\usepackage{graphicx}",
    "\\usepackage{longtable}",
    "\\DeclareGraphicsExtensions{.pdf}",
    "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
    "\\usepackage{caption}",
    "\\captionsetup[table]{labelformat=empty}",
    "\\begin{document}", sep="\n")

print.xtable(xt,
        display = c("s","f","s","f"),
        include.rownames = FALSE,
        sanitize.text.function = function(x){x},
        caption.placement = "top",
        hline.after = c(-1, 0, sum(to.print[,2] > 0), nrow(to.print)))

cat("\\end{document}", sep="\n")

sink()
