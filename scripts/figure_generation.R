#!/usr/bin/env Rscript 


################## FIGURE 1 and 2 --> Stacked Bar Plot of Parasitemia and Volcano Plot with GOI Labeled #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/Fig1_and_2.RData")
# ----- Pre-processing. All of these defined in RData file
library(stringr)
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR/")

# Create Parasitemia dataframe formated for stacked barplot
parasitemia=read.csv("parasitemia_proxy.csv")
parasitemia$Sample_name=str_sub(parasitemia$Sample_name, 1,5)
parasitemia=subset(parasitemia, select=c(Sample_name, Colobus_Reads_gene, Hepatocystis_Reads_gene, Total_Reads_gene))
parasitemia$percent.colobus=parasitemia$Colobus_Reads_gene/parasitemia$Total_Reads_gene
parasitemia$percent.hepatocystis=parasitemia$Hepatocystis_Reads_gene/parasitemia$Total_Reads_gene

parasitemia.colobus=subset(parasitemia, select=c(Sample_name))
parasitemia.colobus$count=parasitemia$Colobus_Reads_gene
parasitemia.colobus$percent=parasitemia$percent.colobus
parasitemia.colobus$type="Colobus"
parasitemia.hepatocystis=subset(parasitemia, select=c(Sample_name))
parasitemia.hepatocystis$count=parasitemia$Hepatocystis_Reads_gene
parasitemia.hepatocystis$percent=parasitemia$percent.hepatocystis
parasitemia.hepatocystis$type="Hepatocystis"

stacked = rbind(parasitemia.colobus, parasitemia.hepatocystis)

# Fig. 1 Stacked Bar plot
hepato = data.frame(Sample_name   = parasitemia$Sample_name,
                    hepato_reads  = parasitemia$Hepatocystis_Reads_gene,
                    colobus_reads = parasitemia$Colobus_Reads_gene)

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

library(ggplot2)
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

ggsave(p, file="../Sub2_figures/Figure_1.pdf",
       height=4, width=4)

# Format Gene lists for volcano plot
colobus_genes = read.csv("Results_gene/colobus_bg_genes.csv")

library(stringr)
colobus_genes = colobus_genes[which(str_sub(colobus_genes$GeneID,1,3) != "LOC"),]

for (i in 1:length(colobus_genes$GeneID)){
  if(abs(colobus_genes$logFC[i])>2 && abs(-log10(colobus_genes$FDR[i]))>1) {colobus_genes$color[i] = "blue"} 
  else if (abs(colobus_genes$logFC[i])>2) {colobus_genes$color[i] = "black"} 
  else if (abs(-log10(colobus_genes$FDR[i]))>1) {colobus_genes$color[i] = "red"} 
  else {colobus_genes$color[i] = "gray"}}

genesnow = c("PAPPA2", "ABCB11", "FBXW11", "TMEM132D", "TBC1D9B", "HRH2", "NAGA")
for (i in 1:length(colobus_genes$GeneID)){
  if(colobus_genes$GeneID[i] %in% genesnow) {colobus_genes$Name[i] = colobus_genes$GeneID[i]} 
  else {colobus_genes$Name[i] = ""}}

bloodgenes = c("RHAG", "SPTA1", "KLF1", "ABCB6", "SLC4A1", "SLC2A1", "STEAP3", "JAK2")

# ----- Plot
# Fig. 2 Volcano plot
p = ggplot(colobus_genes,
           aes(logFC, -log10(FDR),
               size  = -log10(FDR))) +
  geom_point(aes(logFC, -log10(FDR), color = color), size=1) +
  geom_point(data=colobus_genes[colobus_genes$GeneID %in% bloodgenes,],
             aes(logFC, -log10(FDR)), pch=23, col="#4D0000", fill="white", size=1.5) +
  ggrepel::geom_text_repel(data=colobus_genes, aes(label = Name), size = 5, color="black", force = 10, min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf) +
  xlab("Fold Change") +
  ylab("-log10(adjusted p-value)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("#CCCCCC", "#CC0033", "black", "grey45")) +
  scale_shape_manual(values = c(19, 8))

ggsave(p, file="../Sub2_figures/Figure_2.pdf",
       height=3, width=4)

################## FIGURE 3, 4, and S5 --> Up- and Down-regulated genes Genes of Interest for colobus and ACKR1 (DARC) #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/Figure3_and_4_andS5.RData")
# ----- Pre-processing. All of these defined in RData file
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR/")
load("Results_gene/colobus.fc.Rdata")

Colobusfc = colobus.fc

# Normalize and remove low count genes
fc.dge = DGEList(counts=Colobusfc$counts, genes=Colobusfc$annotation)

# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Get CPM just normalized by library composition
cpm = cpm(fc.dge.norm)

# Load parasitemia info
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

# up-regulated - Figure 3
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

# Down-regulated - Figure 4
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

# ACKR1 - Figure S5
cpm.gene_of_interest.ACKR1 = data.frame(cpm[row.names(cpm) == "ACKR1",])
colnames(cpm.gene_of_interest.ACKR1)[1] = "CPM"
cpm.gene_of_interest.ACKR1$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.ACKR1$gene = "ACKR1"
cpm.gene_of_interest.ACKR1$Individuals = rownames(cpm.gene_of_interest.ACKR1)
cpm.gene_of_interest.ACKR1$Individuals=str_sub(cpm.gene_of_interest.ACKR1$Individuals,1,5)

# ----- Plot
# Fig. 3 Up-regulated
up = ggplot(cpm.gene_of_interest, aes(parasitemia, CPM)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  facet_wrap(~ gene, scales="free") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(up, file="../Sub2_figures/Figure_3.pdf",
       height=5, width=5)

# Fig. 4 Down-regulated
dn = ggplot(cpm.gene_of_interest_down, aes(parasitemia, CPM)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  facet_wrap(~ gene, scales="free") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(dn, file="../Sub2_figures/Figure_4.pdf",
       height=3, width=9)

# Fig. S5 ACKR1 - With outlier
S5.out = ggplot(cpm.gene_of_interest.ACKR1, aes(parasitemia, CPM)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(S5.out, file="../Sub2_figures/Figure_S5out.pdf",
       height=4, width=4)


# Fig. S5 ACKR1 - withoutoutlier
S5.without = ggplot(cpm.gene_of_interest.ACKR1, aes(parasitemia, CPM)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  ylim(0, 7.5) +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  ggtitle("ACKR1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(S5.without, file="../Sub2_figures/Figure_S5without.pdf",
       height=4, width=4)

################## FIGURE 5 and S6 --> Functional Enrichment Bar plot, and Functional Enrichment including HP terms #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/Figure5_and_S6.RData")
library(dplyr)
# ----- Pre-processing. All of these defined in RData file
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR")

# Overrep
gmt.up1=read.table("Results_gene/malaria_GMT_overrep.upreg_mod1.txt")
colnames(gmt.up1)=c("GMT", "P_Val")
gmt.up1$log_p_val=-log10(gmt.up1$P_Val)
gmt.up1 <- gmt.up1[order(gmt.up1$log_p_val), ]
gmt.up1$GMT <- factor(gmt.up1$GMT, levels = gmt.up1$GMT[order(gmt.up1$log_p_val)])

# Functional Enrichment 
Functional_enrichment = read.csv("Results_gene/Functional_enrichment_terms.csv")

Fig_5 = subset(Functional_enrichment, Type!="HP")

Fig_5 <- Fig_5[order(Fig_5$adj.P.Val), ]
Fig_5$neg_log_adj.P.Val = -log(Fig_5$adj.P.Val)
Fig_5$Function <- factor(Fig_5$Function, levels = Fig_5$Function[order(Fig_5$neg_log_adj.P.Val)])
Fig_5$color = "blank"

Fig_5 = mutate(Fig_5, color = ifelse(Type %in% "CC", "#CC6666",
                                     ifelse(Type %in% "REAC", "#993333",
                                            ifelse(Type %in% "BP", "#990000",
                                                   ifelse(Type %in% "MF", "#660000", "blank")))))

# Functional Enrichment including HP terms 
Fig_S6 = subset(Functional_enrichment)
Fig_S6$neg_log_adj.P.Val = -log(Functional_enrichment$adj.P.Val)
Fig_S6$Function <- factor(Fig_S6$Function, levels = Fig_S6$Function[order(Fig_S6$neg_log_adj.P.Val)])


# ----- Plot
# Fig 5. Functional enrichment
cols <- c("CC" = "#CC6666", "REAC" = "#993333", "BP" = "#990000", "MF" = "#000000")

FE = ggplot(Fig_5, aes(Function, neg_log_adj.P.Val, fill=Type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cols) +
  ylab("-log(adj. P-val)") +
  coord_flip() +
  theme_classic() 

ggsave(FE, file="../Sub2_figures/Figure_5.pdf",
       height=2, width=8)

# Fig S6 Functional enrichment including HP terms
cols <- c("CC" = "#CC6666", "REAC" = "#993333", "BP" = "#990000", "MF" = "#000000", "HP" = "#666666")

S6 = ggplot(Fig_S6, aes (Function, neg_log_adj.P.Val, fill=Type)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cols) +
  ylab("-log(adj. P-val)") +
  coord_flip() +
  theme_classic() 

ggsave(S12, file="../Sub2_figures/FigureS6.pdf",
       height=2, width=8)


################## FIGURE 6a and b --> Immune composition #####################
rm(list=ls())
load("Results_gene/immune_comp_stopping_point.Rdata")
# ----- Pre-processing. All of these defined in RData file
library(stringr)
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR/")

#--> Cell type proportions against parasitemia
# Monocytes - Figure 6a

Mono = ggplot(parasitemia.immune, aes(parasitemia.proxy, Monocytes)) +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  xlab("Parasitemia Proxy") +
  ylab("Monocyte Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(Mono, file="../Sub2_figures/Figure_6a.pdf",
       height=2.5, width=2.5)

# Memory Activated CD4+ T cells - Figure 6b

MemAct = ggplot(parasitemia.immune, aes(parasitemia.proxy, MemactCD4plus)) +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  xlab("Parasitemia Proxy") +
  ylab("Memory Activated\nCD4+ T-Cell Proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45"))

ggsave(MemAct, file="../Sub2_figures/Figure_6b.pdf",
       height=2.5, width=2.5)


################## FIGURE S1a-b --> Comparing published reads to ours #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/FigureS1a_b.RData")
# ----- Pre-processing. All of these defined in RData file

# Compare hepato read counts for Aunin and EdgeR (gene) - Figure S1a
options(scipen=0)
edgeR_gene_read_fit = lm(log10(aunin_hep_ct) ~ log10(edgeR_gene_hep_count), data = comparasitemia)
summary(edgeR_gene_read_fit)

a = ggplot(comparasitemia, aes(x=edgeR_gene_hep_count, y=aunin_hep_ct)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  ggtitle("EdgeR Hepatocystis Read Counts") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45")) +
  xlab("Our Read Counts") +
  ylab("Published Read Counts")

ggsave(a, file="FigureS1a.pdf",
       height=4, width=4)

# Compare proxy for Aunin and EdgeR (gene) - Figure S1b
edgeR_gene_proxy_fit = lm(Aunin_edgeR_gene_proxy ~ edgeR_gene_proxy, data = comparasitemia)
summary(edgeR_gene_proxy_fit)

b = ggplot(comparasitemia, aes(x=edgeR_gene_proxy, y=Aunin_edgeR_gene_proxy)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  ggtitle("EdgeR Inferred Parasitemia") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45")) + 
  xlab("Our Inferred Parasitemia") +
  ylab("Aunin Inferred Parasitemia")

ggsave(b, file="FigureS1c.pdf",
       height=4, width=4)


################## FIGURE S2, S3, S4, S8 --> Histograms of DE genes for Models 1 and 2 #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/Figure_S2_S3_S4_S8.RData")
# ----- Pre-processing. All of these defined in RData file

# Plot histograms using ggplot
# Figure S2
top.table.df1 = data.frame(tt1$table)
hist1 = ggplot(tt1$table, aes(x=tt1$table$FDR)) + 
  geom_histogram(color="black", fill="#CC0033", bins=20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.78, 0.93),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) +
  xlab("Adj. P-value (FDR)") +
  ylab("Count") +
  ggtitle("Model 1")

ggsave(hist1, file="FigureS2.pdf",
       height=4, width=4)

# Figure S3
hist2 = ggplot(tt2$table, aes(x=tt2$table$FDR)) + 
  geom_histogram(color="black", fill="#CC0033", bins=20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.78, 0.93),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) +
  xlab("Adj. P-value (FDR)") +
  ylab("Count") +
  ggtitle("Model 2")

ggsave(hist2, file="FigureS3.pdf",
       height=4, width=4)

# Figure S4
hist3 = ggplot(tt2.sex$table, aes(x=tt2.sex$table$FDR)) + 
  geom_histogram(color="black", fill="#CC0033", bins=20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.78, 0.93),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) +
  xlab("Adj. P-value (FDR)") +
  ylab("Count") +
  ggtitle("Model 2 - Sex")

ggsave(hist3, file="FigureS4.pdf",
       height=4, width=4)

# Figure S8
hist4 = ggplot(tt3.PC1$table, aes(x=tt3.PC1$table$FDR)) + 
  geom_histogram(color="black", fill="#CC0033", bins=20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.78, 0.93),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) +
  xlab("Adj. P-value (FDR)") +
  ylab("Count") +
  ggtitle("Model 3 - PC1")

ggsave(hist4, file="../Sub2_figures/FigureS8.pdf",
       height=4, width=4)

################## Figure S7 --> PCA #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/Figure_S7.RData")
# ----- Pre-processing. All of these defined in RData file

#Graph of individuals
fviz_pca_ind(pca_res,
             col.ind = "#CC0033", # Color by the contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             geom="point",
             repel = TRUE # Avoid text overlapping
) +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.78, 0.93),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) 


#by quality of cells
PCA_cell = fviz_pca_ind(pca_cell,
                        col.ind = "contrib", # Color by the quality of representation
                        gradient.cols = c("grey45", "black", "#CC0033"),
                        repel = TRUE # Avoid text overlapping
) +
  labs(x = "PC1", y = "PC2", title = "Cells - PCA") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text.align = 0) 

ggsave(PCA_cell, file="../Sub2_figures/FigureS7.pdf",
       height=6, width=6)  

################## FIGURE S9 --> Up- and Down-regulated genes of interest for Hepatocystis #####################
load("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/Sub2_figures/FigureS9.RData")
# ----- Pre-processing. All of these defined in RData file
setwd("~/Desktop/Dissertation/Coding/Colobus_project/Colobus_paper/Third_analysis/edgeR")
load("Results_gene/hepato.fc.Rdata")

fc = hepato.fc

# Load parasitemia info
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

parasitemia$parasitemia.proxy =
  parasitemia$Hepatocystis_Reads_Mapped /
  parasitemia$Colobus_Reads_Mapped

# Normalize and remove low count genes
fc.dge = DGEList(counts=fc$counts, genes=fc$annotation)

# Filter for genes with low counts across conditions
keep = rowSums(cpm(fc.dge) > 5) >= 2

fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

# Normalize for library compositional bias
fc.dge.norm  = calcNormFactors(fc.dge)
fc.dge.norm$samples

# Get CPM just normalized by library composition
cpm = cpm(fc.dge.norm)

# Look at genes of interest - up and down
cpm.gene_of_interest.HEP_00133100 = data.frame(cpm[row.names(cpm) == "HEP_00133100",])
colnames(cpm.gene_of_interest.HEP_00133100)[1] = "CPM"
cpm.gene_of_interest.HEP_00133100$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00133100$gene = "HEP_00133100"

cpm.gene_of_interest.HEP_00502900 = data.frame(cpm[row.names(cpm) == "HEP_00502900",])
colnames(cpm.gene_of_interest.HEP_00502900)[1] = "CPM"
cpm.gene_of_interest.HEP_00502900$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00502900$gene = "HEP_00502900"

cpm.gene_of_interest.HEP_00486000 = data.frame(cpm[row.names(cpm) == "HEP_00486000",])
colnames(cpm.gene_of_interest.HEP_00486000)[1] = "CPM"
cpm.gene_of_interest.HEP_00486000$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00486000$gene = "HEP_00486000"

cpm.gene_of_interest.HEP_00457900 = data.frame(cpm[row.names(cpm) == "HEP_00457900",])
colnames(cpm.gene_of_interest.HEP_00457900)[1] = "CPM"
cpm.gene_of_interest.HEP_00457900$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00457900$gene = "HEP_00457900"

cpm.gene_of_interest.HEP_00452800 = data.frame(cpm[row.names(cpm) == "HEP_00452800",])
colnames(cpm.gene_of_interest.HEP_00452800)[1] = "CPM"
cpm.gene_of_interest.HEP_00452800$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00452800$gene = "HEP_00452800"

cpm.gene_of_interest.HEP_00534450 = data.frame(cpm[row.names(cpm) == "HEP_00534450",])
colnames(cpm.gene_of_interest.HEP_00534450)[1] = "CPM"
cpm.gene_of_interest.HEP_00534450$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00534450$gene = "HEP_00534450"

cpm.gene_of_interest.HEP_00376300 = data.frame(cpm[row.names(cpm) == "HEP_00376300",])
colnames(cpm.gene_of_interest.HEP_00376300)[1] = "CPM"
cpm.gene_of_interest.HEP_00376300$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00376300$gene = "HEP_00376300"

cpm.gene_of_interest.HEP_00525100 = data.frame(cpm[row.names(cpm) == "HEP_00525100",])
colnames(cpm.gene_of_interest.HEP_00525100)[1] = "CPM"
cpm.gene_of_interest.HEP_00525100$parasitemia = parasitemia$parasitemia.proxy
cpm.gene_of_interest.HEP_00525100$gene = "HEP_00525100"

cpm.gene_of_interest = rbind(cpm.gene_of_interest.HEP_00133100, cpm.gene_of_interest.HEP_00502900, cpm.gene_of_interest.HEP_00486000, cpm.gene_of_interest.HEP_00457900, cpm.gene_of_interest.HEP_00452800, cpm.gene_of_interest.HEP_00534450, cpm.gene_of_interest.HEP_00376300, cpm.gene_of_interest.HEP_00525100)
cpm.gene_of_interest$Individual = rownames(cpm.gene_of_interest)
cpm.gene_of_interest$Individual=str_sub(cpm.gene_of_interest$Individual,1,5)

# ----- Plot 
# Up- and Down-regulated GOI for Hepatocystis - Figure S9
hep = ggplot(cpm.gene_of_interest, aes(parasitemia, CPM)) +
  geom_point(pch=21, fill="#CC0033", col="black") +
  geom_smooth(method="lm", se=FALSE, col="#4D0000", lwd=0.25) +
  facet_wrap(~ gene, scales="free") +
  xlab("Inferred parasitemia proxy") +
  ylab("Gene Expression\n(Normalized count per million reads)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#4D0000"),
        axis.text.x = element_text(color = "grey45"),
        axis.text.y = element_text(color = "grey45")) +
  geom_text(data=subset(cpm.gene_of_interest, Individual == "RC127"),
            aes(parasitemia,CPM,label=Individual), vjust=2, hjust=0)

ggsave(hep, file="../Sub2_figures/FigureS9.pdf",
       height=6, width=6)
