#!/usr/bin/env Rscript

# Load Rdata object

hepato=data.frame(Sample_name = parasitemia$Sample_name, hepato_reads = parasitemia$Aunin_hep_ct, colobus_reads = parasitemia$Colobus_Reads_Mapped_Hepatocystis_mapping)
hepato$col_hunthou = hepato$colobus_reads / 100000
hepato$hep_per_hunthoucol = hepato$hepato_reads / hepato$col_hunthou


#--> Scatter Plot of parasitemia
ggplot(hepato, aes(x=Sample_name, y=hep_per_hunthoucol)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Sample Name") +
  ylab("Hepatocystis Reads/100,000 Colobus Reads")
  
#--> Volcano Plot
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sig.heat_ctl = tt[["table"]][abs(log10(tt[["table"]]$FDR)) > 50,]

ttnew = data.frame(tt[["table"]])

library(stringr)
ttnew = subset(ttnew, str_sub(ttnew$GeneID,1,3) != "LOC")

for (i in 1:length(ttnew$GeneID)){
  if(abs(ttnew$logFC[i])>2 && abs(-log10(ttnew$FDR[i]))>1) {ttnew$color[i] = "blue"} else if (abs(ttnew$logFC[i])>2) {ttnew$color[i] = "black"} else if (abs(-log10(ttnew$FDR[i]))>1) {ttnew$color[i] = "red"} else {ttnew$color[i] = "gray"}}

genesnow = c("UBE2K", "LSM14A", "PP2D1", "APOBEC2", "ACKR1")
for (i in 1:length(ttnew$GeneID)){
  if(ttnew$GeneID[i] %in% genesnow) {ttnew$Name[i] = ttnew$GeneID[i]} else {ttnew$Name[i] = ""}}

bloodgenes = c("RHAG","SPTA1","KLF1","ABCB6","SLC4A1", "SLC2A1", "STEAP3", "JAK2")
for (i in 1:length(ttnew$GeneID)){
  if(ttnew$GeneID[i] %in% bloodgenes) {ttnew$shape[i] = "2"} else {ttnew$shape[i] = "19"}
}
                
p = ggplot(ttnew,
           aes(logFC, -log10(FDR),
               size  = -log10(FDR))) +
  geom_point(aes(logFC, -log10(FDR), color = color, shape = shape), size=2) +
  ggrepel::geom_text_repel(data=ttnew, aes(label = Name), size = 5, color="black", force = 10, min.segment.length = 0) +
  xlab("Fold Change") +
  ylab("-log10(adjusted p-value)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("black", "blue", "grey", "red")) +
  scale_shape_manual(values = c(19, 8))
    
p

#--> Plot genes of interest against parasitemia
# After running DE analysis
# --- Plot genes of interest
  
  plot.gene.expr = function(gene.name) {
    
    cpm.disp.gene_of_interest = cpm.disp[row.names(cpm.disp) == gene.name,]
    ct = data.frame(counts      = cpm.disp.gene_of_interest,
                    parasitemia = parasitemia$parasitemia.proxy.aunin_med)
    
    p = ggplot(ct, aes(parasitemia, counts)) +
      geom_point() +
      geom_smooth(method="lm") +
      xlab("Inferred parasitemia proxy") +
      ylab("Normalized count per million reads") +
      theme_bw()
    ggsave(p, file=paste0("reports/cpm_by_parasitemia.", gene.name, ".", suffix= "pdf"))
  }
  
  lapply(c("ACKR1", "UBE2K", "LSM14A", "APOBEC2", "PP2D1"), plot.gene.expr)

#-> Bar plot with side by side counts
attach(hepato)
stacked_hepato = rbind(data.frame(Sample_name,"count" = colobus_reads, "Species" = "Colobus"),
                      data.frame(Sample_name,"count" = hepato_reads, "Species" = "Hepatocystis"))

ggplot(stacked_hepato, aes(x=Sample_name, y=count, fill=Species)) +
  geom_bar(stat='identity', position='dodge', color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  xlab("Sample Name") +
  ylab("Total Read Count")
