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
  
plasfit <- lm(hep_per_hunthoucol ~ col_hunthou, data = hepato)
summary(plasfit)
