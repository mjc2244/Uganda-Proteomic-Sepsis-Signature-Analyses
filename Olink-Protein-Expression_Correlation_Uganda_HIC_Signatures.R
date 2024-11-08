#Clear R environment
rm(list = ls())

#Import data
reserve <- read.csv(file.choose(), header=TRUE)

#Proteins to fixed row identifiers
rownames(reserve) = reserve$protein
reserve$protein = NULL

#Correlation between protein expression
library(stats)
library(RVAideMemoire)

#Uganda Signature vs. Hyper/Hypo-inflammatory
cor.test(reserve$logFC_ug, reserve$logFC_hyper, method = "spearman")
spearman.ci(reserve$logFC_ug, reserve$logFC_hyper, nrep = 10000, conf.level = 0.95)

#Uganda Signature vs. Reactive/Uninflamed
cor.test(reserve$logFC_ug, reserve$logFC_reactive, method = "spearman")
spearman.ci(reserve$logFC_ug, reserve$logFC_reactive, nrep = 10000, conf.level = 0.95)

#Uganda Signature vs. SRS
cor.test(reserve$logFC_ug, reserve$logFC_srs, method = "spearman")
spearman.ci(reserve$logFC_ug, reserve$logFC_srs, nrep = 10000, conf.level = 0.95)

#Plots 
library(ggplot2)
#Uganda Signature vs. Hyper/Hypo-inflammatory
ug_vs_hyper_hypo <- ggplot(reserve, aes(x = logFC_ug, y = logFC_hyper)) + 
  geom_point(size=3, color="#0073C2FF") + 
  theme_bw() + 
  ggtitle("USS vs. Hyper/Hypo-inflammatory") +
  ylab("Hyper vs. Hypo Log2FC") + xlab("USS-2 vs. USS-1 Log2FC") +
  annotate("text", x=0, y=2, label= "\u03c1 = 0.89 (0.84-0.92)", size = 5) +
  theme(text=element_text(size=13), axis.text=element_text(size=13))
ug_vs_hyper_hypo

#Uganda Signature vs. Reactive/Uninflamed
ug_vs_reactive_uninflamed <- ggplot(reserve, aes(x = logFC_ug, y = logFC_reactive)) + 
  geom_point(size=3, color="#EFC000FF") + 
  theme_bw() + 
  ggtitle("USS vs. Reactive/Uninflamed") +
  ylab("Reactive vs. Uninflamed Log2FC") + xlab("USS-2 vs. USS-1 Log2FC") +
  annotate("text", x=0, y=3, label= "\u03c1 = 0.84 (0.78-0.88)", size = 5) +
  theme(text=element_text(size=13), axis.text=element_text(size=13))
ug_vs_reactive_uninflamed

#Uganda Signature vs. SRS
ug_vs_srs <- ggplot(reserve, aes(x = logFC_ug, y = logFC_srs)) + 
  geom_point(size=3, color="#868686FF") + 
  theme_bw() + 
  ggtitle("USS vs. SRS") +
  ylab("SRS1 vs. SRS2/3 Log2FC") + xlab("USS-2 vs. USS-1 Log2FC") +
  annotate("text", x=0, y=1.5, label= "\u03c1 = 0.40 (0.24-0.54)", size = 5) +
  theme(text=element_text(size=13), axis.text=element_text(size=13))
ug_vs_srs

#Combine plots
library(ggpubr)

top_row = ggarrange(ug_vs_hyper_hypo, ug_vs_reactive_uninflamed, ncol = 2, labels = c("", ""))
bottom_row = ggarrange(NULL, ug_vs_srs, NULL, ncol = 3, labels = c("", "", ""), widths = c(1,2,1))
final_plot = ggarrange(top_row, bottom_row, ncol = 1)
final_plot