#Clear R environment
rm(list = ls())

#Load necessary libraries
library(limma)
library(ggplot2)
library(dplyr)

#Specify cohort
cohort="RESERVE-U-2-TOR" 

#Read data file
biomarkers <- read.csv(file.choose(), header = TRUE)

#Reformat column names
colnames(biomarkers) = c('pid', base::toupper(colnames(biomarkers[,-c(1,ncol(biomarkers))])), 'hcc2.clust')
colnames(biomarkers)[colnames(biomarkers) == "IL8"] = "IL-8"
colnames(biomarkers)[colnames(biomarkers) == "IL7"] = "IL-7"
colnames(biomarkers)[colnames(biomarkers) == "IL6"] = "IL-6"
colnames(biomarkers)[colnames(biomarkers) == "IL18"] = "IL-18"
colnames(biomarkers)[colnames(biomarkers) == "IL15"] = "IL-15"
colnames(biomarkers)[colnames(biomarkers) == "IL5"] = "IL-5"
colnames(biomarkers)[colnames(biomarkers) == "IL10"] = "IL-10"
colnames(biomarkers)[colnames(biomarkers) == "IFNGAMMA"] = "IFN-\u03b3"
colnames(biomarkers)[colnames(biomarkers) == "IL12RB1"] = "IL-12RB1"
colnames(biomarkers)[colnames(biomarkers) == "IL12"] = "IL-12"
colnames(biomarkers)[colnames(biomarkers) == "IL7R"] = "IL-7R"
colnames(biomarkers)[colnames(biomarkers) == "PDGFSUBUNITB"] = "PDGFB"
colnames(biomarkers)[colnames(biomarkers) == "LAPTGFBETA1"] = "LAPTGFB1"
colnames(biomarkers)[colnames(biomarkers) == "PDL1"] = "PD-L1"
colnames(biomarkers)[colnames(biomarkers) == "PDL2"] = "PD-L2"

data = biomarkers

#Further format column names
data$pid = paste0('PID', data$pid, '')
data$hcc2.clust = paste0('C', data$hcc2.clust, '')
#Transpose the data so that proteins are rows and samples are columns
expression_data <- t(data[, -c(1,ncol(data))])  # remove the 'clust'/'PID' columns and transpose
colnames(expression_data) = data$pid

group_assignments = factor(data$hcc2.clust, levels = c("C1", "C2"))
  
#Create the design matrix for the linear model
design <- model.matrix(~ 0 + group_assignments)
colnames(design) <- levels(group_assignments)

#Run limma
fit <- lmFit(expression_data, design)

#Apply empirical Bayes statistics
fit <- eBayes(fit)

#Create a contrast matrix for the comparison of C2 vs C1
contrast.matrix <- makeContrasts(C2vsC1 = C2 - C1, levels = design)

#Fit the contrast model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Get the table of differentially expressed proteins
results_tor <- topTable(fit2, coef="C2vsC1", adjust.method="fdr", sort.by="P", number=Inf)

#Volcano plot 
library(EnhancedVolcano)
#Create a data frame for the volcano plot
results_tor$Prot <- toupper(rownames(results_tor))  # Create a column for gene names if not already present
names(results_tor$logFC) <- toupper(results_tor$Prot)
names(results_tor$P.Value) <- toupper(results_tor$Prot)

results_tor$selectedLab <- ifelse(results_tor$adj.P.Val < 0.05 & abs(results_tor$logFC) > 1, as.character(results_tor$Prot), NA)
EnhancedVolcano(
  results_tor,
  lab = NA,
  x = 'logFC',
  y = 'adj.P.Val',
  xlim = c(-2.8,2.8),
  title = paste('Volcano Plot using FDR correction -',cohort, sep = " "),
  pCutoff = 0.05,
  FCcutoff = 1,  
  pointSize = 4,
  labSize = 6,
  legendLabSize = 15.0,
  axisLabSize = 25.0,
  max.overlaps = 150,
  drawConnectors = TRUE
)+ geom_text_repel(
  data = subset(results_tor, !is.na(selectedLab)),  # Only use rows where labels are not NA
  aes(label = selectedLab, x = logFC, y = -log10(adj.P.Val)),  # Ensure 'label' aesthetic is correctly mapped
  size = 6,  # Ensure size matches labSize in EnhancedVolcano
  box.padding = unit(0.3, "lines"),  # Adjust padding to avoid overlapping
  point.padding = unit(0.2, "lines"),  # Increase point padding
  max.overlaps = 150,  # Allow more overlaps
  #  segment.color = NA
)

