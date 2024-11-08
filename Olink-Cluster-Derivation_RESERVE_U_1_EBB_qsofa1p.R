#Clear R environment
rm(list = ls())

#Import RESERVE-U-1-EBB dataset 
reserve <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve) = reserve$pid
reserve$pid = NULL

#Drop PIDs 361 and 176 - extreme outliers 
reserve <- reserve[row.names(reserve) != "361",]
reserve <- reserve[row.names(reserve) != "176",]

#Restrict to patients with qSOFA score >=1
reserve <- subset(reserve, qsofa1p==1)

#Select and view the log2-transformed Olink data
biomarkers <- reserve[1:242, c(225:408)] 
names(biomarkers)

#Drop biomarkers with <20% of NPX values above each panel’s estimated limits of detection in either of the RESERVE-U cohorts 
library(dplyr)
biomarkers <- dplyr::select(biomarkers, -c("il1alpha",
                                          "il2", 
                                           "il33", 
                                           "il4", 
                                           "il13",
                                           "prcp",
                                           "ltbp2",
                                           "sod1",
                                           "itgam",
                                           "fap",
                                           "mfap5"))  

library(skimr)
skim(biomarkers)

#Reformat column names
colnames(biomarkers) <- base::toupper(colnames(biomarkers))
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
                                        
#Scale and center Olink data
biomscaled <- scale(biomarkers, center=TRUE, scale=TRUE)
biomscaled <- as.data.frame(biomscaled)
skim(biomscaled)

#Hierarchical clustering w/ k-means consolidation
library(FactoMineR)
library (factoextra)

res.hc <- HCPC(biomscaled, 
               nb.clust = -1,      
               consol=TRUE,        
               iter.max = 10,      
               min = 2,            
               max = 12,           
               metric="euclidean", 
               method="ward",      
               graph = TRUE)

#Between group sum of sq before and after k-means consolidation 
res.hc$call$bw.before.consol
res.hc$call$bw.after.consol

#Visualize HC dendrogram 
library(ggsci)
library (factoextra)
dendro <- fviz_dend(res.hc, 
                    cex = 0.7,       
                    palette = "jama",
                    show_labels = FALSE,
                    rect = TRUE, 
                    rect_fill = TRUE,
                    rect_border = "jama",
                    lower_rect = -0.05) + ggtitle(label='') 
dendro

#Display the original input data with cluster assignments
head(res.hc$data.clust, 10)

#Display the standardized means of input variables across clusters 
res.hc$desc.var$quanti

#Label individual patients with their cluster assignments
hcc2 <- res.hc[["data.clust"]]
hcc2

#Merge cluster assignments back into original data 
hcc2 <- data.frame(hcc2$clust)
rownames(hcc2) <- rownames(reserve)
olinkclustassign <- merge(reserve,hcc2,all=T,by='row.names')
rownames(olinkclustassign) = olinkclustassign$Row.names
olinkclustassign$Row.names = NULL

#Drop biomarkers with <20% of NPX values above each panel’s estimated limits of detection in either of the RESERVE-U cohorts
library(dplyr)
olinkclustassign <- dplyr::select(olinkclustassign, -c("il1alpha",
                                           "il2", 
                                           "il33", 
                                           "il4", 
                                           "il13",
                                           "prcp",
                                           "ltbp2",
                                           "sod1",
                                           "itgam",
                                           "fap",
                                           "mfap5"))  

#Subset a dataframe with cluster assignments and original (non-scaled) biomarker concentrations
olinkclustassign_proteins_clusters <- olinkclustassign[1:242, c(225:398)] 
names(olinkclustassign_proteins_clusters)

#Dataframes for each cluster
olinkclustassignclust1 <- subset(olinkclustassign, hcc2.clust==1)
olinkclustassignclust2 <- subset(olinkclustassign, hcc2.clust==2)

#Validate cluster partition using NbClust
library(NbClust)
library(factoextra)

set.seed(12345)
nb <- NbClust(biomscaled, distance = "euclidean", min.nc = 2,
              max.nc = 12, method = "ward.D2", index = "all")

nb$All.index
nb$Best.nc
nb$All.CriticalValues
nb$Best.partition

#Validate robustness of cluster separation across measured proteome via PCA
olinkpcamatrix <- data.matrix(olinkclustassign[1:242, c(225:397)], rownames.force = NA)
head(olinkpcamatrix)
names(olinkpcamatrix)

#Reformat column names
colnames(olinkpcamatrix) <- base::toupper(colnames(olinkpcamatrix))
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL8"] = "IL-8"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL7"] = "IL-7"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL6"] = "IL-6"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL18"] = "IL-18"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL15"] = "IL-15"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL5"] = "IL-5"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL10"] = "IL-10"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IFNGAMMA"] = "IFN-\u03b3"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL12RB1"] = "IL-12RB1"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL12"] = "IL-12"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "IL7R"] = "IL-7R"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "PDGFSUBUNITB"] = "PDGFB"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "LAPTGFBETA1"] = "LAPTGFB1"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "PDL1"] = "PD-L1"
colnames(olinkpcamatrix)[colnames(olinkpcamatrix) == "PDL2"] = "PD-L2"

library("FactoMineR")
library("factoextra")
res.pca <- PCA(olinkpcamatrix, ncp=5, scale.unit = TRUE, graph = TRUE)

print(res.pca)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

#Determine top variables driving separation in PCs 
library(corrplot)
library(tidyverse)
library(tibble)
library(dplyr)
var <- get_pca_var(res.pca)
head(var$cos2, 25)
head(var$contrib, 25)

vars <- fviz_pca_var(res.pca, col.var = "purple4",
                     select.var = list(contrib = 25), labelsize = 6, 
                     xlab = "PC1: 25% expl.var",
                     ylab = "PC2: 9% expl.var",
                     repel = TRUE, geom.var = c("point", "text"))
vars <- vars + theme(text = element_text(size = 16),
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 18))
vars <- vars + ggtitle("")
vars

#First 2 PCs stratified by cluster
res.pca <- PCA(olinkclustassign[1:242, c(225:397)], scale.unit=TRUE, graph = TRUE)
fviz_pca_ind(res.pca)

clust <- fviz_pca_ind(res.pca,
                      geom.ind = "point", # show points only (but not "text")
                      mean.point = FALSE,
                      pointsize = 2.5,
                      xlab = "PC1: 25% expl.var",
                      ylab = "PC2: 9% expl.var",
                      col.ind = as.factor(olinkclustassign$hcc2), 
                      pointshape = 16,
                      addEllipses = FALSE, 
                      legend.title = "")
clust <- clust + theme(text = element_text(size = 16),
                       axis.title = element_text(size = 16),
                       axis.text = element_text(size = 16),
                       legend.text=element_text(size=16))
clust <- clust + ggtitle("") + scale_color_manual(labels = c("USS-1", "USS-2"), values= c("#00468B99", "#ED000099"))
clust

#Compare clinical data and mortality across clusters 
library(gmodels)
library(dplyr)
library(skimr)
library(gtsummary)

#Table - Clinical characteristics by signature
olinkclustassign_table_df <- olinkclustassign[, c("sex", 
                                                  "age",
                                                  "illnessduration_enroll",
                                                  "nightsweats",
                                                  "cough",
                                                  "shortnessofbreath",
                                                  "headache",
                                                  "painwithurination",
                                                  "diarrhea",
                                                  "vomiting",
                                                  "skinrash",
                                                  "tempmax",
                                                  "heartrate3", 
                                                  "resprate3", 
                                                  "sbp3",
                                                  "o2sat3", 
                                                  "ams", 
                                                  "qsofa2p", 
                                                  "qsofa1p", 
                                                  "MEWS_score",
                                                  "UVA_score", 
                                                  "hivrdtresult", 
                                                  "hivstage34", 
                                                  "artprior",
                                                  "malariardtresult", 
                                                  "microtbdx", 
                                                  "urinetblamresult",
                                                  "influenzapcrresult", 
                                                  "ivfluid",
                                                  "supplementaloxygen",
                                                  "bloodtransfusion",
                                                  "antimalarials",
                                                  "antibiotics",
                                                  "antitbdrugs",
                                                  "hospoutcome", 
                                                  "dcperf16andgr", 
                                                  "kps70less",
                                                  "death30d",
                                                  "hcc2.clust")]

theme_gtsummary_compact()

olinkclustassign_table_df %>%
  tbl_summary(
    by = hcc2.clust,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(sex ~ "dichotomous", 
                age ~ "continuous",
                illnessduration_enroll ~ "continuous",
                cough ~ "dichotomous",
                shortnessofbreath ~ "dichotomous",
                headache ~ "dichotomous",
                painwithurination ~ "dichotomous",
                diarrhea ~ "dichotomous",
                vomiting ~ "dichotomous",
                skinrash ~ "dichotomous",
                tempmax ~ "continuous",
                heartrate3 ~ "continuous", 
                resprate3 ~ "continuous", 
                sbp3 ~ "continuous",
                o2sat3 ~ "continuous", 
                ams ~ "dichotomous",   
                qsofa2p ~ "dichotomous",   
                qsofa1p ~ "dichotomous",   
                MEWS_score ~ "continuous",
                UVA_score ~ "continuous", 
                hivrdtresult ~ "dichotomous", 
                hivstage34 ~ "dichotomous", 
                artprior ~ "dichotomous",
                malariardtresult ~ "dichotomous", 
                microtbdx ~ "dichotomous", 
                urinetblamresult ~ "dichotomous",
                influenzapcrresult ~ "dichotomous", 
                ivfluid ~ "dichotomous",
                supplementaloxygen ~ "dichotomous",
                bloodtransfusion ~ "dichotomous",
                antimalarials ~ "dichotomous",
                antibiotics ~ "dichotomous",
                antitbdrugs ~ "dichotomous",
                hospoutcome ~ "categorical", 
                dcperf16andgr ~ "continuous",
                kps70less ~ "dichotomous",
                death30d ~ "dichotomous"),
    digits = list(all_continuous() ~ 0, all_dichotomous() ~ c(0,0,1), all_categorical() ~ c(0,0,1), tempmax ~ c(1, 1)),
    missing_text = "(Missing)",
    label = list(sex ~ "Sex", 
                 age ~ "Age",
                 illnessduration_enroll ~ "Illness duration prior to enrollment, days",
                 cough ~ "Cough",
                 shortnessofbreath ~ "Shortness of breath",
                 headache ~ "Headache",
                 painwithurination ~ "Dysuria",
                 diarrhea ~ "Diarrhea",
                 vomiting ~ "Vomiting",
                 skinrash ~ "Skin rash or discoloration",
                 tempmax ~ "Temperature, \u00b0C",
                 heartrate3 ~ "Heart rate, beats/min", 
                 resprate3 ~ "Respiratory rate, breaths/min", 
                 sbp3 ~ "Systolic blood pressure, mmHg",
                 o2sat3 ~ "Oxygen saturation, %", 
                 ams ~ "Altered mental status",   
                 qsofa2p ~ "qSOFA \u2265 2",   
                 qsofa1p ~ "qSOFA \u2265 1",   
                 MEWS_score ~ "Modified Early Warning Score",
                 UVA_score ~ "Universal Vital Assessment score", 
                 hivrdtresult ~ "Person living with HIV (PLWH)", 
                 hivstage34 ~ "WHO HIV clinical stage 3 or 4", 
                 artprior ~ "Receiving ART prior to hospitalization",
                 malariardtresult ~ "Malaria RDT positive", 
                 microtbdx ~ "Microbiological TB positive", 
                 urinetblamresult ~ "Urine TB-LAM positive if PLWH",
                 influenzapcrresult ~ "Influenza PCR result", 
                 ivfluid ~ "Received intravenous fluids",
                 supplementaloxygen ~ "Received oxygen therapy",
                 bloodtransfusion ~ "Received red blood cell transfusion",
                 antimalarials ~ "Received anti-malarial agent(s)",
                 antibiotics ~ "Received anti-bacterial agent(s)",
                 antitbdrugs ~ "Received anti-TB agent(s)",
                 hospoutcome ~ "Hospital outcome", 
                 dcperf16andgr ~ "KPS at alive discharge or transfer", 
                 kps70less ~ "KPS ≤ 70 at alive discharge or transfer",
                 death30d ~ "Dead at 30 days post-discharge")) %>% add_overall() %>% add_n() %>% as_gt() 

#Compare 30d mortality by signature in PLWH vs. without HIV
library(gmodels)

with(olinkclustassign[olinkclustassign$hivrdtresult==1, ], CrossTable(death30d, hcc2.clust))
with(olinkclustassign[olinkclustassign$hivrdtresult==0, ], CrossTable(death30d, hcc2.clust))

#Patient characteristics by signature
library(gmodels)
library(catstats)
library(boot)

#Demographics and pre-enrollment illness duration 
#Age
skim(subset(olinkclustassign$age, olinkclustassign$hcc2.clust==1))
skim(subset(olinkclustassign$age, olinkclustassign$hcc2.clust==2))

group_by(olinkclustassign, hcc2.clust) %>%
  summarise(
    count = n(),
    median = median(age, na.rm = TRUE),
    IQR = IQR(age), na.rm = TRUE)

wilcox.test(age ~ hcc2.clust, data = olinkclustassign, exact = FALSE)

boot_diff_median_age <- function(olinkclustassign, i){
  diff(tapply(olinkclustassign$age[i], olinkclustassign$hcc2.clust[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_age <- boot(olinkclustassign, statistic = boot_diff_median_age, R = 10000)

median(b_age$t)
t(quantile(b_age$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Sex
CrossTable(olinkclustassign$sex, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(sex ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Illness duration before enrollment
skim(subset(olinkclustassign$illnessduration_enroll, olinkclustassign$hcc2.clust==1))
skim(subset(olinkclustassign$illnessduration_enroll, olinkclustassign$hcc2.clust==2))

group_by(olinkclustassign, hcc2.clust) %>%
  summarise(
    count = n(),
    median = median(illnessduration_enroll, na.rm = TRUE),
    IQR = IQR(illnessduration_enroll), na.rm = TRUE)

wilcox.test(illnessduration_enroll ~ hcc2.clust, data = olinkclustassign, exact = FALSE)

boot_diff_median_illnessduration_enroll <- function(olinkclustassign, i){
  diff(tapply(olinkclustassign$illnessduration_enroll[i], olinkclustassign$hcc2.clust[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_illnessduration_enroll <- boot(olinkclustassign, statistic = boot_diff_median_illnessduration_enroll, R = 10000)

median(b_illnessduration_enroll$t)
t(quantile(b_illnessduration_enroll$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Microbiology 
#HIV
CrossTable(olinkclustassign$hivrdtresult, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(hivrdtresult ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#HIV stage 3/4
CrossTable(olinkclustassign$hivstage34, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(hivstage34 ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#TB
CrossTable(olinkclustassign$microtbdx, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(microtbdx ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#TB LAM in PLWH
#Subset to HIV only for TB-LAM proportion
olinkclustassign_hiv <- subset(olinkclustassign, hivrdtresult==1)

CrossTable(olinkclustassign_hiv$urinetblamresult, olinkclustassign_hiv$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(urinetblamresult ~ hcc2.clust,
                            data = olinkclustassign_hiv,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Malaria
CrossTable(olinkclustassign$malariardtresult, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(malariardtresult ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Influenza 
CrossTable(olinkclustassign$influenzapcrresult, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(influenzapcrresult ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#MEWS
library(boot)
skim(subset(olinkclustassign$MEWS_score, olinkclustassign$hcc2.clust==1))
skim(subset(olinkclustassign$MEWS_score, olinkclustassign$hcc2.clust==2))

group_by(olinkclustassign, hcc2.clust) %>%
  summarise(
    count = n(),
    median = median(MEWS_score, na.rm = TRUE),
    IQR = IQR(MEWS_score), na.rm = TRUE)

wilcox.test(MEWS_score ~ hcc2.clust, data = olinkclustassign, exact = FALSE)

boot_diff_median_mews <- function(olinkclustassign, i){
  diff(tapply(olinkclustassign$MEWS_score[i], olinkclustassign$hcc2.clust[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_mews <- boot(olinkclustassign, statistic = boot_diff_median_mews, R = 10000)

median(b_mews$t)
t(quantile(b_mews$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#UVA
skim(subset(olinkclustassign$UVA_score, olinkclustassign$hcc2.clust==1))
skim(subset(olinkclustassign$UVA_score, olinkclustassign$hcc2.clust==2))

group_by(olinkclustassign, hcc2.clust) %>%
  summarise(
    count = n(),
    median = median(UVA_score, na.rm = TRUE),
    IQR = IQR(UVA_score), na.rm = TRUE)

wilcox.test(UVA_score ~ hcc2.clust, data = olinkclustassign, exact = FALSE)

boot_diff_median_uva <- function(olinkclustassign, i){
  diff(tapply(olinkclustassign$UVA_score[i], olinkclustassign$hcc2.clust[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_uva <- boot(olinkclustassign, statistic = boot_diff_median_uva, R = 10000)

median(b_uva$t)
t(quantile(b_uva$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#qSOFA >=2
CrossTable(olinkclustassign$qsofa2p, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(qsofa2p ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Death in hospital
#Create variable
olinkclustassign[olinkclustassign$hospoutcome == "Death", "deathhosp"] <-1
olinkclustassign <- mutate_at(olinkclustassign, c("deathhosp"), ~replace(., is.na(.),0))

CrossTable(olinkclustassign$deathhosp, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(deathhosp ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#KPS at discharge/transfer 
skim(subset(olinkclustassign$dcperf16andgr, olinkclustassign$hcc2.clust==1))
skim(subset(olinkclustassign$dcperf16andgr, olinkclustassign$hcc2.clust==2))

group_by(olinkclustassign, hcc2.clust) %>%
  summarise(
    count = n(),
    median = median(dcperf16andgr, na.rm = TRUE),
    IQR = IQR(dcperf16andgr), na.rm = TRUE)

wilcox.test(dcperf16andgr ~ hcc2.clust, data = olinkclustassign, exact = FALSE)

boot_diff_median_dcperf16andgr <- function(olinkclustassign, i){
  diff(tapply(olinkclustassign$dcperf16andgr[i], olinkclustassign$hcc2.clust[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_dcperf16andgr <- boot(olinkclustassign, statistic = boot_diff_median_dcperf16andgr, R = 10000)

median(b_dcperf16andgr$t)
t(quantile(b_dcperf16andgr$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Death at 30d
CrossTable(olinkclustassign$death30d, olinkclustassign$hcc2.clust, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(death30d ~ hcc2.clust,
                            data = olinkclustassign,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Logistic regression for association between signature and in-hospital/30d mortality risk 
#In-hospital 
#Univariable
signaturehosp_uv <- glm(deathhosp ~ hcc2.clust, data = olinkclustassign, family = "binomial")
summary(signaturehosp_uv)
exp(cbind(OR = coef(signaturehosp_uv), confint(signaturehosp_uv)))

#Multivariable 
signaturehosp_mv <- glm(deathhosp ~ hcc2.clust + age + sex + illnessduration_enroll + hivrdtresult, data = olinkclustassign, family = "binomial")
summary(signaturehosp_mv)
exp(cbind(OR = coef(signaturehosp_mv), confint(signaturehosp_mv)))

#30d 
#Univariable
signature30d_uv <- glm(death30d ~ hcc2.clust, data = olinkclustassign, family = "binomial")
summary(signature30d_uv)
exp(cbind(OR = coef(signature30d_uv), confint(signature30d_uv)))

#Interaction between signature and HIV for 30d mortality risk 
signature30d_interact_hiv <- glm(death30d ~ hcc2.clust*hivrdtresult, data = olinkclustassign, family = "binomial")
summary(signature30d_interact_hiv)

#Multivariable 
signature30d_mv <- glm(death30d ~ hcc2.clust + age + sex + illnessduration_enroll + hivrdtresult, data = olinkclustassign, family = "binomial")
summary(signature30d_mv)
exp(cbind(OR = coef(signature30d_mv), confint(signature30d_mv)))

#Generate protein classifier to discriminate signatures 
#Random forest
library(tidyverse)
library(randomForest)
library(RColorBrewer)
library(caret)
library(mlr)
library(xgboost)

#Recode signature assignment variable
olinkclustassign_proteins_clusters$signature2 <- NA
olinkclustassign_proteins_clusters[olinkclustassign_proteins_clusters$hcc2.clust==2, "signature2"] <- 1
olinkclustassign_proteins_clusters[olinkclustassign_proteins_clusters$hcc2.clust==1, "signature2"] <- 0

#Extract outcome and biomarkers
data_for_model <- data.frame(olinkclustassign_proteins_clusters$signature2, olinkclustassign_proteins_clusters[,c(1:173)] )
colnames(data_for_model)[1] <- "signature2"

#Convert the response variable to a factor
data_for_model$signature2 = as.factor(data_for_model$signature2)

#Set
response_var <- data_for_model$signature2
predictor_vars <- data.matrix(data_for_model[,c(2:ncol(data_for_model))])
length(response_var)

#Random Forest
set.seed(123)
oob.values <- vector(length=10)

for(i in 1:10) {
  temp.model <- randomForest(signature2 ~ . , data=data_for_model, mtry=i, maxnodes=10, ntree=1000)
  oob.values[i] <- temp.model$err.rate[nrow(temp.model$err.rate),1]
}

print(oob.values)

min(oob.values)

which(oob.values == min(oob.values))

rf <- randomForest(signature2 ~ . , data=data_for_model, maxnodes=10, ntree=1000, 
                   mtry=which(oob.values == min(oob.values)))
print(rf)

oob.error.data <- data.frame(
  Trees=rep(1:nrow(rf$err.rate), times=3),
  Type=rep(c("OOB", "0", "1"), each=nrow(rf$err.rate)),
  Error=c(rf$err.rate[,"OOB"],
          rf$err.rate[,"0"],
          rf$err.rate[,"1"]))

#Plot
ggplot(data=oob.error.data, aes(x= Trees, y=Error)) +
  geom_line(aes(color=Type))

rfimportance <- data.frame(imp = as.numeric(rf$importance),
                           feature = row.names(rf$importance)) %>% top_n(20, rf$importance) %>%
  arrange(-imp) %>% ggplot(aes(x=imp, y=reorder(feature, imp))) + 
  geom_col(fill = "aquamarine3", show.legend = FALSE) + 
  xlab("Random Forest Importance (Gini Impurity)") +
  ylab("Protein")
rfimportance

#RF model
#Model
signature5protein <- glm(signature2 ~ cd40 + 
                         pgf + 
                         timp1 + 
                         cd46 + 
                         cst3,
                       data = olinkclustassign_proteins_clusters, family = "binomial")
summary(signature5protein)
exp(signature5protein$coefficients)
exp(confint.default(signature5protein))

#ROC
library(pROC)

predicted <- predict(signature5protein, olinkclustassign_proteins_clusters, type="response")
auc(olinkclustassign_proteins_clusters$signature2, predicted)


rocobj <- plot.roc(olinkclustassign_proteins_clusters$signature2, signature5protein$fitted.values,
                   percent=FALSE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE,
                   print.auc.x = 0.6, 
                   print.auc.y = 0.05,
                   legacy.axes = TRUE,
                   print.auc.pattern = "AUROC %.3f (%.3f-%.3f)")           # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj, boot.n=10000,                         # CI of sensitivity
               specificities = seq(0, 1, 0.05)) # over a select set of specificities
plot(ciobj, type = "shape", col = "#8da0cb")     # plot as a blue shape

#Now, compute probabilities and Youden index for RF model using 100x repeated 10-fold CV
library(caret)

#Set CV parameters and seed
train.control <- trainControl(method="repeatedcv", number=10, repeats = 100, 
                              savePredictions = "all", classProbs = TRUE)                

set.seed(12345)

#Recode signature assignment to factor
olinkclustassign_proteins_clusters[olinkclustassign_proteins_clusters$signature2==1, "signature2_factor"] <- "yes"
olinkclustassign_proteins_clusters[olinkclustassign_proteins_clusters$signature2==0, "signature2_factor"] <- "no"
olinkclustassign_proteins_clusters$signature2_factor <- as.factor(olinkclustassign_proteins_clusters$signature2_factor)
levels(olinkclustassign_proteins_clusters$signature2_factor)
olinkclustassign_proteins_clusters$signature2_factor <- relevel(olinkclustassign_proteins_clusters$signature2_factor, "yes")
levels(olinkclustassign_proteins_clusters$signature2_factor)

#Logistic model
model_rf <- caret::train(signature2_factor ~ cd40 + 
                    pgf + 
                    timp1 + 
                    cd46 + 
                    cst3,
                  data = olinkclustassign_proteins_clusters, method = "glm",
                   family = "binomial", trControl = train.control)

#Model accuracy
print(model_rf)

#Extract probabilities 
model_rf_predict <- predict(model_rf, type = "prob")

#Determine optimal probability threshold
probs <- seq(0.3, 0.7, by = 0.10)

ths <- thresholder(model_rf,
                   threshold = probs,
                   final = TRUE,
                   statistics = "all")

head(ths)

ths %>%
  mutate(prob = probs) %>%
  filter(J == max(J)) %>%
  pull(prob) -> thresh_prob
thresh_prob

#Cluster assignment mapped to clinical PCA
clinicalpcamatrix <- as.data.frame(olinkclustassign[1:242, c(60,66,69,72,78,83,90,98,106,156,157,184)], rownames.force = NA)
head(clinicalpcamatrix)

clinicalpcamatrix <- clinicalpcamatrix %>% mutate_all(~replace(., is.na(.), 0))
clinicalpcamatrix <- clinicalpcamatrix %>% mutate_if(is.character, as.numeric)
clinicalpcamatrix <- clinicalpcamatrix %>% mutate_if(is.factor, as.numeric)

head(clinicalpcamatrix)

#Reformat column names
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "tempmax"] = "Temp."
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "heartrate3"] = "HR"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "resprate3"] = "RR"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "sbp3"] = "SBP"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "o2sat3"] = "SpO2"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "avpu3"] = "AVPU"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "qsofa_score"] = "qSOFA"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "MEWS_score"] = "MEWS"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "UVA_score"] = "UVA"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "malariardtresult"] = "Malaria"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "hivrdtresult"] = "HIV"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "microtbdx"] = "TB"

library("FactoMineR")
library("factoextra")
res.pca_clinical <- PCA(clinicalpcamatrix, ncp=5, scale.unit = TRUE, graph = TRUE)

print(res.pca_clinical)

eigenvalues <- res.pca_clinical$eig
head(eigenvalues[, 1:2])
res.pca_clinical$var$contrib

#First 2 clinical PCs stratified by cluster
clust_clinical <- fviz_pca_ind(res.pca_clinical,
                      geom.ind = "point", # show points only (but not "text")
                      mean.point = FALSE,
                      pointsize = 2.5,
                      xlab = "PC1: 31% expl.var",
                      ylab = "PC2: 12% expl.var",
                      col.ind = as.factor(olinkclustassign$hcc2), 
                      alpha.ind = "cos2", 
                      pointshape = 16,
                      addEllipses = FALSE, 
                      legend.title = "")
clust_clinical <- clust_clinical + theme(text = element_text(size = 16),
                       axis.title = element_text(size = 16), 
                       axis.text = element_text(size = 16),
                       legend.text=element_text(size=16))
clust_clinical <- clust_clinical + ggtitle("") + scale_color_manual(labels = c("USS-1", "USS-2"), values= c("#FF725C", "#374e55"))
clust_clinical

#Variable Plot
library(corrplot)
library(tidyverse)
library(tibble)
library(dplyr)
vars_clinical <- get_pca_var(res.pca_clinical)
head(vars_clinical$cos2, 25)
head(vars_clinical$contrib, 25)

vars_clinical_plot <- fviz_pca_var(res.pca_clinical, col.var="navy",
                     labelsize = 6, 
                     repel = TRUE,
                     xlab = "PC1: 31% expl.var",
                     ylab = "PC2: 12% expl.var")
vars_clinical_plot <- vars_clinical_plot + theme(text = element_text(size = 16),
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 18))
vars_clinical_plot <- vars_clinical_plot + ggtitle("")
vars_clinical_plot

