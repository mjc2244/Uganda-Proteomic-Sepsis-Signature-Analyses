#Import RESERVE-U-2-TOR dataset
reserve_tor <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve_tor) = reserve_tor$pid
reserve_tor$pid = NULL

#Restrict to patients with qSOFA score >=1
reserve_tor <- subset(reserve_tor, qsofa_1p_avpu==1)

#Assignment to signature based on RESERVE-U-1-EBB classifier model 
options(scipen = 999)
reserve_tor$signature2_probs <- predict(signature5protein, reserve_tor, type="response")

#Assign to signature based on probability cutoff >=0.50 (per maximized Youden index in RESERVE-U-1-EBB)
reserve_tor$signature_predicted <- ifelse(reserve_tor$signature2_probs >= 0.5, 2, 1)
reserve_tor$signature_predicted <- as.factor(reserve_tor$signature_predicted)

#Drop biomarkers with <20% of NPX values above each panel’s estimated limits of detection in either of the RESERVE-U cohorts
library(dplyr)
reserve_tor <- dplyr::select(reserve_tor, -c("il1alpha",
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

#Validate robustness of signature separation across measured proteome via PCA
olinkpcamatrix <- data.matrix(reserve_tor[1:253, c(359:531)], rownames.force = NA)
head(olinkpcamatrix)

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
library("corrplot")
var <- get_pca_var(res.pca)
head(var$cos2, 25)
head(var$contrib, 25)

vars <- fviz_pca_var(res.pca, col.var = "purple4",
             select.var = list(contrib = 25), labelsize = 6, 
             xlab = "PC1: 33% expl.var",
             ylab = "PC2: 10% expl.var",
             repel = TRUE, geom.var = c("point", "text"))
vars <- vars + theme(text = element_text(size = 16),
                     axis.title = element_text(size = 16),
                     axis.text = element_text(size = 16))
vars <- vars + ggtitle("")
vars

#First 2 PCs stratified by cluster
res.pca <- PCA(reserve_tor[1:253, c(359:531)], scale.unit=TRUE, graph = TRUE)
fviz_pca_ind(res.pca)

clust <- fviz_pca_ind(res.pca,
                      geom.ind = "point", # show points only (but not "text")
                      mean.point = FALSE,
                      pointsize = 2.5,
                      xlab = "PC1: 33% expl.var",
                      ylab = "PC2: 10% expl.var",
                      col.ind = as.factor(reserve_tor$signature_predicted), 
                      pointshape = 16,
                      addEllipses = FALSE, 
                      legend.title = "")
clust <- clust + theme(text = element_text(size = 16),
                       axis.title = element_text(size = 16),
                       axis.text = element_text(size = 16),
                       legend.text=element_text(size=16))
clust <- clust + ggtitle("") + scale_color_manual(labels = c("USS-1", "USS-2"), values= c("#00468B99", "#ED000099"))
clust

#Subset a dataframe with signature assignments and original (non-scaled) biomarker concentrations
olinkclustassign_proteins_clusters_tor <- reserve_tor[1:253, c(359:531,533)] 
names(olinkclustassign_proteins_clusters_tor)

#Compare clinical data, mortality across clusters 
library(gmodels)
library(dplyr)
library(skimr)
library(gtsummary)
library(boot)

#Recode variable for discharge KPS only in patients discharged alive or transferred (exclude deaths)
reserve_tor[reserve_tor$hosp_outcome == 2, c("dc_perf_16_and_above")] <- NA

reserve_tor_table_df <- reserve_tor[, c(  
            "age",
            "gender",
            "muac_imputed",
            "illnessduration_to_enroll",
            "admit_perf_16_and_above",
            "assist_oob_walk",
            "cough",
            "sob",
            "headache",
            "dysuria",
            "diarrhea",
            "vomiting",
            "rash",
            "tempmax_imputed",
            "heart_rate_3", 
            "resp_rate_3", 
            "sbp_3",
            "o2_sat_3",
            "ams_avpu", 
            "lactate_result_imputed",
            "lactate4plus",
            "qsofa_2p_avpu", 
            "mews_score",
            "uva_score_avpu", 
            "hiv_rdt_result", 
            "hiv_stage_3_4_all", 
            "art_prior_all", 
            "cd4_imputed",
            "suppressedviralload_imputed",
            "malaria_rdt_result", 
            "microtbdx", 
            "urine_lam_result",
            "cryto_ag_result",
            "flu_pcr_uvri_result", 
            "covid_pcr_uvri_result",
            "supp_o2_hosp",
            "ivf_hosp",
            "antimalarials_hosp",
            "abx_hosp",
            "anti_tb_hosp",                     
            "blood_transfusion_hosp",
            "hosp_outcome", 
            "dc_perf_16_and_above",
            "patient_status_30d", 
            "death30d",
            "kps_30d_16plus",
            "patient_status_60d",
            "death60d",
            "kps_60d_16plus",
            "signature_predicted")]

theme_gtsummary_compact()

reserve_tor_table_df %>%
  tbl_summary(
    by = signature_predicted,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(age ~ "continuous",
                gender ~ "categorical",
                muac_imputed ~ "continuous",
                illnessduration_to_enroll ~ "continuous",
                cough ~ "categorical",
                sob ~ "categorical",
                headache ~ "categorical",
                dysuria ~ "categorical",
                diarrhea ~ "categorical",
                vomiting ~ "categorical",
                rash ~ "categorical",
                admit_perf_16_and_above ~ "continuous",
                assist_oob_walk ~ "dichotomous",
                lactate_result_imputed~ "continuous",
                lactate4plus ~ "dichotomous",
                qsofa_2p_avpu ~ "dichotomous", 
                mews_score ~ "continuous",
                uva_score_avpu ~ "continuous", 
                hiv_rdt_result ~ "dichotomous",
                art_prior_all ~ "dichotomous", 
                cd4_imputed ~ "continuous", 
                suppressedviralload_imputed ~ "dichotomous",
                malaria_rdt_result ~ "dichotomous",
                microtbdx ~ "dichotomous",
                urine_lam_result ~ "dichotomous",
                cryto_ag_result ~ "dichotomous",
                flu_pcr_uvri_result ~ "dichotomous",
                covid_pcr_uvri_result ~ "dichotomous",
                supp_o2_hosp ~ "dichotomous",
                ivf_hosp ~ "dichotomous",
                antimalarials_hosp ~ "dichotomous",
                abx_hosp ~ "dichotomous",
                anti_tb_hosp ~ "dichotomous",                    
                blood_transfusion_hosp ~ "dichotomous",
                hosp_outcome ~ "categorical",
                dc_perf_16_and_above ~ "continuous",
                patient_status_30d ~ "categorical", 
                death30d ~ "dichotomous",
                kps_30d_16plus ~ "continuous",
                patient_status_60d ~ "categorical", 
                death60d ~ "dichotomous",
                kps_60d_16plus ~ "continuous"),
    digits = list(all_continuous() ~ 0, all_dichotomous() ~ c(0,0,1), all_categorical() ~ c(0,0,1), tempmax_imputed ~ c(1, 1), lactate_result_imputed ~ c(1, 1)),
    missing_text = "(Missing)",
    label = list(age ~ "Age",
                 gender ~ "Sex",
                 muac_imputed ~ "Mid-upper arm circumference, cm",
                 illnessduration_to_enroll ~ "Illness duration prior to enrollment, days",
                 admit_perf_16_and_above ~ "KPS at hospital admission",
                 assist_oob_walk ~ "Unable to ambulate without assistance",
                 cough ~ "Cough", 
                 sob ~ "Shortness of breath",
                 headache ~ "Headache",
                 dysuria ~ "Dysuria",
                 diarrhea ~ "Diarrhea",
                 vomiting ~ "Vomiting",
                 rash ~ "Skin rash or discoloration",
                 tempmax_imputed ~ "Temperature, \u00b0C",
                 heart_rate_3 ~ "Heart rate, beats/min", 
                 resp_rate_3 ~ "Respiratory rate, breaths/min", 
                 sbp_3 ~ "Systolic blood pressure, mmHg",
                 o2_sat_3 ~ "Oxygen saturation, %",
                 ams_avpu ~ "Altered mental status", 
                 lactate_result_imputed~ "Whole-blood lactate, mmol/L",
                 lactate4plus ~ "Whole-blood lactate \u2265 4 mmol/L",
                 qsofa_2p_avpu ~ "qSOFA \u2265 2", 
                 mews_score ~ "Modified Early Warning Score",
                 uva_score_avpu ~ "Universal Vital Assessment score", 
                 hiv_rdt_result ~ "PLWH",
                 hiv_stage_3_4_all ~ "WHO HIV clinical stage 3 or 4", 
                 art_prior_all ~ "Receiving ART prior to hospitalization", 
                 cd4_imputed ~ "CD4 count, cells/mm3",
                 suppressedviralload_imputed ~ "HIV-1 viral suppression",
                 malaria_rdt_result ~ "Malaria RDT positive",
                 microtbdx ~ "Microbiological TB positive",
                 urine_lam_result ~ "Urine TB-LAM positive if PLWH",
                 cryto_ag_result ~ "Cryptococcal ag. positive if PLWH",
                 flu_pcr_uvri_result ~ "Influenza PCR positive",
                 covid_pcr_uvri_result ~ "SARS-CoV-2 PCR positive",
                 supp_o2_hosp ~ "Received oxygen therapy",
                 ivf_hosp ~ "Received intravenous fluids",
                 antimalarials_hosp ~ "Received anti-malarial agent(s)",
                 abx_hosp ~ "Received anti-bacterial agent(s)",
                 anti_tb_hosp ~ "Received anti-TB agent(s)",                    
                 blood_transfusion_hosp ~ "Received red blood cell transfusion",
                 hosp_outcome ~ "Hospital outcome",
                 dc_perf_16_and_above ~ "KPS at discharge or transfer",
                 patient_status_30d ~ "Vital status at 30 days post-discharge", 
                 death30d ~ "Death at 30 days post-discharge",
                 kps_30d_16plus ~ "KPS at 30 days post-discharge",
                 patient_status_60d ~ "Vital status at 60 days post-discharge", 
                 death60d ~ "Death at 60 days post-discharge",
                 kps_60d_16plus ~ "KPS at 60 days post-discharge")) %>% add_overall() %>% add_p() %>% add_n() %>% as_gt() 

#Compare 30d mortality by signature in PLWH vs. without HIV
library(gmodels)

with(reserve_tor[reserve_tor$hiv_rdt_result==1, ], CrossTable(death30d, signature_predicted))
with(reserve_tor[reserve_tor$hiv_rdt_result==0, ], CrossTable(death30d, signature_predicted))

#Compare 60d mortality by signature in PLWH vs. without HIV
library(gmodels)

with(reserve_tor[reserve_tor$hiv_rdt_result==1, ], CrossTable(death60d, signature_predicted))
with(reserve_tor[reserve_tor$hiv_rdt_result==0, ], CrossTable(death60d, signature_predicted))

#Patient characteristics by signature
library(gmodels)
library(catstats)
library(boot)

#Demographics and pre-enrollment illness duration 
#Age
skim(subset(reserve_tor$age, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$age, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(age, na.rm = TRUE),
    IQR = IQR(age), na.rm = TRUE)

wilcox.test(age ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_age <- function(reserve_tor, i){
  diff(tapply(reserve_tor$age[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_age <- boot(reserve_tor, statistic = boot_diff_median_age, R = 10000)

median(b_age$t)
t(quantile(b_age$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Sex
CrossTable(reserve_tor$gender, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(gender ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "2",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#MUAC
skim(subset(reserve_tor$muac_imputed, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$muac_imputed, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(muac_imputed, na.rm = TRUE),
    IQR = IQR(muac_imputed), na.rm = TRUE)

wilcox.test(muac_imputed ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_muac_imputed <- function(reserve_tor, i){
  diff(tapply(reserve_tor$muac_imputed[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_muac_imputed <- boot(reserve_tor, statistic = boot_diff_median_muac_imputed, R = 10000)

median(b_muac_imputed$t)
t(quantile(b_muac_imputed$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)) 

#Illness duration before enrollment
skim(subset(reserve_tor$illnessduration_to_enroll, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$illnessduration_to_enroll, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(illnessduration_to_enroll, na.rm = TRUE),
    IQR = IQR(illnessduration_to_enroll), na.rm = TRUE)

wilcox.test(illnessduration_to_enroll ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_illnessduration_to_enroll <- function(reserve_tor, i){
  diff(tapply(reserve_tor$illnessduration_to_enroll[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_illnessduration_to_enroll <- boot(reserve_tor, statistic = boot_diff_median_illnessduration_to_enroll, R = 10000)

median(b_illnessduration_to_enroll$t)
t(quantile(b_illnessduration_to_enroll$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#KPS at admission
skim(subset(reserve_tor$admit_perf_16_and_above, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$admit_perf_16_and_above, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(admit_perf_16_and_above, na.rm = TRUE),
    IQR = IQR(admit_perf_16_and_above), na.rm = TRUE)

wilcox.test(admit_perf_16_and_above ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_admit_perf_16_and_above <- function(reserve_tor, i){
  diff(tapply(reserve_tor$admit_perf_16_and_above[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_admit_perf_16_and_above <- boot(reserve_tor, statistic = boot_diff_median_admit_perf_16_and_above, R = 10000)

median(b_admit_perf_16_and_above$t)
t(quantile(b_admit_perf_16_and_above$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Unable to ambulate
CrossTable(reserve_tor$assist_oob_walk, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(assist_oob_walk ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Lactate 
skim(subset(reserve_tor$lactate_result_imputed, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$lactate_result_imputed, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(lactate_result_imputed, na.rm = TRUE),
    IQR = IQR(lactate_result_imputed), na.rm = TRUE)

wilcox.test(lactate_result_imputed ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_lactate_result_imputed <- function(reserve_tor, i){
  diff(tapply(reserve_tor$lactate_result_imputed[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_lactate_result_imputed <- boot(reserve_tor, statistic = boot_diff_median_lactate_result_imputed, R = 10000)

median(b_lactate_result_imputed$t)
t(quantile(b_lactate_result_imputed$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Lactate >=4
CrossTable(reserve_tor$lactate4plus, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(lactate4plus ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Microbiology 
#HIV
CrossTable(reserve_tor$hiv_rdt_result, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(hiv_rdt_result ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#HIV clinical stage 3/4
CrossTable(reserve_tor$hiv_stage_3_4_all, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(hiv_stage_3_4_all ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#CD4 count 
skim(subset(reserve_tor$cd4_imputed, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$cd4_imputed, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(cd4_imputed, na.rm = TRUE),
    IQR = IQR(cd4_imputed), na.rm = TRUE)

wilcox.test(cd4_imputed ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_cd4_imputed <- function(reserve_tor, i){
  diff(tapply(reserve_tor$cd4_imputed[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_cd4_imputed <- boot(reserve_tor, statistic = boot_diff_median_cd4_imputed, R = 10000)

median(b_cd4_imputed$t)
t(quantile(b_cd4_imputed$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#HIV viral suppression
CrossTable(reserve_tor$suppressedviralload_imputed, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(suppressedviralload_imputed ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#TB
CrossTable(reserve_tor$microtbdx, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(microtbdx ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#TB LAM in PLWH
#Subset to HIV only for TB-LAM proportion
reserve_tor_hiv <- subset(reserve_tor, hiv_rdt_result==1)

CrossTable(reserve_tor_hiv$urine_lam_result, reserve_tor_hiv$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(urine_lam_result ~ signature_predicted,
                            data = reserve_tor_hiv,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#CRAG in PLWH
#Subset to HIV only for TB-LAM proportion
CrossTable(reserve_tor_hiv$cryto_ag_result, reserve_tor_hiv$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(cryto_ag_result ~ signature_predicted,
                            data = reserve_tor_hiv,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Malaria
CrossTable(reserve_tor$malaria_rdt_result, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(malaria_rdt_result ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Influenza 
CrossTable(reserve_tor$flu_pcr_uvri_result, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(flu_pcr_uvri_result ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#SARS-CoV-2
CrossTable(reserve_tor$covid_pcr_uvri_result, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(covid_pcr_uvri_result ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#MEWS
library(boot)
skim(subset(reserve_tor$mews_score, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$mews_score, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(mews_score, na.rm = TRUE),
    IQR = IQR(mews_score), na.rm = TRUE)

wilcox.test(mews_score ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_mews <- function(reserve_tor, i){
  diff(tapply(reserve_tor$mews_score[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_mews <- boot(reserve_tor, statistic = boot_diff_median_mews, R = 10000)

median(b_mews$t)
t(quantile(b_mews$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#UVA
skim(subset(reserve_tor$uva_score_avpu, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$uva_score_avpu, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(uva_score_avpu, na.rm = TRUE),
    IQR = IQR(uva_score_avpu), na.rm = TRUE)

wilcox.test(uva_score_avpu ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_uva <- function(reserve_tor, i){
  diff(tapply(reserve_tor$uva_score_avpu[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_uva <- boot(reserve_tor, statistic = boot_diff_median_uva, R = 10000)

median(b_uva$t)
t(quantile(b_uva$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#qSOFA >=2
CrossTable(reserve_tor$qsofa_2p_avpu, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(qsofa_2p_avpu ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#Death in hospital
#Create variable
reserve_tor[reserve_tor$hosp_outcome == 2, "deathhosp"] <-1
reserve_tor <- mutate_at(reserve_tor, c("deathhosp"), ~replace(., is.na(.),0))

CrossTable(reserve_tor$deathhosp, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(deathhosp ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#KPS at discharge 
skim(subset(reserve_tor$dc_perf_16_and_above, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$dc_perf_16_and_above, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(dc_perf_16_and_above, na.rm = TRUE),
    IQR = IQR(dc_perf_16_and_above), na.rm = TRUE)

wilcox.test(dc_perf_16_and_above ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_dc_perf_16_and_above <- function(reserve_tor, i){
  diff(tapply(reserve_tor$dc_perf_16_and_above[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_dc_perf_16_and_above <- boot(reserve_tor, statistic = boot_diff_median_dc_perf_16_and_above, R = 10000)

median(b_dc_perf_16_and_above$t)
t(quantile(b_dc_perf_16_and_above$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Death at 30d
CrossTable(reserve_tor$death30d, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(death30d ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#KPS at 30d
skim(subset(reserve_tor$kps_30d_16plus, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$kps_30d_16plus, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(kps_30d_16plus, na.rm = TRUE),
    IQR = IQR(kps_30d_16plus), na.rm = TRUE)

wilcox.test(kps_30d_16plus ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_kps_30d_16plus <- function(reserve_tor, i){
  diff(tapply(reserve_tor$kps_30d_16plus[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_kps_30d_16plus <- boot(reserve_tor, statistic = boot_diff_median_kps_30d_16plus, R = 10000)

median(b_kps_30d_16plus$t)
t(quantile(b_kps_30d_16plus$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Death at 60d
CrossTable(reserve_tor$death60d, reserve_tor$signature_predicted, digits=3, chisq = TRUE)

two_proportion_bootstrap_CI(death60d ~ signature_predicted,
                            data = reserve_tor,
                            first_in_subtraction = "2",
                            response_value_numerator = "1",
                            confidence_level = 0.95,
                            number_repetitions = 10000)

#KPS at 60d
skim(subset(reserve_tor$kps_60d_16plus, reserve_tor$signature_predicted==1))
skim(subset(reserve_tor$kps_60d_16plus, reserve_tor$signature_predicted==2))

group_by(reserve_tor, signature_predicted) %>%
  summarise(
    count = n(),
    median = median(kps_60d_16plus, na.rm = TRUE),
    IQR = IQR(kps_60d_16plus), na.rm = TRUE)

wilcox.test(kps_60d_16plus ~ signature_predicted, data = reserve_tor, exact = FALSE)

boot_diff_median_kps_60d_16plus <- function(reserve_tor, i){
  diff(tapply(reserve_tor$kps_60d_16plus[i], reserve_tor$signature_predicted[i], FUN = median, na.rm = TRUE))}

set.seed(12345)
b_kps_60d_16plus <- boot(reserve_tor, statistic = boot_diff_median_kps_60d_16plus, R = 10000)

median(b_kps_60d_16plus$t)
t(quantile(b_kps_60d_16plus$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

#Logistic regression for association between signature and in-hospital/30d mortality risk 
#In-hospital 
#Univariable
signaturehosp_uv <- glm(deathhosp ~ signature_predicted, data = reserve_tor, family = "binomial")
summary(signaturehosp_uv)
exp(cbind(OR = coef(signaturehosp_uv), confint(signaturehosp_uv)))

#Multivariable 
signaturehosp_mv <- glm(deathhosp ~ signature_predicted + age + gender + illnessduration_to_enroll + hiv_rdt_result, data = reserve_tor, family = "binomial")
summary(signaturehosp_mv)
exp(cbind(OR = coef(signaturehosp_mv), confint(signaturehosp_mv)))

#30d 
#Univariable
signature30d_uv <- glm(death30d ~ signature_predicted, data = reserve_tor, family = "binomial")
summary(signature30d_uv)
exp(cbind(OR = coef(signature30d_uv), confint(signature30d_uv)))

#Interaction between signature and HIV for 30d mortality risk 
signature30d_interact_hiv <- glm(death30d ~ signature_predicted*hiv_rdt_result, data = reserve_tor, family = "binomial")
summary(signature30d_interact_hiv)

#Multivariable 
signature30d_mv <- glm(death30d ~ signature_predicted + age + gender + illnessduration_to_enroll + hiv_rdt_result, data = reserve_tor, family = "binomial")
summary(signature30d_mv)
exp(cbind(OR = coef(signature30d_mv), confint(signature30d_mv)))

#60d 
#Univariable
signature60d_uv <- glm(death60d ~ signature_predicted, data = reserve_tor, family = "binomial")
summary(signature60d_uv)
exp(cbind(OR = coef(signature60d_uv), confint(signature60d_uv)))

#Interaction between signature and HIV for 60d mortality risk 
signature60d_interact_hiv <- glm(death60d ~ signature_predicted*hiv_rdt_result, data = reserve_tor, family = "binomial")
summary(signature60d_interact_hiv)

#Multivariable 
signature60d_mv <- glm(death60d ~ signature_predicted + age + gender + illnessduration_to_enroll + hiv_rdt_result, data = reserve_tor, family = "binomial")
summary(signature60d_mv)
exp(cbind(OR = coef(signature60d_mv), confint(signature60d_mv)))

#Cluster assignment mapped to clinical PCA
clinicalpcamatrix <- as.data.frame(reserve_tor[1:253, c(21,105,312,313,134,137,301,303,304,305,306,309,334,345,354,356)], rownames.force = NA)
head(clinicalpcamatrix)

clinicalpcamatrix <- clinicalpcamatrix %>% mutate_all(~replace(., is.na(.), 0))
clinicalpcamatrix <- clinicalpcamatrix %>% mutate_if(is.character, as.numeric)
clinicalpcamatrix <- clinicalpcamatrix %>% mutate_if(is.factor, as.numeric)

head(clinicalpcamatrix)

#Reformat column names
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "admit_perf_16_and_above"] = "KPS"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "muac_imputed"] = "MUAC"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "lactate_result_imputed"] = "Lactate"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "assist_oob_walk"] = "Unable to stand"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "malaria_rdt_result"] = "Malaria"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "hiv_rdt_result"] = "HIV"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "tempmax"] = "Temp."
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "heart_rate_3"] = "HR"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "resp_rate_3"] = "RR"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "sbp_3"] = "SBP"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "o2_sat_3"] = "SpO2"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "avpu_3"] = "AVPU"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "qsofa_score_avpu"] = "qSOFA"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "mews_score"] = "MEWS"
colnames(clinicalpcamatrix)[colnames(clinicalpcamatrix) == "uva_score_avpu"] = "UVA"
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
                               xlab = "PC1: 24% expl.var",
                               ylab = "PC2: 11% expl.var",
                               col.ind = as.factor(reserve_tor$signature_predicted), 
                               pointshape = 16,
                               alpha.ind = "cos2",
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
                                   xlab = "PC1: 24% expl.var",
                                   ylab = "PC2: 11% expl.var")
vars_clinical_plot <- vars_clinical_plot + theme(text = element_text(size = 16),
                                                 axis.title = element_text(size = 18),
                                                 axis.text = element_text(size = 18))
vars_clinical_plot <- vars_clinical_plot + ggtitle("")
vars_clinical_plot
