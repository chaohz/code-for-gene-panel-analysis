#============================================================
#get the clinical data of cisplatin treated patiencs of ovarian cancer from TCGA
#============================================================
OV_clinic_data = read.csv("ationwidechildrens.org_clinical_drug_ov.txt", header=T, sep="\t", stringsAsFactors=F)

patient_id_cisplatin = OV_clinic_data[which((OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "CCDP") |
                                      (OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "Ciplastin") |
                                      (OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "Cisplatin") |
                                      (OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "Cisplatin #2-7") |
                                      (OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "Cisplatin/Gemzar") |
                                      (OV_clinic_data[, "pharmaceutical_therapy_drug_name"] == "Cisplatinum"))                                  
                                      , "bcr_patient_barcode"]

#============================================================
#get the gene expression data of patiencs of ovarian cancer from TCGA
#============================================================
load("E:/oncology/A3B/DATA/OV/transcription.Rdata")

#============================================================
#calculate the PCA score of the gene panel
#============================================================
SNPs = t(trans_GBM1[c("ENSG00000132382.13", "ENSG00000108559.10", "ENSG00000179409.9", "ENSG00000141456.13", "ENSG00000108561.7", "ENSG00000141499.15", "ENSG00000005100.11", "ENSG00000167721.9"),])
pcaSNPs=prcomp(SNPs, scale=TRUE)
str(pcaSNPs$rotation)
summary(pcaSNPs)
a = scale(t(trans_GBM1[c("ENSG00000132382.13", "ENSG00000108559.10", "ENSG00000179409.9", "ENSG00000141456.13", "ENSG00000108561.7", "ENSG00000141499.15", "ENSG00000005100.11", "ENSG00000167721.9"),])) %*% pcaSNPs$rotation[,1]
#If the thresholds in pcaSNPs$rotation are negative, add minus sign in the following assignment statement.
trans_gene = -t(a[,1])                        
names(trans_gene) = rownames(a)
names(trans_gene) = strtrim(names(trans_gene),12)

#============================================================
#select patients with both expression data and survival data
#============================================================
patient_id_cisplatin_common = intersect(patient_id_cisplatin, names(trans_gene))
#gene_cisplatin represents the gene panel expression for all patients treated with cisplatin in OV
gene_cisplatin = trans_gene[patient_id_cisplatin_common]

#============================================================
#survival data processing
#============================================================
OV_clinic_data_survival = read.csv("ationwidechildrens.org_clinical_follow_up_v1.0_ov.txt", header=T, sep="\t", stringsAsFactors=F)
patient_id_survival = OV_clinic_data_survival[which((OV_clinic_data_survival[, "vital_status"] == "Dead") |
                                                    (OV_clinic_data_survival[, "vital_status"] == "Alive")) 
                                                    , "bcr_patient_barcode"]

patient_id_survival_single = names(table(patient_id_survival)[which(table(patient_id_survival) == 1)])
patient_id_survival_multiple = names(table(patient_id_survival)[which(table(patient_id_survival) != 1)])

survival_array = as.data.frame(array(0, dim=c((length(patient_id_survival_single) + length(patient_id_survival_multiple)),5)))
names(survival_array) = c("Status", "Days", "Expression", "Group", "HR_Status")
rownames(survival_array) = c(patient_id_survival_single, patient_id_survival_multiple)

#Dead 1; Alive 0
for (i in patient_id_survival_single)
{
  survival_array[i, "Status"] = ifelse(OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i), "vital_status"] == "Dead", 1, 0) 
  survival_array[i, "Days"] = ifelse(survival_array[i, "Status"] == 1, OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i), "death_days_to"], OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i), "last_contact_days_to"])
}

for (i in patient_id_survival_multiple)
{
  num = which(OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i), "vital_status"] == "Dead")[1]
  survival_array[i, "Status"] =  ifelse(is.na(num), 0, 1)
  survival_array[i, "Days"] = ifelse(is.na(num), max(OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i), "last_contact_days_to"]), OV_clinic_data_survival[which(OV_clinic_data_survival[, "bcr_patient_barcode"] == i)[num], "death_days_to"])
}

#delete patients without survival data
survival_array = survival_array[-which(survival_array[, "Days"] == "[Not Available]"), ] 
common_patients = intersect(names(gene_cisplatin), rownames(survival_array))
gene_cisplatin = gene_cisplatin[common_patients]
survival_array = survival_array[common_patients,]
survival_array[, "Expression"] = t(gene_cisplatin)[1,]

#The signature 3 contribution were calculated by MutationalPatterns R package
load("HR_signature.Rdata")
common_patients_HR = intersect(rownames(survival_array), names(HR_sig_num))
survival_array = survival_array[common_patients_HR,]
for (i in rownames(survival_array))
{
  if (HR_sig_num[i] >= fivenum(HR_sig_num[common_patients_HR])[4])
    survival_array[i, "HR_Status"] = "HRD"
  if (HR_sig_num[i] <= fivenum(HR_sig_num[common_patients_HR])[2])
    survival_array[i, "HR_Status"] = "non-HRD"
}
#survival_array = survival_array[which(survival_array[,"HR_Status"] == "non-HRD"), ]
survival_array = survival_array[which((survival_array[,"HR_Status"] == "non-HRD") | (survival_array[,"HR_Status"] == "HRD")), ]

#classify patients according to their expression levels of the gene panel
for (i in rownames(survival_array))
{
  if (survival_array[i,"Expression"] <= fivenum(survival_array[,"Expression"])[3])
     survival_array[i, "Group"] = "Low"
  if (survival_array[i,"Expression"] >= fivenum(survival_array[,"Expression"])[3])
     survival_array[i, "Group"] = "High"
}
survival_array = survival_array[which((survival_array[,"Group"] == "High") | (survival_array[,"Group"] == "Low")), ]
survival_array[, "Days"] = as.numeric(survival_array[, "Days"])

#Overall survival analysis
library(ggplot2)
library(survival)
library(survminer)
####################################
da = survival_array
table(da$Status)
#da$Status = ifelse(da$Status == "Dead", 1, 0)
survf = survfit(Surv(Days, Status)~Group+HR_Status, data=da)
ggsurvplot(survf, conf.int = F, pval = T)
################################
#cox analysis
library("survival")
library("survminer")
res.cox <- coxph(Surv(Days, Status) ~ Group+HR_Status, data = survival_array)
#ftest <- cox.zph(res.cox)
summary(res.cox)
