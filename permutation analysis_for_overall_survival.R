library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(survival)
library(survminer)

#============================================================
#get the gene expression data of patiencs of ovarian cancer from TCGA
#============================================================
load("OV/transcription.Rdata")
rownames(trans_GBM1) = substr(rownames(trans_GBM1), 1, 15)
library("genefilter")
trans_GBM1 = trans_GBM1[-(which(rowSds(data.matrix(trans_GBM1))== 0)),]

#============================================================
#Get protein coding genes
#============================================================
#############################################################
#Note: 1 gene_protein.csv was generated from the "hsapiens_gene_ensembl" dataset of "ENSEMBL_MART_ENSEMBL" Mart, filtering by biotype "protein_coding" 
#      2 Only protein_coding genes were involved in gene_protein.csv
#      3 Gene stable ID	and the corresponding Gene name were included in gene_protein.csv
#############################################################
cod_gene_sym = read.csv("gene_protein.csv", header=T, sep="\t", stringsAsFactors=F, encoding="UTF-8")                                
cod_gene_ccle = intersect(cod_gene_sym[,1], rownames(trans_GBM1))
trans_GBM1 = trans_GBM1[cod_gene_ccle, ] 

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

OV_clinic_data_survival = read.csv("ationwidechildrens.org_clinical_follow_up_v1.0_ov.txt", header=T, sep="\t", stringsAsFactors=F)
patient_id_survival = OV_clinic_data_survival[which((OV_clinic_data_survival[, "vital_status"] == "Dead") |
                                                    (OV_clinic_data_survival[, "vital_status"] == "Alive")) 
                                                    , "bcr_patient_barcode"]

#Set the times for permutation analysis
n = 1000

results = as.data.frame(array(0, dim=c(n,12)))
names(results) = c("logrank.p.value", "hazard.ratio", "lower95", "upper95", "gene 1", "gene 2", "gene 3" ,"gene 4", "gene 5", "gene 6", "gene 7", "gene 8")

#counter for permutation analysis
zc = 1
while (zc<=n)
{  
 
  #select random lists of 8 genes
  m1_gene = rownames(trans_GBM1)
  x=1:length(m1_gene)
  bio_gene = m1_gene[sample(x=x,size=8)][1:8]
  
  gene_ids <- bitr(bio_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  if(length(which(duplicated(gene_ids[,1])==TRUE)) != 0)
  {
    gene_ids = gene_ids[-which(duplicated(gene_ids[,1])==TRUE),]
  }
  #if duplicated genes exist, then re-select genes
  if(dim(gene_ids)[1]<8)
  {
    next
  }  
  rownames(gene_ids) = t(gene_ids["SYMBOL"])
  index = which(gene_ids$ENSEMBL %in% rownames(trans_GBM1))
 
  #calculate the PCA scores for patients from ovarian cancer by these random selected genes
  SNPs = t(trans_GBM1[bio_gene,])
  pcaSNPs=prcomp(SNPs, scale=TRUE)
  str(pcaSNPs$rotation)
  summary(pcaSNPs)

  a = scale(t(trans_GBM1[bio_gene,])) %*% pcaSNPs$rotation[,1]
  trans_gene = t(a[,1])
                        
  names(trans_gene) = rownames(a)
  names(trans_gene) = strtrim(names(trans_gene),12) 


  patient_id_cisplatin_common = intersect(patient_id_cisplatin, names(trans_gene))
  gene_cisplatin = trans_gene[patient_id_cisplatin_common]

  patient_id_survival_single = names(table(patient_id_survival)[which(table(patient_id_survival) == 1)])
  patient_id_survival_multiple = names(table(patient_id_survival)[which(table(patient_id_survival) != 1)])


  #============================================================
  #survival data processing
  #============================================================
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

  survival_array = survival_array[-which(survival_array[, "Days"] == "[Not Available]"), ]
  #################################

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
  survival_array = survival_array[which(survival_array[,"HR_Status"] == "non-HRD"), ]

  for (i in rownames(survival_array))
  {
    if (survival_array[i,"Expression"] >= fivenum(survival_array[,"Expression"])[3])
       survival_array[i, "Group"] = "Low"
    if (survival_array[i,"Expression"] <= fivenum(survival_array[,"Expression"])[3])
       survival_array[i, "Group"] = "High"
  }
  survival_array = survival_array[which((survival_array[,"Group"] == "High") | (survival_array[,"Group"] == "Low")), ]
  survival_array[, "Days"] = as.numeric(survival_array[, "Days"])

  #overall survival analysis
  da = survival_array
  table(da$Status)
  #da$Status = ifelse(da$Status == "Dead", 1, 0)
  survf = survfit(Surv(Days, Status)~Group, data=da)
  ################################
  #cox
  library("survival")
  library("survminer")
  res.cox <- coxph(Surv(Days, Status) ~ Group, data = survival_array)
  #ftest <- cox.zph(res.cox)
  summary(res.cox)
  fit = summary(res.cox) 
  
  hazard.ratio = fit$conf.int[1] 
  lower95 = fit$conf.int[3] 
  upper95 = fit$conf.int[4]
  logrank.p.value = fit$sctest[3]
  
  results[zc,"logrank.p.value"] = logrank.p.value
  results[zc,"hazard.ratio"] = hazard.ratio
  results[zc,"lower95"] = lower95
  results[zc,"upper95"] = upper95
  results[zc,"gene 1"] = gene_ids$SYMBOL[index[1]]
  results[zc,"gene 2"] = gene_ids$SYMBOL[index[2]]
  results[zc,"gene 3"] = gene_ids$SYMBOL[index[3]]
  results[zc,"gene 4"] = gene_ids$SYMBOL[index[4]]
  results[zc,"gene 5"] = gene_ids$SYMBOL[index[5]]
  results[zc,"gene 6"] = gene_ids$SYMBOL[index[6]]
  results[zc,"gene 7"] = gene_ids$SYMBOL[index[7]]
  results[zc,"gene 8"] = gene_ids$SYMBOL[index[8]]

  zc = zc+1

}

write.csv(results, file="1000_times_results_fullgene_nolowHRD.csv")
results = read.csv("1000_times_results_fullgene_nolowHRD.csv", header=T, sep=",", stringsAsFactors=F)
save(results, file = "results.RData")

hazard.ratio = results$hazard.ratio
hist(hazard.ratio, xlim = c(0, 14), col = "black", breaks =100 )
abline(v=4.3336, col = "blue", lwd =2)
sum(results$hazard.ratio > 4.3336)/1000 # one-tailed test
