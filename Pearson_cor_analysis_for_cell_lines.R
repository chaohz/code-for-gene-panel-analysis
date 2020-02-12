#SANGER Drug data
#data source: https://www.cancerrxgene.org/downloads
#CCLE raw expression data
#data source: https://portals.broadinstitute.org/ccle/data

#============================================================
#Input the CCLE expression data
#============================================================
#input ccle expression data
CCLE_expres = read.csv("CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct", header=T, sep="\t", stringsAsFactors=F)
 

#Get cell names without tissue type info
ori_names = names(CCLE_expres)[-c(1:2)]
cell_names = character()
for (i in c(1:length(ori_names)))
{
  cell_names[i] = strsplit(ori_names, "[_.]")[[i]][1]
}
names(CCLE_expres)[-c(1:2)] = cell_names

#Get gene ensembl id
resem_id = substr(CCLE_expres$Name, 1, 15)
rownames(CCLE_expres) = resem_id

#============================================================
#Get protein coding genes
#============================================================
#############################################################
#Note: 1 gene_protein.csv was generated from the "hsapiens_gene_ensembl" dataset of "ENSEMBL_MART_ENSEMBL" Mart, filtering by biotype "protein_coding" 
#      2 Only protein_coding genes were involved in gene_protein.csv
#      3 Gene stable ID	and the corresponding Gene name were included in gene_protein.csv
#############################################################
cod_gene_sym = read.csv("gene_protein.csv", header=T, sep="\t", stringsAsFactors=F, encoding="UTF-8")                               
cod_gene_ccle = intersect(cod_gene_sym[,1], resem_id)
CCLE_expres = CCLE_expres[cod_gene_ccle, ] 

#============================================================
#Calculate the sum expression for genes with duplicated records
#============================================================
CCLE_expres$Description = factor(CCLE_expres$Description) 
CCLE_expres_c = as.data.frame(array(0, dim=c(length(levels(CCLE_expres$Description)),dim(CCLE_expres[-c(1:2)])[2]))) 
names(CCLE_expres_c) = names(CCLE_expres)[-c(1:2)]
rownames(CCLE_expres_c) = levels(CCLE_expres$Description)
 
for (i in c(1:dim(CCLE_expres[-c(1:2)])[2]))
{
  CCLE_expres_c[,i]=tapply(CCLE_expres[-c(1:2)][,i], CCLE_expres$Description, sum)
}

#============================================================
#Input the drug sensitivity data
#============================================================
CCLE_sens_o = read.csv("v17.3_fitted_dose_response.txt", header=T, sep="\t", stringsAsFactors=F)

#=================================================================
#select drug 
#parp inhibitor: olaparib-1
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1495"), ]
#-----------------------------------------------------------------
#DNA alkylating agent: Temozolomide
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1375"), ]
#-----------------------------------------------------------------
#Olaparib: 
CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1017"), ]
#-----------------------------------------------------------------
#cisplatin: 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1005"), ] 
#-----------------------------------------------------------------
#5-Fluorouracil: 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "179"), ] 
#-----------------------------------------------------------------
#Rucaparib 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1175"), ]
#-----------------------------------------------------------------
#Talazoparib 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1259"), ]
#-----------------------------------------------------------------
#Mitomycin-C
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "136"), ]
#-----------------------------------------------------------------
#JQ12
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "164"), ]
#-----------------------------------------------------------------
#Bleomycin
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "190"), ]
#-----------------------------------------------------------------
#Etoposide
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "134"), ]
#-----------------------------------------------------------------
#SN-38
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1494"), ]
#-----------------------------------------------------------------

#Delete cells with duplicated names
dup_cell = which(duplicated(gsub("-","", CCLE_sens[ ,"CELL_LINE_NAME"])) == TRUE) 
if (length(dup_cell) != 0)
{
  CCLE_sens = CCLE_sens[-dup_cell,]
}
rownames(CCLE_sens) = gsub("-","", CCLE_sens[ ,"CELL_LINE_NAME"]) 
 
#=================================================================
#Get cell lines with both expression data and drug sensitivity data
#=================================================================
CCLE_sens = CCLE_sens[,c("LN_IC50", "AUC")] 
com_cell = intersect(rownames(CCLE_sens), names(CCLE_expres_c))
CCLE_sens = CCLE_sens[com_cell,] 
CCLE_expres_c = CCLE_expres_c[,com_cell]

#Delete all genes with null expression
CCLE_expres_c = CCLE_expres_c[!apply(CCLE_expres_c,1,function(x) all(x==0)),]
CCLE_sens = CCLE_sens[names(CCLE_expres_c),] 
CCLE_expres_f = CCLE_expres_c 

rm(CCLE_expres, CCLE_sens_o, CCLE_expres_c)

#=================================================================
#Get the tissue sources of cell lines
#=================================================================
CCLE_cell_name = read.csv("v21.meta.per_cell_line.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(CCLE_cell_name) = CCLE_cell_name[,"ccl_name"]
CCLE_info = as.data.frame(array(0, dim=c(dim(CCLE_sens)[1],4)))
names(CCLE_info) = c("ccle_primary_site", "ccle_primary_hist", "ccle_hist_subtype_1", "AUC")
CCLE_info$ccle_primary_site = CCLE_cell_name[rownames(CCLE_sens), "ccle_primary_site"]
CCLE_info$ccle_primary_hist = CCLE_cell_name[rownames(CCLE_sens), "ccle_primary_hist"]
CCLE_info$ccle_hist_subtype_1 = CCLE_cell_name[rownames(CCLE_sens), "ccle_hist_subtype_1"]
CCLE_info$sen = CCLE_sens[,"LN_IC50"]

#Delete cells without tissue sources info
del_cell = which(is.na(CCLE_info$ccle_primary_site)) 
CCLE_info = CCLE_info[-del_cell,] 
CCLE_expres_f = CCLE_expres_f[,-del_cell] 
CCLE_sens = CCLE_sens[-del_cell,]

#Delete genes whose standard deviation estimate for all cell lines is zero
library("genefilter")
if (length(which(rowSds(CCLE_expres_f)== 0)) != 0)
{
  CCLE_expres_f = CCLE_expres_f[-(which(rowSds(CCLE_expres_f)== 0)),]
}

#=================================================================
#VST transformation
#=================================================================
library('DESeq2')
colData <- data.frame(condition=factor(names(CCLE_expres_f)))
dds <- DESeqDataSetFromMatrix(data.matrix(CCLE_expres_f), colData, formula(~ condition))
VST_CCLE_expres_f = varianceStabilizingTransformation(dds)
#glm.fit: algorithm did not converge
CCLE_expres_f1 = assay(VST_CCLE_expres_f)
#save(list=ls(), file="Veliparib/vst.Rdata") 
CCLE_expres_f = CCLE_expres_f1

#=================================================================
#Regression of drug sensitivity by category variance to remove tissue type effects
#=================================================================
reg = lm(sen~ccle_primary_site+ccle_primary_hist+ccle_hist_subtype_1,data=CCLE_info)
#Use the Residuals: the tissue type effects have been removed by using the residuals
regressed_sen = reg$residuals
names(regressed_sen) = rownames(CCLE_sens)
#Z scale transformation
z_sregressed_sen = scale(regressed_sen)

#=================================================================
#Regression of gene expressions by category variance to remove tissue type effects
#=================================================================
CCLE_info$sen = z_sregressed_sen
CCLE_sens_LN = z_sregressed_sen
CCLE_expres_f1 = as.data.frame(array(0,dim=c(length(rownames(CCLE_expres_f)),dim(CCLE_expres_f)[2]))) 
for (i in c(1:length(rownames(CCLE_expres_f))))
{
  CCLE_info$sen = CCLE_expres_f[i,]
  reg = lm(sen~ccle_primary_site+ccle_primary_hist+ccle_hist_subtype_1,data=CCLE_info)
  CCLE_expres_f1[i,] = reg$residuals
}
rownames(CCLE_expres_f1) = rownames(CCLE_expres_f)
CCLE_expres_f = CCLE_expres_f1

#=================================================================
#calculate pearson correlation
#=================================================================
Cor_result = as.data.frame(array(0,dim=c(length(rownames(CCLE_expres_f)),3))) 
for (i in c(1:length(rownames(CCLE_expres_f))))
{
  result = cor.test(as.numeric(CCLE_expres_f[i,]), CCLE_sens_LN, alternative = "two.sided", method = "pearson",conf.level = 0.95)
  Cor_result[i,1] = result$estimate
  Cor_result[i,2] = result$p.value
}

names(Cor_result) = c("R_value", "P_value", "Q_value")
rownames(Cor_result) = rownames(CCLE_expres_f) 

z_cor_result = scale(Cor_result$R_value)
 
rownames(z_cor_result) = rownames(Cor_result)

#calculate the FDR based on the p values derived from Pearson correlation
library(qvalue)
Cor_result[,3] = qvalue(p=Cor_result$P_value)$qvalue

write.csv(data.frame(z_cor_result, Cor_result$P_value, Cor_result$Q_value), file="Pearson_results.csv")
