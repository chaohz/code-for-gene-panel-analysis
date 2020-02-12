#============================================================
#Input the CCLE expression data
#============================================================
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

#Calculate the sum expression for genes with duplicated records
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
#Olaparib: 
CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1017"), ]
#-----------------------------------------------------------------
#cisplatin: 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1005"), ] 
#-----------------------------------------------------------------
#Rucaparib 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1175"), ]
#-----------------------------------------------------------------
#Talazoparib 
#CCLE_sens = CCLE_sens_o[which(CCLE_sens_o[,"DRUG_ID"] == "1259"), ]
#-----------------------------------------------------------------

#Delete cells with duplicated names
dup_cell = which(duplicated(gsub("-","", CCLE_sens[ ,"CELL_LINE_NAME"])) == TRUE) 
if (length(dup_cell) != 0)
{
  CCLE_sens = CCLE_sens[-dup_cell,]
}
rownames(CCLE_sens) = gsub("-","", CCLE_sens[ ,"CELL_LINE_NAME"]) 
 
CCLE_sens = CCLE_sens[,c("LN_IC50", "AUC")] 
com_cell = intersect(rownames(CCLE_sens), names(CCLE_expres_c))
CCLE_sens = CCLE_sens[com_cell,] 
CCLE_expres_c = CCLE_expres_c[,com_cell]

#Delete all genes with null expression
CCLE_expres_c = CCLE_expres_c[!apply(CCLE_expres_c,1,function(x) all(x==0)),] 
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
CCLE_expres_f = CCLE_expres_f[-(which(rowSds(CCLE_expres_f)== 0)),]

#=================================================================
#VST transformation
#=================================================================
library('DESeq2')
colData <- data.frame(condition=factor(names(CCLE_expres_f)))
dds <- DESeqDataSetFromMatrix(data.matrix(CCLE_expres_f), colData, formula(~ condition))
VST_CCLE_expres_f = varianceStabilizingTransformation(dds)
CCLE_expres_f1 = assay(VST_CCLE_expres_f)
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

CCLE_info$sen = z_sregressed_sen
CCLE_sens_LN = z_sregressed_sen

#=================================================================
#Regression of gene expressions by category variance to remove tissue type effects
#=================================================================
CCLE_expres_f1 = as.data.frame(array(0,dim=c(length(rownames(CCLE_expres_f)),dim(CCLE_expres_f)[2]))) 
for (i in c(1:length(rownames(CCLE_expres_f))))
{
  CCLE_info$AUC = CCLE_expres_f[i,]
  reg = lm(AUC~ccle_primary_site+ccle_primary_hist+ccle_hist_subtype_1,data=CCLE_info)
  CCLE_expres_f1[i,] = reg$residuals
}
rownames(CCLE_expres_f1) = rownames(CCLE_expres_f)
CCLE_expres_f = CCLE_expres_f1

#============================================================
#Calculate PCA score for the gene list
#============================================================
trans_GBM1 = t(CCLE_expres_f[c("MYBBP1A" ,"NUP88", "GEMIN4", "PELP1","C1QBP", "DHX33", "WRAP53","TSR1"),])
SNPs = trans_GBM1
pcaSNPs=prcomp(SNPs, scale=TRUE)
str(pcaSNPs$rotation)
summary(pcaSNPs)

a = scale(trans_GBM1) %*% pcaSNPs$rotation[,1]

trans_gene = a[,1]
                        
names(trans_gene) = rownames(CCLE_sens_LN)


cbind(trans_gene, CCLE_sens_LN)

write.csv(cbind(trans_gene, CCLE_sens_LN), file="result.csv")
