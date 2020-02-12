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
#read the drug sensitivity data
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

#============================================================
#generate cell lines catalog for specific cancer type, i.e. BREAST cancer
#============================================================
cell_info = cbind(ori_names,cell_names)
rownames(cell_info) = cell_names
write.csv(grep(pattern = "BREAST", x = cell_info[com_cell,1],value = TRUE), file="BRCA_o.txt")
cell_info_cancer_type = read.csv("BRCA_o.txt", header=T, sep=",", stringsAsFactors=F)
rownames(cell_info_cancer_type) = cell_info_cancer_type[,1]

com_cell = intersect(rownames(cell_info_cancer_type), names(CCLE_expres_f))
CCLE_sens = CCLE_sens[com_cell,] 
CCLE_expres_f = CCLE_expres_f[,com_cell]

#=================================================================
#VST transformation
#=================================================================
library('DESeq2')
colData <- data.frame(condition=factor(names(CCLE_expres_f)))
dds <- DESeqDataSetFromMatrix(data.matrix(CCLE_expres_f), colData, formula(~ condition))
VST_CCLE_expres_f = varianceStabilizingTransformation(dds)
CCLE_expres_f1 = assay(VST_CCLE_expres_f)
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
names(trans_gene) = rownames(a)

cbind(trans_gene, cell_info_cancer_type[,2], CCLE_sens[,"LN_IC50"])

write.csv(cbind(trans_gene, cell_info_cancer_type[,2], CCLE_sens[,"LN_IC50"]), file="BRCA_cell_o.csv")
