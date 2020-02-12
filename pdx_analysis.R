#============================================================
#Input the Model expression data
#============================================================
cell_expres = read.csv("ExpressionModels.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(cell_expres) = cell_expres[,1]
cell_expres = cell_expres[,-1]

#============================================================
#Calculate PCA score for the gene list
#============================================================
trans_GBM1 = t(cell_expres[c("MYBBP1A" ,"NUP88", "GEMIN4", "PELP1","C1QBP", "DHX33", "WRAP53","TSR1"),])

SNPs = trans_GBM1
pcaSNPs=prcomp(SNPs, scale=TRUE)
str(pcaSNPs$rotation)
summary(pcaSNPs)

a = scale(trans_GBM1) %*% pcaSNPs$rotation[,1]

trans_gene = a[,1]
                        
names(trans_gene) = rownames(a)
write.csv(trans_gene, file="trans_gene_for_models.csv")
