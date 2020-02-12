#data source: https://www.cancerrxgene.org/downloads


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
CCLE_expres_f1 = assay(VST_CCLE_expres_f)
CCLE_expres_f = CCLE_expres_f1

#=================================================================
#Regression of drug sensitivity by category variance to remove tissue type effects
#=================================================================
reg = lm(sen~ccle_primary_site+ccle_primary_hist+ccle_hist_subtype_1,data=CCLE_info)
anova(reg)
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
  CCLE_info$sen = CCLE_expres_f[i,]
  reg = lm(sen~ccle_primary_site+ccle_primary_hist+ccle_hist_subtype_1,data=CCLE_info)
  CCLE_expres_f1[i,] = reg$residuals
}
rownames(CCLE_expres_f1) = rownames(CCLE_expres_f)
CCLE_expres_f = CCLE_expres_f1

##################################################################
#WGCNA analysis
##################################################################

#=================================================================
#CCLE_sens_LN (drug sensitivity) and CCLE_expres_f (gene expression) as input for WGCNA analysis
#=================================================================
genenames = rownames(CCLE_expres_f)
CCLE_expres_f = as.data.frame(t(CCLE_expres_f))
names(CCLE_expres_f) = genenames

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Remove bad genes
gsg = goodSamplesGenes(CCLE_expres_f, verbose = 3);
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(CCLE_expres_f)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(CCLE_expres_f)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  CCLE_expres_f = CCLE_expres_f[gsg$goodSamples, gsg$goodGenes]
}

datExpr = CCLE_expres_f

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
save.image("WGCNA_soft_result.RData")

#=================================================================
#Select usable soft-thresholding power
#=================================================================
# Plot the results:
sizeGrWindow(9, 10)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#=================================================================
#auto module construction
#=================================================================
net = blockwiseModules(datExpr, power = 4,
TOMType = "signed", minModuleSize = 30,maxBlockSize = dim(datExpr)[2],
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "femaleMouseTOM",
verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];


# open a graphics window
sizeGrWindow(12, 9)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

#=================================================================
# Relating modules to external traits (datTraits: drug sensitivity data)
#=================================================================
datTraits = CCLE_sens_LN
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels and find interesting modules
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,12)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

#=================================================================
#Select the co-expressed gene module that most negatively correlated with drug sensitivity data
#=================================================================
sort(moduleTraitCor)

ADJ1=abs(cor(CCLE_expres_f,use="p"))^4
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)

#=================================================================
#Output the genes in specific gene modules by the descending order of their intra-module degrees
#=================================================================
darkturquoise_dgree=Alldegrees1[which(moduleColors=="darkturquoise"),"kWithin"]
names(darkturquoise_dgree)= rownames(Alldegrees1)[which(moduleColors=="darkturquoise")]
darkturquoise_dgree=sort(darkturquoise_dgree,decreasing = TRUE)
#head(darkturquoise_dgree)
write(names(darkturquoise_dgree), file='darkturquoise_degree.txt',ncolumns=1)


###################################################
#export the figure to Cytoscape
###################################################
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 4);
# Select modules
modules = c("yellowgreen");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = modProbes;
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
altNodeNames = modGenes,
nodeAttr = moduleColors[inModule]);
