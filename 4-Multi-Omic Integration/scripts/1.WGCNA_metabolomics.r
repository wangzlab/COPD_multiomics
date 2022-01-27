# the required input file: 
# Metabolomic feature abundance (metabolome.txt): An abundance table with compound ID in rows and samples in columns.

# output files generated:
# Co-occurring module abundance tables with prefix being "1_metaB" stored in the "Output" directory



library(WGCNA)
library(flashClust)
allowWGCNAThreads()
enableWGCNAThreads()


min.ModuleSize = 10
cut.Height = 0.1
output.prefix = "1_metaB"
outputDir = "Output"

if(!dir.exists(outputDir)) dir.create(outputDir)

dataExpr <- read.table("metabolome.txt", 
                       sep='\t', row.names=1, header=T, quote="", comment="", check.names=F)
dim(dataExpr)


#### filter those mad < 0.01 ###############
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
#dataExpr <- as.data.frame(t(dataExpr))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)


powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType="signed", verbose=5)

######### plot power curves ####################
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf(paste(outputDir,"/",output.prefix,".power_curve.pdf",sep = ""))

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
############ plot power curves ################

softPower <- sft$powerEstimate

############ scale free examination #############
#k <- softConnectivity(datE=dataExpr,power=softPower) 
#sizeGrWindow(10, 5)
#par(mfrow=c(1,2))
#hist(k)
#pdf("2_scale_free_plot.pdf")
#scaleFreePlot(k,main="Check Scale free topology\n")
#dev.off()
########### scale free examination ##############

adjacency = adjacency (dataExpr, power = softPower, type = "signed", corFnc = "bicor", corOptions = "use = 'pairwise.complete.obs'")
#adjacency = adjacency (dataExpr, power = softPower, type = "signed", corFnc = "cor", corOptions = "use = 'p', method = 'spearman'")

TOM = TOMsimilarity(adjacency);

dissTOM = 1-TOM
geneTree = flashClust (as.dist (dissTOM), method = "complete")

#geneTree = hclust(as.dist(dissTOM), method = "average");


# minModuleSize = 10;

# try both pamRespectsDendro=F and T, see which one is better #
moduleLabels1 = cutreeDynamic (dendro = geneTree, distM = dissTOM, method = "hybrid", deepSplit = 2, pamRespectsDendro = F, minClusterSize = 5)
moduleLabels1 = labels2colors (moduleLabels1)
table(moduleLabels1)

#mergeCutHeight: 
#################### merge modules based on cutHeight #################
merge = mergeCloseModules (dataExpr, moduleLabels1, corFnc = "cor", corOptions = list (use = 'p', method = 'spearman'), cutHeight = cut.Height)

moduleLabels2 = merge$colors


MEList = moduleEigengenes(dataExpr, colors = moduleLabels2)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "complete");

pdf(paste(outputDir,"/",output.prefix,".eigengene_heatmap.pdf", sep = "") )
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), plotDendrograms = TRUE, xLabelsAngle = 90) 
dev.off()
pdf(paste(outputDir,"/",output.prefix,".gene_cluster.pdf",sep = "") )
plotDendroAndColors(geneTree, cbind(moduleLabels1, moduleLabels2), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColorsMeta = moduleLabels2
names (moduleColorsMeta) = colnames (dataExpr)
MEsMeta = orderMEs (MEs)
rownames (MEsMeta) = rownames (dataExpr)

###################save results #####################
write.table(moduleColorsMeta, paste(outputDir,"/",output.prefix, ".module_assign.txt", sep = "") , sep="\t", append=FALSE, quote=FALSE)
write.table(MEsMeta, paste(outputDir,"/",output.prefix, ".module_eigengene.txt", sep = ""), sep="\t", append=FALSE, quote=FALSE)

save(TOM,file=paste(outputDir,"/",output.prefix,".tom.rda",sep = "") )
