





#Import normalized gene expression data and fix names
data=read.table("/home/obig/Projects/Colon_Cancer_Resistance/exon_array_data/EnsEMBL_Genes_v53_iterplier_backgroundCorrected_sketchQuantiles.summary.clean.txt", header=TRUE)
names(data) = c("probeset_id","ensembl_g_id","MIP_R1","MIPFU_R1","MIP_R2","MIPFU_R2","MIP_R3","MIPFU_R3","MIP_R6","MIPFU_R6","RKO_R1","RKOFU_R1","RKO_R2","RKO_R3","RKO_R4","RKOFU_R2","RKOFU_R3","RKOFU_R4","HCT_R1","HCT_R2","HCT_R3","HCTFU_R1", "HCTFU_R2", "HCTFU_R3", "MIP101+5FU_R1", "MIP101+5FU_R2","MIP101+5FU_R3","MIP5FUR+5FU_R1","MIP5FUR+5FU_R2","MIP5FUR+5FU_R3")

#Define the 6 classes
mip=c(3,5,7)
mipfu=c(4,6,8)
rko=c(13,14,15)
rkofu=c(16,17,18)
hct=c(19,20,21)
hctfu=c(22,23,24)
all=c(mip, mipfu, rko, rkofu, hct, hctfu)

#Calculate correlations for the data matrix
r = cor(data[,all], method="pearson")
d = 1-r

#Use MDS to create spatial similarity plots based on the input expression values
mds = cmdscale(d, k=2, eig=TRUE)

cols=rainbow(12)
col1=c(rep("black",length(data[1,]))) 
col1[c(mip)]=cols[1]
col1[c(mipfu)]=cols[2]
col1[c(rko)]=cols[5]
col1[c(rkofu)]=cols[7]
col1[c(hct)]=cols[10]
col1[c(hctfu)]=cols[12]

setwd("/home/malachig/")
tiff(filename="AffyExonArray_MDS_SimilarityPlot.tiff", width = 480, height = 480, compression = c("none"), bg="white", type="cairo")
par(font.main=2, font.axis=2, font.lab=2)
fixed_names = c("MIP101","","","","","MIP/5FU","RKO","","","RKO/5FU","","","HCT116","","","HCT/5FU","","")
plot(mds$points, type="n", xlab="Similarity (x)", ylab="Similarity (y)", main="Gene expression profile similarity (6 cell lines, 3 replicates each)")
points(mds$points[,1], mds$points[,2], col=col1[all], pch=20, cex=2)
#text(mds$points[,1], mds$points[,2], fixed_names, col=col1[all], cex=1, font=2)

text (-0.55, -0.02, "MIP101", col=cols[1], cex=1.2, font=2);  
text (-0.55, +0.02, "MIP/5FU", col=cols[2], cex=1.2, font=2);  

text (0.2, -0.1, "RKO", col=cols[5], cex=1.2, font=2);  
text (0.2, -0.075, "RKO/5FU", col=cols[7], cex=1.2, font=2);  

text (0.2, 0.11, "HCT116", col=cols[10], cex=1.2, font=2);  
text (0.2, 0.07, "HCT/5FU", col=cols[12], cex=1.2, font=2);  

dev.off()


