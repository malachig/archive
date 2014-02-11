
#Load neccessary libraries
library(gcrma)
library(preprocessCore)

#First load a junction expression matrix for known junctions only.
data_file1 = "/projects/alexa2/hmmSplicer/Summary/Compendium/JunctionExpressionMatrix.DA.tsv";
summary_file1 = "/projects/alexa2/hmmSplicer/Summary/Compendium/LibrarySummary.tsv";

#Set output dir
setwd("/projects/alexa2/hmmSplicer/Summary/Compendium/MLL2/")

#Set size of PDF output files
pdfh=12
pdfw=20

#Libraries containing an MLL2 mutation
mll2_mut = c("HS0637","HS0640","HS0641","HS0644","HS0646","HS0647","HS0649","HS0650","HS0656","HS0747","HS0798","HS0804","HS0806","HS0841","HS0842","HS0900","HS0927","HS0930","HS0932","HS0935","HS0939","HS0942","HS0943","HS1135","HS1137","HS1163","HS1181","HS1199","HS1200","HS1203","HS1204","HS1361","HS1462","HS1508","HS1598","HS1978","HS1983","HS2047","HS2048","HS2049","HS2051","HS2055","HS2056","HS2060","HS2250","HS2252","HS2604","HS2606")

######################################################################################################################################################################
#Compendium data first
data = read.table(file=data_file1, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings="na")
summary = read.table(file=summary_file1, header=TRUE, sep="\t", as.is=c(1), na.strings="na")

#Create a normalized expression matrix to allow display of the levels from each library on the same scale - try quantiles normalized values
x = as.matrix(data[,13:(ncol(data)-2)])
x_norm = normalize.quantiles(x)

#MLL2 canonical junctions
mll2_j = which(data[,"Gene_Name"]=="MLL2")

#Calculate expression level of MLL2 by combining all junctions together
mll2_u = apply(x_norm[mll2_j,], 2, median)

lib_ids = colnames(data[,13:(ncol(data)-2)])
mll2 = data.frame(lib_ids, summary[,"Project"], mll2_u)
names(mll2) = c("Library_ID","Project","MLL2_U")

o=order(mll2[,"MLL2_U"])
write.table(file="MLL2_Compendium.tsv", mll2[o,], sep="\t", quote=FALSE);


#Plot expression of MLL2 across all libraries as barplot
pdf(file="MLL2_Compendium.pdf", height=pdfh, width=pdfw)
par(font=2)
bp = barplot(height=mll2[o,"MLL2_U"], col="blue", border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2,
                      xlab="Library", ylab="MLL2 expression (Median normalized read count of 53 known junctions)", main="Expression level of MLL2 across 576 RNA-seq libraries",
                      names.arg=mll2[o,"Project"], cex.names=0.6)
cols=factor(mll2[,"Project"], label=1:12)
rcols=rainbow(12)
mll2[,"color_code"]=cols
mll2[,"color"]=rcols[mll2[,"color_code"]]
points(bp,rep(-1, length(bp)), pch=15, col=mll2[o,"color"], cex=0.50)

#Label DLBCL libs
dlbcl = which(mll2[o,"Project"] == "DLBCL" | mll2[o,"Project"] == "DLBCL_LINE")
points(bp[dlbcl], (mll2[o,"MLL2_U"][dlbcl])+5, pch=8, col="black")

#Label DLBCL libs with an MLL2 mutation 
dlbcl_mut_i = which(mll2[o,"Library_ID"] %in% mll2_mut)
dlbcl_mut_ni = which((!mll2[o,"Library_ID"] %in% mll2_mut) & (mll2[o,"Project"] == "DLBCL" | mll2[o,"Project"] == "DLBCL_LINE"))
#dlbcl_mut_ni = which((!mll2[o,"Library_ID"] %in% mll2_mut) & (mll2[o,"Project"] == "DLBCL"))

points(bp[dlbcl_mut_i], (mll2[o,"MLL2_U"][dlbcl_mut_i])+5, pch=8, col="red")

ltext=levels(mll2[,"Project"])
legend("topleft", ltext, col=rainbow(12), pch=15)
dev.off()

#Compare expression level by t.test
z=t.test(x = mll2[o,"MLL2_U"][dlbcl_mut_i], y = mll2[o,"MLL2_U"][dlbcl_mut_ni], alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
#p-value = 0.03859


#compare MLL2 mutated lib expression versus non-mutated
pdf(file="MLL2_Mutated_vs.NotMutated.pdf", height=7.5, width=7.5)
x = list(mll2[o,"MLL2_U"][dlbcl_mut_i], mll2[o,"MLL2_U"][dlbcl_mut_ni])
names(x) = c("MLL2 Mutated (n=46)","Not Mutated (n=64)")
boxplot(x, col=c("red","gray90"), ylab="MLL2 expression (Median normalized read count of 53 known junctions)", main="Distribution of MLL2 expression in mutated vs. non-mutated libraries")
pv_text = paste("T-test p-value = ", round(z$p.value, digits=4), sep="")
text(x=1.5, y=80, pv_text)
dev.off()


#Plot boxplot of MLL2 expression versus all expression
jpeg(filename="MLL2_Dist.jpg", height=500, width=500)
par(font=2)
tmp=list(log2(x_norm[mll2_j,]+1), log2(x_norm+1))
names(tmp)=c("MLL2 junctions","All junctions")
boxplot(tmp, col=c("red","blue"), ylab="Junction expression (log2 of normalized read counts + 1)")
dev.off();

