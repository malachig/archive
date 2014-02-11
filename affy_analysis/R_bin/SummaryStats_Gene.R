
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Affymetrix_Exon_Arrays/Sensitive_vs_Resistant/ExpressionConsoleResults"
dir(datadir)
setwd(datadir)

gene_plier = read.table("iterplier-gene-CORE_EnsEMBL_45_36g.summary.txt", header=T, quote="", sep="\t", comment.char="#")
names(gene_plier) = c("alexa_gene_id","MIP101_Rep1","MIP5FUR_Rep1","MIP101_Rep2","MIP5FUR_Rep2","MIP101_Rep3","MIP5FUR_Rep3","MIP101_Rep6","MIP5FUR_Rep6","RKO_Rep1","RKO5FUR_Rep1")

#determine DE of (RKO - RKO-5FUR)
#*** NOTE these are not background corrected, this is just a priliminary test!!***

gene_plier[,"DE_RKO_vs_RKO5FUR"] = log2(gene_plier[,"RKO_Rep1"] + 16) - log2(gene_plier[,"RKO5FUR_Rep1"] + 16)



#Summarize gene values - already log2 scale

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_5-29-2007/genelevel_ensembl_CORE-quant-antibgp-Rep6_only/Rep6_mRNA/import/normalize/summarize/"
dir(datadir)
setwd(datadir)
data_norm_plier_transcript = read.table("plier.summary.log2.Rep1-3_6.txt", header=T, quote="", sep="\t", comment.char="#")

summary(data_norm_plier_transcript[,"mip101.1_log2"])
summary(data_norm_plier_transcript[,"mip101.2_log2"])
summary(data_norm_plier_transcript[,"mip101.3_log2"])
summary(data_norm_plier_transcript[,"mip101.6_log2"])

summary(data_norm_plier_transcript[,"mip5fur.1_log2"])
summary(data_norm_plier_transcript[,"mip5fur.2_log2"])
summary(data_norm_plier_transcript[,"mip5fur.3_log2"])
summary(data_norm_plier_transcript[,"mip5fur.6_log2"])

summary(data_norm_plier_transcript[,"mip_vs_5fur.1_DE"])
summary(data_norm_plier_transcript[,"mip_vs_5fur.2_DE"])
summary(data_norm_plier_transcript[,"mip_vs_5fur.3_DE"])
summary(data_norm_plier_transcript[,"mip_vs_5fur.6_DE"])

par(mfrow=c(4,1))
hist(data_norm_plier_transcript[,"mip_vs_5fur.1_DE"], breaks=110, density=100, col="dark green", main="Gene Level DE MIP101-1 versus 5FUR-1", xlab="Log2 DE Intensity", ylim=c(0,3100), xlim=c(-5,5))
hist(data_norm_plier_transcript[,"mip_vs_5fur.2_DE"], breaks=110, density=100, col="dark green", main="Gene Level DE MIP101-2 versus 5FUR-2", xlab="Log2 DE Intensity", ylim=c(0,3100), xlim=c(-5,5))
hist(data_norm_plier_transcript[,"mip_vs_5fur.3_DE"], breaks=110, density=100, col="dark green", main="Gene Level DE MIP101-3 versus 5FUR-3", xlab="Log2 DE Intensity", ylim=c(0,3100), xlim=c(-5,5))
hist(data_norm_plier_transcript[,"mip_vs_5fur.6_DE"], breaks=110, density=100, col="dark green", main="Gene Level DE MIP101-3 versus 5FUR-3", xlab="Log2 DE Intensity", ylim=c(0,3100), xlim=c(-5,5))

length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.1_DE"]) > 1))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.2_DE"]) > 1))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.3_DE"]) > 1))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.6_DE"]) > 1))

length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.1_DE"]) > 2))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.2_DE"]) > 2))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.3_DE"]) > 2))
length(which(abs(data_norm_plier_transcript[,"mip_vs_5fur.6_DE"]) > 2))





