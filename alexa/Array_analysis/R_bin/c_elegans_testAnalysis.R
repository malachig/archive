#Test of normalization strategies on NimbleGen data
#Use data from 7 stages of C. elegans development, provided by Kim Wong

#1.) NORMALIZATION

#A.) RAW DATA
datadir = "/home/malachig/AlternativeSplicing/data_sources/C_elegans_NimbleGenData/RawData";
setwd(datadir)

#Read in the raw data files, which have the following format: 
#X       Y       SEQ_ID  PROBE_ID        X_PIXEL Y_PIXEL HEIGHT  WIDTH   FGD_PIX SIGNAL_MEAN     SIGNAL_STDEV

ma_1_38385 = read.table("38385_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_2_38486 = read.table("38486_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_3_38491 = read.table("38491_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_4_38733 = read.table("38733_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_5_38739 = read.table("38739_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_6_38774 = read.table("38774_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)
ma_7_39127 = read.table("39127_532_clean.ftr", header=T, quote="", sep="\t", as.is=3:4)

#Now combine the raw intensities for all seven experiments into a single dataframe 
#With labels
ma_combined = data.frame(ma_1_38385[,c("SEQ_ID","PROBE_ID","SIGNAL_MEAN")], ma_2_38486[,"SIGNAL_MEAN"], ma_3_38491[,"SIGNAL_MEAN"], ma_4_38733[,"SIGNAL_MEAN"], ma_5_38739[,"SIGNAL_MEAN"], ma_6_38774[,"SIGNAL_MEAN"], ma_7_39127[,"SIGNAL_MEAN"])

#Without labels
raw_data_combined = data.frame(ma_1_38385[,"SIGNAL_MEAN"], ma_2_38486[,"SIGNAL_MEAN"], ma_3_38491[,"SIGNAL_MEAN"], ma_4_38733[,"SIGNAL_MEAN"], ma_5_38739[,"SIGNAL_MEAN"], ma_6_38774[,"SIGNAL_MEAN"], ma_7_39127[,"SIGNAL_MEAN"])

#Rename the columns
names = c("ma_1","ma_2","ma_3","ma_4","ma_5","ma_6","ma_7")
dimnames(raw_data_combined)[[2]] = names

#B.) NIMBLEGEN NORMALIZED

#Compare these data to the data normalized by NimbleGen
datadir = "/home/malachig/AlternativeSplicing/data_sources/C_elegans_NimbleGenData/NormalizedData";
setwd(datadir)
ma_nimble_norm = read.table("All_norm_pair.txt", header=T, quote="", sep="\t", as.is=1:3)
ma_nimble_norm_clean = data.frame(ma_nimble_norm[,"X38385_532"], ma_nimble_norm[,"X38486_532"], ma_nimble_norm[,"X38491_532"], ma_nimble_norm[,"X38733_532"], ma_nimble_norm[,"X38739_532"], ma_nimble_norm[,"X38774_532"], ma_nimble_norm[,"X39127_532"])

names = c("ma_1","ma_2","ma_3","ma_4","ma_5","ma_6","ma_7")
dimnames(ma_nimble_norm_clean)[[2]] = names

#Plot the raw intensities before and after NimbleGen Normalization
#Change to a suitable directory
datadir = "/home/malachig/AlternativeSplicing/Array_analysis/c_elegans_test";
setwd(datadir)

jpeg(file = "f01_raw_intensities_bplot.jpg", quality=100, res=600)
boxplot(raw_data_combined, main = "Figure 1: Raw intensities for 7 hybridizations", ylab="Raw Intensity", xlab="Experiment")
dev.off()

jpeg(file = "f02_nimblegen_norm_intensities_bplot.jpg", quality=100, res=600)
boxplot(ma_nimble_norm_clean, main = "Figure 2: NimbleGen normalized intensities for 7 hybridizations", ylab="Normalized Intensity", xlab="Experiment")
dev.off()


#C.) MANUAL NORMALIZATION USING QUANTILES

#BACKGROUND CORRRECTION - Note: I think that NimbleGen uses the same background correction as RMA (bg.correct.rma()).  However, since this uses an Affy object, which 
#I can not create for NimbleGen data (at least so far I haven't figured out a way) I need to use a generic version of the background correction.  I believe bg.correct.rma()
#uses the generic function bg.adjust(), which accepts a simple matrix as input
x = as.matrix(raw_data_combined)
bg_corrected_data = bg.adjust(x)

#Now attempt to a normalization of the bg.corrected raw data using Quantiles Normalization 
library(affy)
raw_data_quant_norm = normalize.quantiles(bg_corrected_data)
raw_data_quant_norm = as.data.frame(raw_data_quant_norm)

names = c("ma_1","ma_2","ma_3","ma_4","ma_5","ma_6","ma_7")
dimnames(raw_data_quant_norm)[[2]] = names

jpeg(file = "f03_quantiles_norm_intensities_bplot.jpg", quality=100, res=600)
boxplot(raw_data_quant_norm, main = "Figure 3: Quantiles normalized intensities for 7 hybridizations", ylab="Normalized Intensity", xlab="Experiment")
dev.off()

#D.) NORMALIZATION USING VSN
library(vsn)
library(Biobase)
x = as.matrix(raw_data_combined)
eSet = new("exprSet", exprs=x)
raw_data_vsn_norm = vsn(eSet)  #This will take a while!
raw_data_vsn_norm_clean = exprs(raw_data_vsn_norm)
raw_data_vsn_norm_clean = as.data.frame(raw_data_vsn_norm_clean)
jpeg(file = "f04_vsn_norm_intensities_bplot.jpg", quality=100, res=600)
boxplot(raw_data_vsn_norm_clean, main = "Figure 4: VSN Normalized intensities for 7 Hybridizations", ylab="Normalized Intensity", xlab="Experiment")
dev.off()


#2.) EXAMPLE COMPARISONS OF TWO ARRAYS USING ALTERNATE NORMALIZATION METHODS
#Try some M versus A plots to show a data comparison before and after normalization

#i.) Before normalization
M_2vs7 = log2(raw_data_combined[,2] / raw_data_combined[,7])
A_2vs7 = log2(sqrt(raw_data_combined[,2] * raw_data_combined[,7]))
jpeg(file = "f05_rawdata_MvsA_plot.jpg", quality=100, res=600)
plot (A_2vs7, M_2vs7, main="Figure 5: Before Normalization", xlab="A", ylab="M")
dev.off()
#lines(loess.smooth(A_2vs7, M_2vs7), lwd=2, col="red")

#ii.) NimbleGen normalization
M_2vs7 = log2(ma_nimble_norm_clean[,2] / ma_nimble_norm_clean[,7])
A_2vs7 = log2(sqrt(ma_nimble_norm_clean[,2] * ma_nimble_norm_clean[,7]))
jpeg(file = "f06_NimbleGeneNorm_MvsA_plot.jpg", quality=100, res=600)
plot (A_2vs7, M_2vs7, main="Figure6: NimbleGen Normalization", xlab="A", ylab="M")
dev.off()
#lines(loess.smooth(A_2vs7, M_2vs7), lwd=2, col="red")

#iii.) Quantiles normalization
M_2vs7 = log2(raw_data_quant_norm[,2] / raw_data_quant_norm[,7])
A_2vs7 = log2(sqrt(raw_data_quant_norm[,2] * raw_data_quant_norm[,7]))
jpeg(file = "f07_QuantilesNorm_MvsA_plot.jpg", quality=100, res=600)
plot (A_2vs7, M_2vs7, main="Figure 7: Quantiles Normalization", xlab="A", ylab="M")
dev.off()
#lines(loess.smooth(A_2vs7, M_2vs7), lwd=2, col="red")

#iv.) 
M_2vs7 = exprs(raw_data_vsn_norm)[,2] - exprs(raw_data_vsn_norm)[,7]
A_2vs7 = rowMeans(exprs(raw_data_vsn_norm))
jpeg(file = "f08_VsnNorm_MvsA_plot.jpg", quality=100, res=600)
plot (A_2vs7, M_2vs7, main="Figure 8: VSN Normalization", xlab="A", ylab="M")
dev.off()

#Create scatterplots to compare the log ratios of each of the three methods
raw_2vs7 = log2(raw_data_combined[,2] / raw_data_combined[,7])
nimblegen_2vs7 = log2(ma_nimble_norm_clean[,2] / ma_nimble_norm_clean[,7])
quantiles_2vs7 = log2(raw_data_quant_norm[,2] / raw_data_quant_norm[,7])
vsn_2vs7 = exprs(raw_data_vsn_norm)[,2] - exprs(raw_data_vsn_norm)[,7]

jpeg(file = "f09-12_scatter_ratios_plot.jpg", quality=100, res=600)
par(mfrow=c(2,2))
plot (raw_2vs7, main="Figure9: Raw - Log Ratio 2 vs 7", ylab="Log2 Ratio")
plot (nimblegen_2vs7, main="Figure 10: NimbleGen - Log Ratio 2 vs 7", ylab="Log2 Ratio")
plot (quantiles_2vs7, main="Figure 11: Quantiles - Log Ratio 2 vs 7", ylab="Log2 Ratio")
plot (vsn_2vs7, main="Figure 12: VSN - VSN Differences 2 vs 7", ylab="VSN Ratio")
dev.off()

