#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Compare Affymetrix Exon Array,NimbleGen Array and Illumina paired-end data for HOUSEKEEPING GENES AND NEGATIVE CONTROL PROBES
#START BY USING RAW PROBE-LEVEL DATA or AVERAGE COVERAGE DATA EXPRESSION ESTIMATES
#THE SAME ~100 GENES USED AS CONTROLS BY AFFY WERE INCLUDED IN THE NIMBLEGEN DESIGN (93 WHICH MAP UNAMBIGUOUSLY) AND CAN BE ASSESSED IN ILLUMINA

#1.) RAW DATA ANALYSIS - NIMBLEGEN VERSUS AFFY VERSUS ILLUMINA
#Head-to-head comparison.  Make sure the same list of genes is being compared in both cases
#Get this from the list of Affy control data 
# -this is the overlap between Alexa and Affy Hk genes that have at least some probes in both platforms


#Import housekeeping expression data
datadir = "/projects/malachig/solexa/cross_platform_comparisons/housekeeping_genes/"
dir(datadir)
setwd(datadir)

#1-A.) Import Affymetrix Exon array data - Hybidization intensities for Intron, Exon and Random sequence probes ('antigenomic' control probes)
data_affy_hk = read.table("Affymetrix_data_HKandNC_ControlsOnly.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(2,6), na.strings='na')
master_hk_gene_list = unique(data_affy_hk[,"AlexaGene_ID"])
master_hk_gene_list = master_hk_gene_list[1:93]

#1-B.) Import NimbleGen data - Hybridization intensities for Intron, Exon and Random sequence probes ('purely random sequences')
data_nimble_complete = read.table("NimbleGen_ALEXA_data_Complete.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,14,17), na.strings='na')

#Create a dataframe which has only the negative control gene probes
data_nimble_nc = data_nimble_complete[which(data_nimble_complete[,"Probe_Type"] == "Control-Negative"),]

#data_nimble_hk = data_nimble_complete[which(data_nimble_complete[,"AlexaGene_ID"] %in% master_hk_gene_list),]
#write.table (data_nimble_hk, sep="\t", file="NimbleGen_ALEXA_data_HKandNC_ControlsOnly.txt", quote=FALSE, row.names=FALSE)
#Create a dataframe which has only the overlaping housekeeping gene probes
data_nimble_hk = read.table("NimbleGen_ALEXA_data_HKandNC_ControlsOnly.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,14,17), na.strings='na')

#1-C.) Import Illumina data - Average sequence coverage values for Intron, Exon and Intergenic (Silent) regions
#Note that intergenic sequences are not a fair comparison to random probes 
#One estimates genomic DNA contamination, the other estimates non-specific hybridization noise
#The design of both array platforms does not allow for estimation of genomic background...

#Use data for MIP101 library only
illumina_hk_exon1 = read.table("HS04391_Illumina_HK_Exons.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3,4,6), na.strings='NA')
illumina_hk_intron1 = read.table("HS04391_Illumina_HK_Introns.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3,4), na.strings='NA')
illumina_hk_intergenic1 = read.table("HS04391_Illumina_HK_Intergenics.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3), na.strings='NA')
illumina_hk_exon2 = read.table("HS04401_Illumina_HK_Exons.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3,4,6), na.strings='NA')
illumina_hk_intron2 = read.table("HS04401_Illumina_HK_Introns.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3,4), na.strings='NA')
illumina_hk_intergenic2 = read.table("HS04401_Illumina_HK_Intergenics.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(1,2,3), na.strings='NA')


#2.) Eliminate expression estimates from Illumina that are based on very small exons (i.e. some exon have 0 unmasked bases!)
#These will not generate good expression estimates and correspond to regions that would be skipped in a microarray design
illumina_hk_exon1 = illumina_hk_exon1[which(illumina_hk_exon1[,"Exon_Region_Length"] > 50),]
illumina_hk_intron1 = illumina_hk_intron1[which(illumina_hk_intron1[,"UnMasked_Base_Count"] > 50),]
illumina_hk_intergenic1 = illumina_hk_intergenic1[which(illumina_hk_intergenic1[,"UnMasked_Base_Count"] > 50),]
illumina_hk_exon2 = illumina_hk_exon2[which(illumina_hk_exon2[,"Exon_Region_Length"] > 50),]
illumina_hk_intron2 = illumina_hk_intron2[which(illumina_hk_intron2[,"UnMasked_Base_Count"] > 50),]
illumina_hk_intergenic2 = illumina_hk_intergenic2[which(illumina_hk_intergenic2[,"UnMasked_Base_Count"] > 50),]

#3.) NORMALIZE REPLICATES WITHIN EACH ARRAY PLATFORM
library(gcrma)
names = c("MIP_C","MIP_GH_A","MIP_EF_B","FUR_C","FUR_GH_A","FUR_EF_B")
data_combined = data.frame(data_affy_hk[,"MIP_C"],data_affy_hk[,"MIP_GH_A"],data_affy_hk[,"MIP_EF_B"],data_affy_hk[,"FUR_C"],data_affy_hk[,"FUR_GH_A"],data_affy_hk[,"FUR_EF_B"])
names(data_combined) = names
x = as.matrix(data_combined)
data_quant_norm = normalize.quantiles(x)
data_quant_norm = as.data.frame(data_quant_norm)
dimnames(data_quant_norm)[[2]] = names
data_affy_hk = data.frame(data_affy_hk[,"Probe_Type"], data_affy_hk[,"Probe_Tm"], data_quant_norm)
names(data_affy_hk) = c("Probe_Type", "Probe_Tm", names)

names = c("MIP_EF_A","MIP_GH_A","MIP_EF_B","FUR_EF_A","FUR_GH_A","FUR_EF_B")
data_combined = data.frame(data_nimble_hk[,"MIP_EF_A"],data_nimble_hk[,"MIP_GH_A"],data_nimble_hk[,"MIP_EF_B"],data_nimble_hk[,"FUR_EF_A"],data_nimble_hk[,"FUR_GH_A"],data_nimble_hk[,"FUR_EF_B"])
names(data_combined) = names
x = as.matrix(data_combined)
data_quant_norm = normalize.quantiles(x)
data_quant_norm = as.data.frame(data_quant_norm)
dimnames(data_quant_norm)[[2]] = names
data_nimble_hk = data.frame(data_nimble_hk[,"Probe_Type"], data_nimble_hk[,"Probe_Tm"], data_quant_norm)
names(data_nimble_hk) = c("Probe_Type", "Probe_Tm", names)

names = c("MIP_EF_A","MIP_GH_A","MIP_EF_B","FUR_EF_A","FUR_GH_A","FUR_EF_B")
data_combined = data.frame(data_nimble_nc[,"MIP_EF_A"],data_nimble_nc[,"MIP_GH_A"],data_nimble_nc[,"MIP_EF_B"],data_nimble_nc[,"FUR_EF_A"],data_nimble_nc[,"FUR_GH_A"],data_nimble_nc[,"FUR_EF_B"])
names(data_combined) = names
x = as.matrix(data_combined)
data_quant_norm = normalize.quantiles(x)
data_quant_norm = as.data.frame(data_quant_norm)
dimnames(data_quant_norm)[[2]] = names
data_nimble_nc = data.frame(data_nimble_nc[,"Probe_Type"], data_nimble_nc[,"Probe_Tm"], data_quant_norm)
names(data_nimble_nc) = c("Probe_Type", "Probe_Tm", names)


#CHANGE TO AN OUTPUT DIRECTORY
datadir = "/projects/malachig/solexa/cross_platform_comparisons/housekeeping_genes/results/"
dir(datadir)
setwd(datadir)

#4.) Summarize the distribution of expression estimates for Exons, Introns, and Controls for all three platforms
#Create dataframes which have only the housekeeping EXON probes
nimble_hk_exon = data_nimble_hk[which(data_nimble_hk[,"Probe_Type"] == "Exon"),]
affy_hk_exon = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Exon"),]

#Create a dataframe which has only the housekeeping INTRON probes
nimble_hk_intron = data_nimble_hk[which(data_nimble_hk[,"Probe_Type"] == "Intron"),]
affy_hk_intron = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Intron"),]

#Create a dataframe which has only the negative control gene probes
nimble_nc = data_nimble_nc[which(data_nimble_nc[,"Probe_Type"] == "Control-Negative"),]
affy_nc = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Control-Negative"),]

#Summarize the number of data points for each array platform
length(affy_hk_exon[,1])
length(nimble_hk_exon[,1])
length(illumina_hk_exon1[,1])

length(affy_hk_intron[,1])
length(nimble_hk_intron[,1])
length(illumina_hk_intron1[,1])

length(affy_nc[,1])
length(nimble_nc[,1])
length(illumina_hk_intergenic1[,1])


#Create box plots of the same comparisons
box_hk_nimble = function(exp){
  title = paste ("NimbleGen - Exon and intron expression estimates for 100 house-keeping genes")
  probe_list = list(log2(nimble_hk_exon[,exp]+1), log2(nimble_hk_intron[,exp]+1), log2(nimble_nc[,exp]+1))
  names(probe_list) = c("Exon","Intron","Random")
  boxplot(probe_list, main=title, ylab="Log2 Intensity", ylim=c(0,16), col=c("blue","green","red"), col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.2)
  return(1)
}
box_hk_affy = function(exp){
  title = paste ("Affymetrix - Exon and intron expression estimates for 100 house-keeping genes")
  probe_list = list(log2(affy_hk_exon[,exp]+1), log2(affy_hk_intron[,exp]+1), log2(affy_nc[,exp]+1))
  names(probe_list) = c("Exon","Intron","Random")
  boxplot(probe_list, main=title, ylab="Log2 Intensity", ylim=c(0,16), col=c("blue","green","red"), col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.2)
  return(1)
}
box_hk_illumina = function(exp){
  title = paste ("Illumina - Exon and intron expression estimates for 100 house-keeping genes")
  probe_list = list(log2(illumina_hk_exon1[,exp]+1), log2(illumina_hk_intron1[,exp]+1), log2(illumina_hk_intergenic1[,exp]+1))
  names(probe_list) = c("Exon","Intron","Intergenic")
  boxplot(probe_list, main=title, ylab="Log2 Intensity", ylim=c(0,16), col=c("blue","green","red"), col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.2)
  return(1)
}

tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_RepA.tiff", width=600, height=600, compression="none")
par(mfrow=c(3,1), font.main = 2, font.lab = 2, font.axis=2)
 box_hk_affy("MIP_EF_B")
 box_hk_nimble("MIP_EF_B")
 box_hk_illumina("Average_Coverage_RAW")
dev.off()

tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_RepB.tiff", width=600, height=600, compression="none")
par(mfrow=c(3,1), font.main = 2, font.lab = 2)
 box_hk_affy("MIP_GH_A")
 box_hk_nimble("MIP_GH_A")
 box_hk_illumina("Average_Coverage_RAW")
dev.off()

tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_RepC.tiff", width=600, height=600, compression="none")
par(mfrow=c(3,1), font.main = 2, font.lab = 2)
 box_hk_affy("MIP_C")
 box_hk_nimble("MIP_EF_A")
 box_hk_illumina("Average_Coverage_RAW")
dev.off()

#Now combine data from all three platforms into one big boxplot
box_hk_all = function(affy_name, nimble_name, illumina_name){
  title = paste ("Exon and intron expression estimates for 100 house-keeping genes")
  probe_list = list(log2(affy_hk_exon[,affy_name]+1), log2(affy_hk_intron[,affy_name]+1), log2(affy_nc[,affy_name]+1),
                    log2(nimble_hk_exon[,nimble_name]+1), log2(nimble_hk_intron[,nimble_name]+1), log2(nimble_nc[,nimble_name]+1),
                    log2(illumina_hk_exon1[,illumina_name]+1), log2(illumina_hk_intron1[,illumina_name]+1), log2(illumina_hk_intergenic1[,illumina_name]+1)
                    )
  names(probe_list) = c("Exon","Intron","Random","Exon","Intron","Random","Exon","Intron","Intergenic")
  boxplot(probe_list, main=title, ylab="Log2(Expression+1)", ylim=c(0,17), col=c("blue","green","red"), col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4, las=2)
  abline(v=3.5, lwd=1, col="black")
  abline(v=6.5, lwd=1, col="black")
  return(1)
}
tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_COMBO_RepA.tiff", width=800, height=600, compression="none")
par(font.main=2, font.lab=2, mar=(c(7, 4, 4, 2)+0.1))
box_hk_all("MIP_EF_B","MIP_EF_B","Average_Coverage_RAW")
text(2,17, "Affymetrix", cex=1.4, font=2)
text(5,17, "NimbleGen", cex=1.4, font=2)
text(8,17, "Illumina", cex=1.4, font=2)
dev.off()

tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_COMBO_RepB.tiff", width=800, height=600, compression="none")
par(font.main=2, font.lab=2, mar=(c(7, 4, 4, 2)+0.1))
box_hk_all("MIP_GH_A","MIP_GH_A","Average_Coverage_RAW")
text(2,17, "Affymetrix", cex=1.4, font=2)
text(5,17, "NimbleGen", cex=1.4, font=2)
text(8,17, "Illumina", cex=1.4, font=2)
dev.off()

tiff(file="HouseKeeping_ExonIntronNC_BoxPlots_3platforms_COMBO_RepC.tiff", width=800, height=600, compression="none")
par(font.main=2, font.lab=2, mar=(c(7, 4, 4, 2)+0.1))
box_hk_all("MIP_C","MIP_EF_A","Average_Coverage_RAW")
text(2,17, "Affymetrix", cex=1.4, font=2)
text(5,17, "NimbleGen", cex=1.4, font=2)
text(8,17, "Illumina", cex=1.4, font=2)
dev.off()


#5.) Calculate the signal-to-noise ratio for each platform

#Calculate the median of all HK Exon probe intensities for all hybes
nimble_exon_medians = as.numeric(lapply (list((nimble_hk_exon[,"MIP_EF_A"]), (nimble_hk_exon[,"MIP_EF_B"]), (nimble_hk_exon[,"MIP_GH_A"]), 
                                   (nimble_hk_exon[,"FUR_EF_A"]), (nimble_hk_exon[,"FUR_EF_B"]), (nimble_hk_exon[,"FUR_GH_A"])), median))

affy_exon_medians = as.numeric(lapply (list((affy_hk_exon[,"MIP_C"]), (affy_hk_exon[,"MIP_EF_B"]), (affy_hk_exon[,"MIP_GH_A"]), 
                                 (affy_hk_exon[,"FUR_C"]), (affy_hk_exon[,"FUR_EF_B"]), (affy_hk_exon[,"FUR_GH_A"])), median))

illumina_exon_medians = as.numeric(lapply (list((illumina_hk_exon1[,"Average_Coverage_RAW"]), illumina_hk_exon2[,"Average_Coverage_RAW"]), median))


#Calculate the median of all HK Intron probe intensities for all hybes
nimble_intron_medians = as.numeric(lapply (list((nimble_hk_intron[,"MIP_EF_A"]), (nimble_hk_intron[,"MIP_EF_B"]), (nimble_hk_intron[,"MIP_GH_A"]), 
                                     (nimble_hk_intron[,"FUR_EF_A"]), (nimble_hk_intron[,"FUR_EF_B"]), (nimble_hk_intron[,"FUR_GH_A"])), median))

affy_intron_medians = as.numeric(lapply (list((affy_hk_intron[,"MIP_C"]), (affy_hk_intron[,"MIP_EF_B"]), (affy_hk_intron[,"MIP_GH_A"]), 
                                   (affy_hk_intron[,"FUR_C"]), (affy_hk_intron[,"FUR_EF_B"]), (affy_hk_intron[,"FUR_GH_A"])), median))

illumina_intron_medians = as.numeric(lapply (list((illumina_hk_intron1[,"Average_Coverage_RAW"]), illumina_hk_intron2[,"Average_Coverage_RAW"]), median))


nimble_ei = (nimble_exon_medians / nimble_intron_medians)
affy_ei = (affy_exon_medians / affy_intron_medians)
illumina_ei = (illumina_exon_medians / illumina_intron_medians)

#Affy
mean(affy_ei)
sqrt(var(affy_ei))

#NimbleGen
mean(nimble_ei)
sqrt(var(nimble_ei))

#Illumina
mean(illumina_ei)
sqrt(var(illumina_ei))


#Sample sizes are too small to prove normality
library(nortest)
ad.test(nimble_ei)
ad.test(affy_ei)

#Are these signal-to-noise ratios significantly different between ALEXA and AFFY??
#t.test(x=nimble_ei, y = affy_ei, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
wilcox.test(x=nimble_ei, y=affy_ei, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, 
	     conf.int = FALSE, conf.level = 0.95)


#6.) Calculate sensitivity/specificity across the entire set of housekeeping genes
#First create a data frame to store results for each experiment/platform
cutoffs = seq(0, 16, 0.0025)

calculate_roc_scores = function(exon_data, intron_data){

  x = array(0, c(length(cutoffs),2), dimnames=list(cutoffs, c("sensitivity","specificity")))
  result = as.data.frame(x) 

  true_pos = length(exon_data)
  true_neg = length(intron_data)

  #Iterate through posible cutoff intensity values
  for(i in 1:length(cutoffs)){

    cutoff = cutoffs[i]

    obs_pos = length(which(log2(exon_data+1) >= cutoff))
    obs_neg = length(which(log2(intron_data+1) <= cutoff))

    #TP = length(which(exon_probes[,exp] >= cutoff))
    #FN = length(which(exon_probes[,exp] < cutoff))
 
    #TN = length(which(intron_probes[,exp] <= cutoff))
    #FP = length(which(intron_probes[,exp] > cutoff))

    result[i,"sensitivity"] = obs_pos/true_pos
    result[i,"specificity"] = obs_neg/true_neg

    #result[i,"sensitivity"] = TP/(TP+FN)
    #result[i,"specificity"] = TP/(TP+FP)

  }
  return(result)
}

#Also Calculate the Area Under the Curve for each line
library(caTools)

#AFFY
affy_roc_mip_ef_b = calculate_roc_scores(affy_hk_exon[,"MIP_EF_B"], affy_hk_intron[,"MIP_EF_B"])
 affy_auc_mip_ef_b = trapz(affy_roc_mip_ef_b[,"specificity"],affy_roc_mip_ef_b[,"sensitivity"])
affy_roc_mip_gh_a = calculate_roc_scores(affy_hk_exon[,"MIP_GH_A"], affy_hk_intron[,"MIP_GH_A"])
 affy_auc_mip_gh_a = trapz(affy_roc_mip_gh_a[,"specificity"],affy_roc_mip_gh_a[,"sensitivity"])
affy_roc_mip_c = calculate_roc_scores(affy_hk_exon[,"MIP_C"], affy_hk_intron[,"MIP_C"])
 affy_auc_mip_c = trapz(affy_roc_mip_c[,"specificity"],affy_roc_mip_c[,"sensitivity"])

affy_roc_fur_ef_b = calculate_roc_scores(affy_hk_exon[,"FUR_EF_B"], affy_hk_intron[,"FUR_EF_B"])
 affy_auc_fur_ef_b = trapz(affy_roc_fur_ef_b[,"specificity"],affy_roc_fur_ef_b[,"sensitivity"])
affy_roc_fur_gh_a = calculate_roc_scores(affy_hk_exon[,"FUR_GH_A"], affy_hk_intron[,"FUR_GH_A"])
 affy_auc_fur_gh_a = trapz(affy_roc_fur_gh_a[,"specificity"],affy_roc_fur_gh_a[,"sensitivity"])
affy_roc_fur_c = calculate_roc_scores(affy_hk_exon[,"FUR_C"], affy_hk_intron[,"FUR_C"])
 affy_auc_fur_c = trapz(affy_roc_fur_c[,"specificity"],affy_roc_fur_c[,"sensitivity"])

#NIMBLEGEN
nimble_roc_mip_ef_b = calculate_roc_scores(nimble_hk_exon[,"MIP_EF_B"], nimble_hk_intron[,"MIP_EF_B"])
 nimble_auc_mip_ef_b = trapz(nimble_roc_mip_ef_b[,"specificity"],nimble_roc_mip_ef_b[,"sensitivity"])
nimble_roc_mip_gh_a = calculate_roc_scores(nimble_hk_exon[,"MIP_GH_A"], nimble_hk_intron[,"MIP_GH_A"])
 nimble_auc_mip_gh_a = trapz(nimble_roc_mip_gh_a[,"specificity"],nimble_roc_mip_gh_a[,"sensitivity"])
nimble_roc_mip_ef_a = calculate_roc_scores(nimble_hk_exon[,"MIP_EF_A"], nimble_hk_intron[,"MIP_EF_A"]) 
 nimble_auc_mip_ef_a = trapz(nimble_roc_mip_ef_a[,"specificity"],nimble_roc_mip_ef_a[,"sensitivity"])

nimble_roc_fur_ef_b = calculate_roc_scores(nimble_hk_exon[,"FUR_EF_B"], nimble_hk_intron[,"FUR_EF_B"])
 nimble_auc_fur_ef_b = trapz(nimble_roc_fur_ef_b[,"specificity"],nimble_roc_fur_ef_b[,"sensitivity"])
nimble_roc_fur_gh_a = calculate_roc_scores(nimble_hk_exon[,"FUR_GH_A"], nimble_hk_intron[,"FUR_GH_A"])
 nimble_auc_fur_gh_a = trapz(nimble_roc_fur_gh_a[,"specificity"],nimble_roc_fur_gh_a[,"sensitivity"])
nimble_roc_fur_ef_a = calculate_roc_scores(nimble_hk_exon[,"FUR_EF_A"], nimble_hk_intron[,"FUR_EF_A"])
 nimble_auc_fur_ef_a = trapz(nimble_roc_fur_ef_a[,"specificity"],nimble_roc_fur_ef_a[,"sensitivity"])

#ILLUMINA
illumina_roc_mip = calculate_roc_scores(illumina_hk_exon1[,"Average_Coverage_RAW"], illumina_hk_intron1[,"Average_Coverage_RAW"])
 illumina_auc_mip = trapz(illumina_roc_mip[,"specificity"],illumina_roc_mip[,"sensitivity"])

illumina_roc_fur = calculate_roc_scores(illumina_hk_exon2[,"Average_Coverage_RAW"], illumina_hk_intron2[,"Average_Coverage_RAW"])
 illumina_auc_fur = trapz(illumina_roc_fur[,"specificity"],illumina_roc_fur[,"sensitivity"])

#Create ROC graphs with the data generated above
#MIP Data
#Create names for each experiment that include the AUC value
tiff(file="HouseKeeping_ExonIntronNC_ROC_3platforms_MIP101_Reps.tiff", width=700, height=700, compression="none")
par(font.main = 2, font.lab = 2)
name_n_EF_A = paste("Nimble_EF_A", "(AUC =", round(nimble_auc_mip_ef_a, digit=4), ")")
name_n_EF_B = paste("Nimble_EF_B", "(AUC =", round(nimble_auc_mip_ef_b, digit=4), ")")
name_n_GH_A = paste("Nimble_GH_A", "(AUC =", round(nimble_auc_mip_gh_a, digit=4), ")")
name_a_C = paste("Affy_C", "(AUC =", round(affy_auc_mip_c, digit=4), ")")
name_a_EF_B = paste("Affy_EF_B", "(AUC =", round(affy_auc_mip_ef_b, digit=4), ")")
name_a_GH_A = paste("Affy_GH_A", "(AUC =", round(affy_auc_mip_gh_a, digit=4), ")")
name_i = paste("Illumina", "(AUC =", round(illumina_auc_mip, digit=4), ")")

plot(1-(nimble_roc_mip_ef_a[,"specificity"]), nimble_roc_mip_ef_a[,"sensitivity"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affy vs. NimbleGen - 3 Hybes Each (Sensitive data)", lty=1, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(1-(nimble_roc_mip_ef_b[,"specificity"]), nimble_roc_mip_ef_b[,"sensitivity"], col="blue", type="l", lty=2, lwd=2)
lines(1-(nimble_roc_mip_gh_a[,"specificity"]), nimble_roc_mip_gh_a[,"sensitivity"], col="blue", type="l", lty=3, lwd=2)
lines(1-(affy_roc_mip_c[,"specificity"]), affy_roc_mip_c[,"sensitivity"], col="red", type="l", lty=1, lwd=2)
lines(1-(affy_roc_mip_ef_b[,"specificity"]), affy_roc_mip_ef_b[,"sensitivity"], col="red", type="l", lty=2, lwd=2)
lines(1-(affy_roc_mip_gh_a[,"specificity"]), affy_roc_mip_gh_a[,"sensitivity"], col="red", type="l", lty=3, lwd=2)
lines(1-(illumina_roc_mip[,"specificity"]), illumina_roc_mip[,"sensitivity"], col="dark green", type="l", lty=1, lwd=2)
legend("bottomright", c(name_n_EF_A,name_n_EF_B,name_n_GH_A,name_a_C,name_a_EF_B,name_a_GH_A,name_i), lty=c(1,2,3,1,2,3,1), col=c("blue","blue","blue","red","red","red","dark green"), lwd=2)
dev.off()

#5FUR Data
#Create names for each experiment that include the AUC value
tiff(file="HouseKeeping_ExonIntronNC_ROC_3platforms_MIP5FUR_Reps.tiff", width=700, height=700, compression="none")
par(font.main = 2, font.lab = 2)
name_n_EF_A = paste("Nimble_EF_A", "(AUC =", round(nimble_auc_fur_ef_a, digit=4), ")")
name_n_EF_B = paste("Nimble_EF_B", "(AUC =", round(nimble_auc_fur_ef_b, digit=4), ")")
name_n_GH_A = paste("Nimble_GH_A", "(AUC =", round(nimble_auc_fur_gh_a, digit=4), ")")
name_a_C = paste("Affy_C", "(AUC =", round(affy_auc_fur_c, digit=4), ")")
name_a_EF_B = paste("Affy_EF_B", "(AUC =", round(affy_auc_fur_ef_b, digit=4), ")")
name_a_GH_A = paste("Affy_GH_A", "(AUC =", round(affy_auc_fur_gh_a, digit=4), ")")
name_i = paste("Illumina", "(AUC =", round(illumina_auc_fur, digit=4), ")")

plot(1-(nimble_roc_fur_ef_a[,"specificity"]), nimble_roc_fur_ef_a[,"sensitivity"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affy vs. NimbleGen - 3 Hybes Each (Resistant data)", lty=1, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(1-(nimble_roc_fur_ef_b[,"specificity"]), nimble_roc_fur_ef_b[,"sensitivity"], col="blue", type="l", lty=2, lwd=2)
lines(1-(nimble_roc_fur_gh_a[,"specificity"]), nimble_roc_fur_gh_a[,"sensitivity"], col="blue", type="l", lty=3, lwd=2)
lines(1-(affy_roc_fur_c[,"specificity"]), affy_roc_fur_c[,"sensitivity"], col="red", type="l", lty=1, lwd=2)
lines(1-(affy_roc_fur_ef_b[,"specificity"]), affy_roc_fur_ef_b[,"sensitivity"], col="red", type="l", lty=2, lwd=2)
lines(1-(affy_roc_fur_gh_a[,"specificity"]), affy_roc_fur_gh_a[,"sensitivity"], col="red", type="l", lty=3, lwd=2)
lines(1-(illumina_roc_fur[,"specificity"]), illumina_roc_fur[,"sensitivity"], col="dark green", type="l", lty=1, lwd=2)
legend("bottomright", c(name_n_EF_A,name_n_EF_B,name_n_GH_A,name_a_C,name_a_EF_B,name_a_GH_A,name_i), lty=c(1,2,3,1,2,3,1), col=c("blue","blue","blue","red","red","red","dark green"), lwd=2)
dev.off()

#Combine all sensitivity/specificity scores (median at each cutoff) for the two platforms
#For each cutoff, also calculate the SD values
#Plot the median ROC curves and corresponding SD curves for both platforms

nimble_roc_all = cbind(nimble_roc_mip_ef_b, nimble_roc_mip_gh_a, nimble_roc_mip_ef_a, nimble_roc_fur_ef_b, nimble_roc_fur_gh_a, nimble_roc_fur_ef_a)
names(nimble_roc_all) = c("sen_mip_ef_b", "spec_mip_ef_b", "sen_mip_gh_a", "spec_mip_gh_a", "sen_mip_ef_a", "spec_mip_ef_a", "sen_fur_ef_b", "spec_fur_ef_b", "sen_fur_gh_a", "spec_fur_gh_a", "sen_fur_ef_a", "spec_fur_ef_a")

affy_roc_all = cbind(affy_roc_mip_ef_b, affy_roc_mip_gh_a, affy_roc_mip_c, affy_roc_fur_ef_b, affy_roc_fur_gh_a, affy_roc_fur_c)
names(affy_roc_all) = c("sen_mip_ef_b", "spec_mip_ef_b", "sen_mip_gh_a", "spec_mip_gh_a", "sen_mip_ef_a", "spec_mip_ef_a", "sen_fur_ef_b", "spec_fur_ef_b", "sen_fur_gh_a", "spec_fur_gh_a", "sen_fur_ef_a", "spec_fur_ef_a")

illumina_roc_all = cbind(illumina_roc_mip, illumina_roc_fur)
names(illumina_roc_all) = c("sen_mip","spec_mip","sen_fur","spec_fur")

calculate_mean_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = mean(sen)
  return(mean_sen)
}
calculate_mean_sens2 = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip"], roc_data["sen_fur"]))
  mean_sen = mean(sen)
  return(mean_sen)
}
nimble_roc_all[,"mean_sens"] = apply(nimble_roc_all, 1, calculate_mean_sens)
affy_roc_all[,"mean_sens"] = apply(affy_roc_all, 1, calculate_mean_sens)
illumina_roc_all[,"mean_sens"] = apply(illumina_roc_all, 1, calculate_mean_sens2)

calculate_mean_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = mean(spec)
  return(mean_spec)
}
calculate_mean_spec2 = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip"], roc_data["spec_fur"]))
  mean_spec = mean(spec)
  return(mean_spec)
}
nimble_roc_all[,"mean_spec"] = apply(nimble_roc_all, 1, calculate_mean_spec)
affy_roc_all[,"mean_spec"] = apply(affy_roc_all, 1, calculate_mean_spec)
illumina_roc_all[,"mean_spec"] = apply(illumina_roc_all, 1, calculate_mean_spec2)

nimble_auc_all = trapz(nimble_roc_all[,"mean_spec"],nimble_roc_all[,"mean_sens"])
affy_auc_all = trapz(affy_roc_all[,"mean_spec"],affy_roc_all[,"mean_sens"])
illumina_auc_all = trapz(illumina_roc_all[,"mean_spec"],illumina_roc_all[,"mean_sens"])

calculate_lower_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  lower_sen = mean_sen - sd_sen
  return(lower_sen)
}
calculate_lower_sens2 = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip"], roc_data["sen_fur"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  lower_sen = mean_sen - sd_sen
  return(lower_sen)
}
nimble_roc_all[,"lower_sens"] = apply(nimble_roc_all, 1, calculate_lower_sens)
affy_roc_all[,"lower_sens"] = apply(affy_roc_all, 1, calculate_lower_sens)
illumina_roc_all[,"lower_sens"] = apply(illumina_roc_all, 1, calculate_lower_sens2)

calculate_upper_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  upper_sen = mean_sen + sd_sen
  return(upper_sen)
}
calculate_upper_sens2 = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip"], roc_data["sen_fur"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  upper_sen = mean_sen + sd_sen
  return(upper_sen)
}
nimble_roc_all[,"upper_sens"] = apply(nimble_roc_all, 1, calculate_upper_sens)
affy_roc_all[,"upper_sens"] = apply(affy_roc_all, 1, calculate_upper_sens)
illumina_roc_all[,"upper_sens"] = apply(illumina_roc_all, 1, calculate_upper_sens2)

calculate_lower_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  lower_spec = mean_spec - sd_spec
  return(lower_spec)
}
calculate_lower_spec2 = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip"], roc_data["spec_fur"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  lower_spec = mean_spec - sd_spec
  return(lower_spec)
}
nimble_roc_all[,"lower_spec"] = apply(nimble_roc_all, 1, calculate_lower_spec)
affy_roc_all[,"lower_spec"] = apply(affy_roc_all, 1, calculate_lower_spec)
illumina_roc_all[,"lower_spec"] = apply(illumina_roc_all, 1, calculate_lower_spec2)

calculate_upper_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  upper_spec = mean_spec + sd_spec
  return(upper_spec)
}
calculate_upper_spec2 = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip"], roc_data["spec_fur"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  upper_spec = mean_spec + sd_spec
  return(upper_spec)
}
nimble_roc_all[,"upper_spec"] = apply(nimble_roc_all, 1, calculate_upper_spec)
affy_roc_all[,"upper_spec"] = apply(affy_roc_all, 1, calculate_upper_spec)
illumina_roc_all[,"upper_spec"] = apply(illumina_roc_all, 1, calculate_upper_spec2)


#Determine the max combined specificity and sensitivity achieved by each platform
affy_roc_all[,"SS_combined"] = affy_roc_all[,"mean_sens"] + affy_roc_all[,"mean_spec"]
max_combo = max(affy_roc_all[,"SS_combined"])
affy_roc_all[which(affy_roc_all[,"SS_combined"]==max_combo),c("mean_sens","mean_spec")]
xy_affy = c((1 - affy_roc_all[which(affy_roc_all[,"SS_combined"]==max_combo),"mean_spec"]), (affy_roc_all[which(affy_roc_all[,"SS_combined"]==max_combo),"mean_sens"]))

nimble_roc_all[,"SS_combined"] = nimble_roc_all[,"mean_sens"] + nimble_roc_all[,"mean_spec"]
max_combo = max(nimble_roc_all[,"SS_combined"])
nimble_roc_all[which(nimble_roc_all[,"SS_combined"]==max_combo),c("mean_sens","mean_spec")]
xy_nimble = c((1 - nimble_roc_all[which(nimble_roc_all[,"SS_combined"]==max_combo),"mean_spec"]), (nimble_roc_all[which(nimble_roc_all[,"SS_combined"]==max_combo),"mean_sens"]))

illumina_roc_all[,"SS_combined"] = illumina_roc_all[,"mean_sens"] + illumina_roc_all[,"mean_spec"]
max_combo = max(illumina_roc_all[,"SS_combined"])
illumina_roc_all[which(illumina_roc_all[,"SS_combined"]==max_combo),c("mean_sens","mean_spec")]
xy_illumina = c((1 - illumina_roc_all[which(illumina_roc_all[,"SS_combined"]==max_combo),"mean_spec"]), (illumina_roc_all[which(illumina_roc_all[,"SS_combined"]==max_combo),"mean_sens"]))


#Create a plot of mean sensitivity/specificity as an ROC curve.  Include the S.D. at each point and the AUC values calculated by the trapz function

#Create names for each experiment that include the AUC value
tiff(file="HouseKeeping_ExonIntronNC_ROC_3platforms_Means.tiff", width=700, height=700, compression="none")
par(font.main = 2, font.lab = 2)
name_alexa = paste("ALEXA", "(AUC =", round(nimble_auc_all, digit=3), ")")
name_affy = paste("Affymetrix", "(AUC =", round(affy_auc_all, digit=3), ")")
name_illumina = paste("Illumina", "(AUC =", round(illumina_auc_all, digit=3), ")")

plot(1-(nimble_roc_all[,"mean_spec"]), nimble_roc_all[,"mean_sens"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affymetrix vs. NimbleGen vs. WTSS", lty=1, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(1-(nimble_roc_all[,"lower_spec"]), nimble_roc_all[,"lower_sens"], col="blue", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
lines(1-(nimble_roc_all[,"upper_spec"]), nimble_roc_all[,"upper_sens"], col="blue", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))

lines(1-(affy_roc_all[,"mean_spec"]), affy_roc_all[,"mean_sens"], col="red", type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,1))
lines(1-(affy_roc_all[,"lower_spec"]), affy_roc_all[,"lower_sens"], col="red", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
lines(1-(affy_roc_all[,"upper_spec"]), affy_roc_all[,"upper_sens"], col="red", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))

lines(1-(illumina_roc_all[,"mean_spec"]), illumina_roc_all[,"mean_sens"], col="dark green", type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,1))
lines(1-(illumina_roc_all[,"lower_spec"]), illumina_roc_all[,"lower_sens"], col="dark green", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
lines(1-(illumina_roc_all[,"upper_spec"]), illumina_roc_all[,"upper_sens"], col="dark green", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
legend("bottomright", c(name_alexa,"S.D.",name_affy,"S.D.",name_illumina,"S.D."), lty=c(1,2,1,2,1,2), col=c("blue","blue","red","red","dark green", "dark green"), lwd=c(2,1,2,1,2,1))
dev.off()



tiff(file="HouseKeeping_ExonIntronNC_ROC_3platforms_Means_Simple.tiff", width=700, height=700, compression="none")
par(font.main = 2, font.lab = 2)
name_alexa ="NimbleGen (ALEXA arrays)"
name_affy = "Affymetrix (exon arrays)"
name_illumina = "Illumina (WTSS)"
name_max = "Point of maximum Sensitivity/Specificity"

plot(1-(nimble_roc_all[,"mean_spec"]), nimble_roc_all[,"mean_sens"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affymetrix vs. NimbleGen vs. Illumina", lty=1, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(1-(affy_roc_all[,"mean_spec"]), affy_roc_all[,"mean_sens"], col="red", type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,1))
lines(1-(illumina_roc_all[,"mean_spec"]), illumina_roc_all[,"mean_sens"], col="dark green", type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,1))

points(xy_affy[1], xy_affy[2], pch=16, col="black")
points(xy_nimble[1], xy_nimble[2], pch=16, col="black")
points(xy_illumina[1], xy_illumina[2], pch=16, col="black")
legend("bottomright", c(name_illumina,name_alexa,name_affy,name_max), lty=c(1,1,1,0), pch=16, col=c("dark green", "blue", "red", "black"), lwd=c(2,2,2,0), cex=1.4)
dev.off()











