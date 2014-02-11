#COMPARE AFFY AND NIMBLEGEN RAW DATA FOR HOUSEKEEPING GENES AND NEGATIVE CONTROL PROBES
#START BY USING RAW PROBE-LEVEL DATA FOR BOTH
#THE SAME ~100 GENES USED AS CONTROLS BY AFFY WERE INCLUDED IN THE NIMBLEGEN DESIGN (93 WHICH MAP UNAMBIGUOUSLY)

#1.) RAW DATA ANALYSIS - NIMBLEGEN

#Specify the working directory 
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw"
dir(datadir)
setwd(datadir)

#Now input the actual intensity measurements measurements 
#File contains data for multiple hybridizations

#IMPORT DATA FILE 
#Note: Convert 'na' values to NA.  Use 'as.is' lines with strings

#Complete file
data_complete = read.table("All_hybes_withProbeInfo.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,12,15), na.strings='na')

#Now import the list of housekeeping genes
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006"
dir(datadir)
setwd(datadir)

hk_gene_list = read.table("affy_alexa_controls_genes_final.txt", header=T, quote="", sep="\t", comment.char="", na.strings='na')

#CHANGE TO AN OUTPUT DIRECTORY
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/Figures"
dir(datadir)
setwd(datadir)

#Create a dataframe which has only the housekeeping gene probes
data_hk_genes = data_complete[which(data_complete[,"AlexaGene_ID"] %in% hk_gene_list[,"alexa_id"]),]

#Create a dataframe which has only the housekeeping EXON probes
data_hk_exon = data_hk_genes[which(data_hk_genes[,"Probe_Type"] == "Exon"),]

#Create a dataframe which has only the housekeeping INTRON probes
data_hk_intron = data_hk_genes[which(data_hk_genes[,"Probe_Type"] == "Intron"),]

#Create a dataframe which has only the negative control gene probes
data_nc_genes = data_complete[which(data_complete[,"Probe_Type"] == "Control-Negative"),]

#Create a scatter plot of housekeeping exon and intron probe intensities and contrast with negative-control probes
#Use MIP_GH_A as an example

draw_hk = function(exp){
  title = paste ("Intensity of Housekeeping Exon and Intron Probes (", exp, ")")
  plot(data_nc_genes[,"Probe_Tm"], log2(data_nc_genes[,exp]), ylim=c(6,18), col="red", ylab="RAW Log2 Probe Intensity",
       xlab="Probe Tm", main=title, pch=1)
  points(data_hk_intron[,"Probe_Tm"], log2(data_hk_intron[,exp]), ylim=c(6,18), col="green", pch=1)
  points(data_hk_exon[,"Probe_Tm"], log2(data_hk_exon[,exp]), ylim=c(6,18), col="blue", pch=1)
  legend(locator(n=1,type="n"), c("Exonic","Intronic","Random Sequences"), pch=c(1,1,1), col=c("blue","green","red"))
  return(1)
}
draw_hk("MIP_GH_A")

#Create boxplots of the same groups
box_hk = function(exp){
  title = paste ("Intensity of Housekeeping Exon and Intron Probes (", exp, ")")
  probe_list = list(log2(data_hk_exon[,exp]), log2(data_hk_intron[,exp]), log2(data_nc_genes[,exp]))
  names(probe_list) = c("Exonic","Intronic","Random")
  boxplot(probe_list, main=title, ylab="Log2 Intensity")
  return(1)
}
par(mfrow=c(2,1))
  box_hk("MIP_EF_A")
  box_hk("FUR_EF_A")
par(mfrow=c(2,1))
  box_hk("MIP_EF_B")
  box_hk("FUR_EF_B")
par(mfrow=c(2,1))
  box_hk("MIP_GH_A")
  box_hk("FUR_GH_A")

#Calculate the mean of all HK Exon log2 probe intensities for all hybes
exon_means = lapply (list(log2(data_hk_exon[,"MIP_EF_A"]), log2(data_hk_exon[,"MIP_EF_B"]), log2(data_hk_exon[,"MIP_GH_A"]), 
                          log2(data_hk_exon[,"FUR_EF_A"]), log2(data_hk_exon[,"FUR_EF_B"]), log2(data_hk_exon[,"FUR_GH_A"])), mean)

#Calculate the mean of all HK Intron log2 probe intensities for all hybes
intron_means = lapply (list(log2(data_hk_intron[,"MIP_EF_A"]), log2(data_hk_intron[,"MIP_EF_B"]), log2(data_hk_intron[,"MIP_GH_A"]), 
                            log2(data_hk_intron[,"FUR_EF_A"]), log2(data_hk_intron[,"FUR_EF_B"]), log2(data_hk_intron[,"FUR_GH_A"])), mean)


#Calculate the median of all HK Exon log2 probe intensities for all hybes
exon_medians = lapply (list(log2(data_hk_exon[,"MIP_EF_A"]), log2(data_hk_exon[,"MIP_EF_B"]), log2(data_hk_exon[,"MIP_GH_A"]), 
                            log2(data_hk_exon[,"FUR_EF_A"]), log2(data_hk_exon[,"FUR_EF_B"]), log2(data_hk_exon[,"FUR_GH_A"])), median)

#Calculate the median of all HK Intron log2 probe intensities for all hybes
intron_medians = lapply (list(log2(data_hk_intron[,"MIP_EF_A"]), log2(data_hk_intron[,"MIP_EF_B"]), log2(data_hk_intron[,"MIP_GH_A"]), 
                              log2(data_hk_intron[,"FUR_EF_A"]), log2(data_hk_intron[,"FUR_EF_B"]), log2(data_hk_intron[,"FUR_GH_A"])), median)


#2.) RAW DATA ANALYSIS - NIMBLEGEN VERSUS AFFY
#Head-to-head comparison.  Make sure the same list of genes is being compared in both cases
#Get this from the list of Affy control data 
# -this is the overlap between Alexa and Affy Hk genes that have at least some probes in both platforms


#Import Affy housekeeping gene data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel"
dir(datadir)
setwd(datadir)

data_affy_hk = read.table("HuEX-1_001-006_RawProbe_HK+NC_ControlsOnly.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(2,6), na.strings='na')
master_hk_gene_list = unique(data_affy_hk[,"AlexaGene_ID"])

#Import NimbleGen data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw/"
dir(datadir)
setwd(datadir)

#Raw data - Before smoothing!!
data_nimble_complete = read.table("All_hybes_withProbeInfo.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,14,17), na.strings='na')

#After NMPP smoothing
#data_nimble_complete = read.table("All_hybes_withProbeInfo_smoothed.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,14,17), na.strings='na')

#Create a dataframe which has only the overlaping housekeeping gene probes
data_nimble_hk = data_nimble_complete[which(data_nimble_complete[,"AlexaGene_ID"] %in% master_hk_gene_list),]

#NORMALIZE DATASETS WITHIN EACH PLATFORM
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



#CHANGE TO AN OUTPUT DIRECTORY
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/Figures/HousekeepingGenes"
dir(datadir)
setwd(datadir)

#Create dataframes which have only the housekeeping EXON probes
nimble_hk_exon = data_nimble_hk[which(data_nimble_hk[,"Probe_Type"] == "Exon"),]
affy_hk_exon = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Exon"),]

#Create a dataframe which has only the housekeeping INTRON probes
nimble_hk_intron = data_nimble_hk[which(data_nimble_hk[,"Probe_Type"] == "Intron"),]
affy_hk_intron = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Intron"),]

#Create a dataframe which has only the negative control gene probes
nimble_nc = data_nimble_hk[which(data_nimble_hk[,"Probe_Type"] == "Control-Negative"),]
affy_nc = data_affy_hk[which(data_affy_hk[,"Probe_Type"] == "Control-Negative"),]


#Create a scatter plot of housekeeping exon and intron probe intensities and contrast with negative-control probes
#Use MIP_GH_A as an example

draw_hk_nimble = function(exp){
  title = paste ("NimbleGen Intensity of HK Exon and Intron Probes (", exp, ")")
  plot(nimble_nc[,"Probe_Tm"], log2(nimble_nc[,exp]), ylim=c(6,18), col="red", ylab="RAW Log2 Probe Intensity",
       xlab="Probe Tm", main=title, pch=1)
  points(nimble_hk_intron[,"Probe_Tm"], log2(nimble_hk_intron[,exp]), ylim=c(6,18), col="green", pch=1)
  points(nimble_hk_exon[,"Probe_Tm"], log2(nimble_hk_exon[,exp]), ylim=c(6,18), col="blue", pch=1)
  legend(locator(n=1,type="n"), c("Exonic","Intronic","Random Sequences"), pch=c(1,1,1), col=c("blue","green","red"))
  return(1)
}

draw_hk_affy = function(exp){
  title = paste ("Affy Intensity of HK Exon and Intron Probes (", exp, ")")
  plot(affy_nc[,"Probe_Tm"], log2(affy_nc[,exp]), ylim=c(4,16), col="red", ylab="RAW Log2 Probe Intensity",
       xlab="Probe Tm", main=title, pch=1)
  points(affy_hk_intron[,"Probe_Tm"], log2(affy_hk_intron[,exp]), ylim=c(4,16), col="green", pch=1)
  points(affy_hk_exon[,"Probe_Tm"], log2(affy_hk_exon[,exp]), ylim=c(4,16), col="blue", pch=1)
  legend(locator(n=1,type="n"), c("Exonic","Intronic","Random Sequences"), pch=c(1,1,1), col=c("blue","green","red"))
  return(1)
}

draw_hk_nimble("MIP_EF_B")
draw_hk_affy("MIP_EF_B")
draw_hk_nimble("FUR_EF_B")
draw_hk_affy("FUR_EF_B")
draw_hk_nimble("MIP_GH_A")
draw_hk_affy("MIP_GH_A")
draw_hk_nimble("FUR_GH_A")
draw_hk_affy("FUR_GH_A")
draw_hk_nimble("MIP_EF_A")
draw_hk_affy("MIP_C")
draw_hk_nimble("FUR_EF_A")
draw_hk_affy("FUR_C")

#Create box plots of the same comparisons
box_hk_nimble = function(exp){
  title = paste ("NimbleGen Intensity of HK Exon and Intron Probes (", exp, ")")
  probe_list = list(log2(nimble_hk_exon[,exp]), log2(nimble_hk_intron[,exp]), log2(nimble_nc[,exp]))
  names(probe_list) = c("Exonic","Intronic","Random")
  boxplot(probe_list, main=title, ylab="Log2 Intensity", ylim=c(4,16), col=c("blue","green","red"))
  return(1)
}
box_hk_affy = function(exp){
  title = paste ("Affy Intensity of HK Exon and Intron Probes (", exp, ")")
  probe_list = list(log2(affy_hk_exon[,exp]), log2(affy_hk_intron[,exp]), log2(affy_nc[,exp]))
  names(probe_list) = c("Exonic","Intronic","Random")
  boxplot(probe_list, main=title, ylab="Log2 Intensity", ylim=c(4,16), col=c("blue","green","red"))
  return(1)
}

par(mfrow=c(2,1))
 box_hk_nimble("MIP_EF_B")
 box_hk_affy("MIP_EF_B")
par(mfrow=c(2,1))
 box_hk_nimble("FUR_EF_B")
 box_hk_affy("FUR_EF_B")
par(mfrow=c(2,1))
 box_hk_nimble("MIP_GH_A")
 box_hk_affy("MIP_GH_A")
par(mfrow=c(2,1))
 box_hk_nimble("FUR_GH_A")
 box_hk_affy("FUR_GH_A")
par(mfrow=c(2,1))
 box_hk_nimble("MIP_EF_A")
 box_hk_affy("MIP_C")
par(mfrow=c(2,1))
 box_hk_nimble("FUR_EF_A")
 box_hk_affy("FUR_C")


#Calculate the median of all HK Exon probe intensities for all hybes
nimble_exon_medians = as.numeric(lapply (list((nimble_hk_exon[,"MIP_EF_A"]), (nimble_hk_exon[,"MIP_EF_B"]), (nimble_hk_exon[,"MIP_GH_A"]), 
                                   (nimble_hk_exon[,"FUR_EF_A"]), (nimble_hk_exon[,"FUR_EF_B"]), (nimble_hk_exon[,"FUR_GH_A"])), median))

affy_exon_medians = as.numeric(lapply (list((affy_hk_exon[,"MIP_C"]), (affy_hk_exon[,"MIP_EF_B"]), (affy_hk_exon[,"MIP_GH_A"]), 
                                 (affy_hk_exon[,"FUR_C"]), (affy_hk_exon[,"FUR_EF_B"]), (affy_hk_exon[,"FUR_GH_A"])), median))

#Calculate the median of all HK Intron probe intensities for all hybes
nimble_intron_medians = as.numeric(lapply (list((nimble_hk_intron[,"MIP_EF_A"]), (nimble_hk_intron[,"MIP_EF_B"]), (nimble_hk_intron[,"MIP_GH_A"]), 
                                     (nimble_hk_intron[,"FUR_EF_A"]), (nimble_hk_intron[,"FUR_EF_B"]), (nimble_hk_intron[,"FUR_GH_A"])), median))

affy_intron_medians = as.numeric(lapply (list((affy_hk_intron[,"MIP_C"]), (affy_hk_intron[,"MIP_EF_B"]), (affy_hk_intron[,"MIP_GH_A"]), 
                                   (affy_hk_intron[,"FUR_C"]), (affy_hk_intron[,"FUR_EF_B"]), (affy_hk_intron[,"FUR_GH_A"])), median))

nimble_ei = (nimble_exon_medians / nimble_intron_medians)
affy_ei = (affy_exon_medians / affy_intron_medians)

mean(nimble_ei)
sqrt(var(nimble_ei))
mean(affy_ei)
sqrt(var(affy_ei))

#Sample sizes are too small to prove normality
library(nortest)
ad.test(nimble_ei)
ad.test(affy_ei)

#Are these signal-to-noise ratios significantly different between ALEXA and AFFY??
#t.test(x=nimble_ei, y = affy_ei, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
wilcox.test(x=nimble_ei, y=affy_ei, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, 
	     conf.int = FALSE, conf.level = 0.95)

#Calculate sensitivity/specificity across the entire set of housekeeping genes

#First create a data frame to store results for each experiment/platform
cutoffs = seq(1, 2^16, 10)

calculate_roc_scores = function(exp, probe_data){

  x = array(0, c(length(cutoffs),2), dimnames=list(cutoffs, c("sensitivity","specificity")))
  result = as.data.frame(x) 

  #Get the exon and intron probes for this experiment 
  exon_probes = probe_data[which(probe_data[,"Probe_Type"] == "Exon"),]
  intron_probes = probe_data[which(probe_data[,"Probe_Type"] == "Intron"),]

  true_pos = length(exon_probes[,exp])
  true_neg = length(intron_probes[,exp])

  #Iterate through posible cutoff intensity values
  for(i in 1:length(cutoffs)){

    cutoff = cutoffs[i]

    obs_pos = length(which(exon_probes[,exp] >= cutoff))
    obs_neg = length(which(intron_probes[,exp] <= cutoff))

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

nimble_roc_mip_ef_b = calculate_roc_scores("MIP_EF_B", data_nimble_hk)
 nimble_auc_mip_ef_b = trapz(nimble_roc_mip_ef_b[,"specificity"],nimble_roc_mip_ef_b[,"sensitivity"])
nimble_roc_mip_gh_a = calculate_roc_scores("MIP_GH_A", data_nimble_hk)
 nimble_auc_mip_gh_a = trapz(nimble_roc_mip_gh_a[,"specificity"],nimble_roc_mip_gh_a[,"sensitivity"])
nimble_roc_mip_ef_a = calculate_roc_scores("MIP_EF_A", data_nimble_hk) 
 nimble_auc_mip_ef_a = trapz(nimble_roc_mip_ef_a[,"specificity"],nimble_roc_mip_ef_a[,"sensitivity"])

affy_roc_mip_ef_b = calculate_roc_scores("MIP_EF_B", data_affy_hk)
 affy_auc_mip_ef_b = trapz(affy_roc_mip_ef_b[,"specificity"],affy_roc_mip_ef_b[,"sensitivity"])
affy_roc_mip_gh_a = calculate_roc_scores("MIP_GH_A", data_affy_hk)
 affy_auc_mip_gh_a = trapz(affy_roc_mip_gh_a[,"specificity"],affy_roc_mip_gh_a[,"sensitivity"])
affy_roc_mip_c = calculate_roc_scores("MIP_C", data_affy_hk) 
 affy_auc_mip_c = trapz(affy_roc_mip_c[,"specificity"],affy_roc_mip_c[,"sensitivity"])

nimble_roc_fur_ef_b = calculate_roc_scores("FUR_EF_B", data_nimble_hk)
 nimble_auc_fur_ef_b = trapz(nimble_roc_fur_ef_b[,"specificity"],nimble_roc_fur_ef_b[,"sensitivity"])
nimble_roc_fur_gh_a = calculate_roc_scores("FUR_GH_A", data_nimble_hk)
 nimble_auc_fur_gh_a = trapz(nimble_roc_fur_gh_a[,"specificity"],nimble_roc_fur_gh_a[,"sensitivity"])
nimble_roc_fur_ef_a = calculate_roc_scores("FUR_EF_A", data_nimble_hk) 
 nimble_auc_fur_ef_a = trapz(nimble_roc_fur_ef_a[,"specificity"],nimble_roc_fur_ef_a[,"sensitivity"])

affy_roc_fur_ef_b = calculate_roc_scores("FUR_EF_B", data_affy_hk)
 affy_auc_fur_ef_b = trapz(affy_roc_fur_ef_b[,"specificity"],affy_roc_fur_ef_b[,"sensitivity"])
affy_roc_fur_gh_a = calculate_roc_scores("FUR_GH_A", data_affy_hk)
 affy_auc_fur_gh_a = trapz(affy_roc_fur_gh_a[,"specificity"],affy_roc_fur_gh_a[,"sensitivity"])
affy_roc_fur_c = calculate_roc_scores("FUR_C", data_affy_hk) 
 affy_auc_fur_c = trapz(affy_roc_fur_c[,"specificity"],affy_roc_fur_c[,"sensitivity"])

#Create ROC graphs with the data generated above

#MIP Data
#Create names for each experiment that include the AUC value
name_n_EF_A = paste("Nimble_EF_A", "(AUC =", round(nimble_auc_mip_ef_a, digit=4), ")")
name_n_EF_B = paste("Nimble_EF_B", "(AUC =", round(nimble_auc_mip_ef_b, digit=4), ")")
name_n_GH_A = paste("Nimble_GH_A", "(AUC =", round(nimble_auc_mip_gh_a, digit=4), ")")
name_a_C = paste("Affy_C", "(AUC =", round(affy_auc_mip_c, digit=4), ")")
name_a_EF_B = paste("Affy_EF_B", "(AUC =", round(affy_auc_mip_ef_b, digit=4), ")")
name_a_GH_A = paste("Affy_GH_A", "(AUC =", round(affy_auc_mip_gh_a, digit=4), ")")

plot(1-(nimble_roc_mip_ef_a[,"specificity"]), nimble_roc_mip_ef_a[,"sensitivity"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affy vs. NimbleGen - 3 Hybes Each (Sensitive data)", lty=1, lwd=2)
lines(1-(nimble_roc_mip_ef_b[,"specificity"]), nimble_roc_mip_ef_b[,"sensitivity"], col="blue", type="l", lty=2, lwd=2)
lines(1-(nimble_roc_mip_gh_a[,"specificity"]), nimble_roc_mip_gh_a[,"sensitivity"], col="blue", type="l", lty=3, lwd=2)
lines(1-(affy_roc_mip_c[,"specificity"]), affy_roc_mip_c[,"sensitivity"], col="red", type="l", lty=1, lwd=2)
lines(1-(affy_roc_mip_ef_b[,"specificity"]), affy_roc_mip_ef_b[,"sensitivity"], col="red", type="l", lty=2, lwd=2)
lines(1-(affy_roc_mip_gh_a[,"specificity"]), affy_roc_mip_gh_a[,"sensitivity"], col="red", type="l", lty=3, lwd=2)
legend(locator(n=1,type="n"), c(name_n_EF_A,name_n_EF_B,name_n_GH_A,name_a_C,name_a_EF_B,name_a_GH_A), lty=c(1,2,3,1,2,3), col=c("blue","blue","blue","red","red","red"), lwd=c(2,2,2,2,2,2))

#5FUR Data
#Create names for each experiment that include the AUC value
name_n_EF_A = paste("Nimble_EF_A", "(AUC =", round(nimble_auc_fur_ef_a, digit=4), ")")
name_n_EF_B = paste("Nimble_EF_B", "(AUC =", round(nimble_auc_fur_ef_b, digit=4), ")")
name_n_GH_A = paste("Nimble_GH_A", "(AUC =", round(nimble_auc_fur_gh_a, digit=4), ")")
name_a_C = paste("Affy_C", "(AUC =", round(affy_auc_fur_c, digit=4), ")")
name_a_EF_B = paste("Affy_EF_B", "(AUC =", round(affy_auc_fur_ef_b, digit=4), ")")
name_a_GH_A = paste("Affy_GH_A", "(AUC =", round(affy_auc_fur_gh_a, digit=4), ")")

plot(1-(nimble_roc_fur_ef_a[,"specificity"]), nimble_roc_fur_ef_a[,"sensitivity"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affy vs. NimbleGen - 3 Hybes Each (Resistant data)", lty=1, lwd=2)
lines(1-(nimble_roc_fur_ef_b[,"specificity"]), nimble_roc_fur_ef_b[,"sensitivity"], col="blue", type="l", lty=2, lwd=2)
lines(1-(nimble_roc_fur_gh_a[,"specificity"]), nimble_roc_fur_gh_a[,"sensitivity"], col="blue", type="l", lty=3, lwd=2)

lines(1-(affy_roc_fur_c[,"specificity"]), affy_roc_fur_c[,"sensitivity"], col="red", type="l", lty=1, lwd=2)
lines(1-(affy_roc_fur_ef_b[,"specificity"]), affy_roc_fur_ef_b[,"sensitivity"], col="red", type="l", lty=2, lwd=2)
lines(1-(affy_roc_fur_gh_a[,"specificity"]), affy_roc_fur_gh_a[,"sensitivity"], col="red", type="l", lty=3, lwd=2)
legend(locator(n=1,type="n"), c(name_n_EF_A,name_n_EF_B,name_n_GH_A,name_a_C,name_a_EF_B,name_a_GH_A), lty=c(1,2,3,1,2,3), col=c("blue","blue","blue","red","red","red"), lwd=c(2,2,2,2,2,2))


#Combine all sensitivity/specificity scores (median at each cutoff) for the two platforms
#For each cutoff, also calculate the SD values
#Plot the median ROC curves and corresponding SD curves for both platforms

nimble_roc_all = cbind(nimble_roc_mip_ef_b, nimble_roc_mip_gh_a, nimble_roc_mip_ef_a, nimble_roc_fur_ef_b, nimble_roc_fur_gh_a, nimble_roc_fur_ef_a)
names(nimble_roc_all) = c("sen_mip_ef_b", "spec_mip_ef_b", "sen_mip_gh_a", "spec_mip_gh_a", "sen_mip_ef_a", "spec_mip_ef_a", "sen_fur_ef_b", "spec_fur_ef_b", "sen_fur_gh_a", "spec_fur_gh_a", "sen_fur_ef_a", "spec_fur_ef_a")

affy_roc_all = cbind(affy_roc_mip_ef_b, affy_roc_mip_gh_a, affy_roc_mip_c, affy_roc_fur_ef_b, affy_roc_fur_gh_a, affy_roc_fur_c)
names(affy_roc_all) = c("sen_mip_ef_b", "spec_mip_ef_b", "sen_mip_gh_a", "spec_mip_gh_a", "sen_mip_ef_a", "spec_mip_ef_a", "sen_fur_ef_b", "spec_fur_ef_b", "sen_fur_gh_a", "spec_fur_gh_a", "sen_fur_ef_a", "spec_fur_ef_a")

calculate_mean_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = mean(sen)
  return(mean_sen)
}
nimble_roc_all[,"mean_sens"] = apply(nimble_roc_all, 1, calculate_mean_sens)
affy_roc_all[,"mean_sens"] = apply(affy_roc_all, 1, calculate_mean_sens)

calculate_mean_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = mean(spec)
  return(mean_spec)
}
nimble_roc_all[,"mean_spec"] = apply(nimble_roc_all, 1, calculate_mean_spec)
affy_roc_all[,"mean_spec"] = apply(affy_roc_all, 1, calculate_mean_spec)

nimble_auc_all = trapz(nimble_roc_all[,"mean_spec"],nimble_roc_all[,"mean_sens"])
affy_auc_all = trapz(affy_roc_all[,"mean_spec"],affy_roc_all[,"mean_sens"])

calculate_lower_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  lower_sen = mean_sen - sd_sen
  return(lower_sen)
}
nimble_roc_all[,"lower_sens"] = apply(nimble_roc_all, 1, calculate_lower_sens)
affy_roc_all[,"lower_sens"] = apply(affy_roc_all, 1, calculate_lower_sens)

calculate_upper_sens = function(roc_data){
  sen = as.numeric(c(roc_data["sen_mip_ef_b"], roc_data["sen_mip_gh_a"], roc_data["sen_mip_ef_a"], roc_data["sen_fur_ef_b"], roc_data["sen_fur_gh_a"], roc_data["sen_fur_ef_a"]))
  mean_sen = roc_data["mean_sens"]
  sd_sen = sd(sen)
  upper_sen = mean_sen + sd_sen
  return(upper_sen)
}
nimble_roc_all[,"upper_sens"] = apply(nimble_roc_all, 1, calculate_upper_sens)
affy_roc_all[,"upper_sens"] = apply(affy_roc_all, 1, calculate_upper_sens)

calculate_lower_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  lower_spec = mean_spec - sd_spec
  return(lower_spec)
}
nimble_roc_all[,"lower_spec"] = apply(nimble_roc_all, 1, calculate_lower_spec)
affy_roc_all[,"lower_spec"] = apply(affy_roc_all, 1, calculate_lower_spec)

calculate_upper_spec = function(roc_data){
  spec = as.numeric(c(roc_data["spec_mip_ef_b"], roc_data["spec_mip_gh_a"], roc_data["spec_mip_ef_a"], roc_data["spec_fur_ef_b"], roc_data["spec_fur_gh_a"], roc_data["spec_fur_ef_a"]))
  mean_spec = roc_data["mean_spec"]
  sd_spec = sd(spec)
  upper_spec = mean_spec + sd_spec
  return(upper_spec)
}
nimble_roc_all[,"upper_spec"] = apply(nimble_roc_all, 1, calculate_upper_spec)
affy_roc_all[,"upper_spec"] = apply(affy_roc_all, 1, calculate_upper_spec)

#Determine the max combined specificity and sensitivity achieved by each platform
nimble_roc_all[,"SS_combined"] = nimble_roc_all[,"mean_sens"] + nimble_roc_all[,"mean_spec"]
max_combo = max(nimble_roc_all[,"SS_combined"])
nimble_roc_all[which(nimble_roc_all[,"SS_combined"]==max_combo),]

affy_roc_all[,"SS_combined"] = affy_roc_all[,"mean_sens"] + affy_roc_all[,"mean_spec"]
max_combo = max(affy_roc_all[,"SS_combined"])
affy_roc_all[which(affy_roc_all[,"SS_combined"]==max_combo),]


#Create names for each experiment that include the AUC value
name_alexa = paste("ALEXA", "(AUC =", round(nimble_auc_all, digit=3), ")")
name_affy = paste("Affymetrix", "(AUC =", round(affy_auc_all, digit=3), ")")

plot(1-(nimble_roc_all[,"mean_spec"]), nimble_roc_all[,"mean_sens"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affymetrix vs. ALEXA", lty=1, lwd=2)
lines(1-(nimble_roc_all[,"lower_spec"]), nimble_roc_all[,"lower_sens"], col="blue", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
lines(1-(nimble_roc_all[,"upper_spec"]), nimble_roc_all[,"upper_sens"], col="blue", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))

lines(1-(affy_roc_all[,"mean_spec"]), affy_roc_all[,"mean_sens"], col="red", type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,1))
lines(1-(affy_roc_all[,"lower_spec"]), affy_roc_all[,"lower_sens"], col="red", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))
lines(1-(affy_roc_all[,"upper_spec"]), affy_roc_all[,"upper_sens"], col="red", type="l", lty=2, lwd=1, ylim=c(0,1), xlim=c(0,1))

legend(locator(n=1,type="n"), c(name_alexa,"S.D.",name_affy,"S.D."), lty=c(1,2,1,2), col=c("blue","blue","red","red"), lwd=c(2,1,2,1))
locator(n=1,type="n")





#NORMALIZED DATA ANALYSIS
library(gcrma)

#Quantiles normalize WITHIN each condition (3 sensitive normalized together, 3 resistant normalized together)
#Then take mean WITHIN each condition for each probe and plot these for both conditions for both platforms

#Import Affy housekeeping gene data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel"
dir(datadir)
setwd(datadir)

data_affy_hk = read.table("HuEX-1_001-006_RawProbe_ControlsOnly.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(2,6), na.strings='na')
master_hk_gene_list = unique(data_affy_hk[,"AlexaGene_ID"])

summary(data_affy_hk[,"MIP_C"])
summary(data_affy_hk[,"MIP_GH_A"])
summary(data_affy_hk[,"MIP_EF_B"])

#Normalize data WITHIN sample groups (i.e. all sensitive together, all resistant together
#MIP
names = c("MIP_C","MIP_GH_A","MIP_EF_B")
mip_data_combined = data.frame(data_affy_hk[,"MIP_C"],data_affy_hk[,"MIP_GH_A"],data_affy_hk[,"MIP_EF_B"])
names(mip_data_combined) = names
x = as.matrix(mip_data_combined)

mip_data_quant_norm = normalize.quantiles(x)
mip_data_quant_norm = as.data.frame(mip_data_quant_norm)
dimnames(mip_data_quant_norm)[[2]] = names

data_affy_hk[,"MIP_C"] = mip_data_quant_norm[,"MIP_C"]
data_affy_hk[,"MIP_GH_A"] = mip_data_quant_norm[,"MIP_GH_A"]
data_affy_hk[,"MIP_EF_B"] = mip_data_quant_norm[,"MIP_EF_B"]
summary(data_affy_hk[,"MIP_C"])
summary(data_affy_hk[,"MIP_GH_A"])
summary(data_affy_hk[,"MIP_EF_B"])

data_affy_hk[,"MIP_MEAN"] = apply(data_affy_hk[,c("MIP_C","MIP_GH_A","MIP_EF_B")], 1, mean)

#FUR
summary(data_affy_hk[,"FUR_C"])
summary(data_affy_hk[,"FUR_GH_A"])
summary(data_affy_hk[,"FUR_EF_B"])

names = c("FUR_C","FUR_GH_A","FUR_EF_B")
fur_data_combined = data.frame(data_affy_hk[,"FUR_C"],data_affy_hk[,"FUR_GH_A"],data_affy_hk[,"FUR_EF_B"])
names(fur_data_combined) = names
x = as.matrix(fur_data_combined)

fur_data_quant_norm = normalize.quantiles(x)
fur_data_quant_norm = as.data.frame(fur_data_quant_norm)
dimnames(fur_data_quant_norm)[[2]] = names

data_affy_hk[,"FUR_C"] = fur_data_quant_norm[,"FUR_C"]
data_affy_hk[,"FUR_GH_A"] = fur_data_quant_norm[,"FUR_GH_A"]
data_affy_hk[,"FUR_EF_B"] = fur_data_quant_norm[,"FUR_EF_B"]
summary(data_affy_hk[,"FUR_C"])
summary(data_affy_hk[,"FUR_GH_A"])
summary(data_affy_hk[,"FUR_EF_B"])

data_affy_hk[,"FUR_MEAN"] = apply(data_affy_hk[,c("FUR_C","FUR_GH_A","FUR_EF_B")], 1, mean)

#Import NimbleGen data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted"
dir(datadir)
setwd(datadir)

#After NMPP smoothing!!
data_nimble_complete = read.table("All_hybes_withProbeInfo_smoothed.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,14,17), na.strings='na')

#Create a dataframe which has only the overlaping housekeeping gene probes
data_nimble_hk = data_nimble_complete[which(data_nimble_complete[,"AlexaGene_ID"] %in% master_hk_gene_list),]

summary(data_nimble_hk[,"MIP_EF_A"])
summary(data_nimble_hk[,"MIP_EF_B"])
summary(data_nimble_hk[,"MIP_GH_A"])

#Normalize data WITHIN sample groups (i.e. all sensitive together, all resistant together
#MIP
names = c("MIP_EF_A","MIP_EF_B","MIP_GH_A")
mip_data_combined = data.frame(data_nimble_hk[,"MIP_EF_A"],data_nimble_hk[,"MIP_EF_B"],data_nimble_hk[,"MIP_GH_A"])
names(mip_data_combined) = names
x = as.matrix(mip_data_combined)

mip_data_quant_norm = normalize.quantiles(x)
mip_data_quant_norm = as.data.frame(mip_data_quant_norm)
dimnames(mip_data_quant_norm)[[2]] = names

data_nimble_hk[,"MIP_EF_A"] = mip_data_quant_norm[,"MIP_EF_A"]
data_nimble_hk[,"MIP_EF_B"] = mip_data_quant_norm[,"MIP_EF_B"]
data_nimble_hk[,"MIP_GH_A"] = mip_data_quant_norm[,"MIP_GH_A"]

summary(data_nimble_hk[,"MIP_EF_A"])
summary(data_nimble_hk[,"MIP_EF_B"])
summary(data_nimble_hk[,"MIP_GH_A"])

data_nimble_hk[,"MIP_MEAN"] = apply(data_nimble_hk[,c("MIP_EF_A","MIP_EF_B","MIP_GH_A")], 1, mean)

#FUR
summary(data_nimble_hk[,"FUR_EF_A"])
summary(data_nimble_hk[,"FUR_EF_B"])
summary(data_nimble_hk[,"FUR_GH_A"])

names = c("FUR_EF_A","FUR_EF_B","FUR_GH_A")
fur_data_combined = data.frame(data_nimble_hk[,"FUR_EF_A"],data_nimble_hk[,"FUR_EF_A"],data_nimble_hk[,"FUR_GH_A"])
names(fur_data_combined) = names
x = as.matrix(fur_data_combined)

fur_data_quant_norm = normalize.quantiles(x)
fur_data_quant_norm = as.data.frame(fur_data_quant_norm)
dimnames(fur_data_quant_norm)[[2]] = names

data_nimble_hk[,"FUR_EF_A"] = fur_data_quant_norm[,"FUR_EF_A"]
data_nimble_hk[,"FUR_EF_B"] = fur_data_quant_norm[,"FUR_EF_B"]
data_nimble_hk[,"FUR_GH_A"] = fur_data_quant_norm[,"FUR_GH_A"]

summary(data_nimble_hk[,"FUR_EF_A"])
summary(data_nimble_hk[,"FUR_EF_B"])
summary(data_nimble_hk[,"FUR_GH_A"])

data_nimble_hk[,"FUR_MEAN"] = apply(data_nimble_hk[,c("FUR_EF_A","FUR_EF_B","FUR_GH_A")], 1, mean)

#CHANGE TO AN OUTPUT DIRECTORY
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/Figures"
dir(datadir)
setwd(datadir)

cutoffs = seq(1, 2^16, 10)

calculate_roc_scores = function(exp, probe_data){

  x = array(0, c(length(cutoffs),2), dimnames=list(cutoffs, c("sensitivity","specificity")))
  result = as.data.frame(x) 

  #Get the exon and intron probes for this experiment 
  exon_probes = probe_data[which(probe_data[,"Probe_Type"] == "Exon"),]
  intron_probes = probe_data[which(probe_data[,"Probe_Type"] == "Intron"),]

  true_pos = length(exon_probes[,exp])
  true_neg = length(intron_probes[,exp])

  #Iterate through posible cutoff intensity values
  for(i in 1:length(cutoffs)){

    cutoff = cutoffs[i]

    obs_pos = length(which(exon_probes[,exp] >= cutoff))
    obs_neg = length(which(intron_probes[,exp] <= cutoff))

    result[i,"sensitivity"] = obs_pos/true_pos
    result[i,"specificity"] = obs_neg/true_neg

  }
  return(result)
}

#Also Calculate the Area Under the Curve for each line
library(caTools)

nimble_roc_mip = calculate_roc_scores("MIP_MEAN", data_nimble_hk)
 nimble_auc_mip = trapz(nimble_roc_mip[,"specificity"],nimble_roc_mip[,"sensitivity"])

affy_roc_mip = calculate_roc_scores("MIP_MEAN", data_affy_hk)
 affy_auc_mip = trapz(affy_roc_mip[,"specificity"],affy_roc_mip[,"sensitivity"])

nimble_roc_fur = calculate_roc_scores("FUR_MEAN", data_nimble_hk)
 nimble_auc_fur = trapz(nimble_roc_fur[,"specificity"],nimble_roc_fur[,"sensitivity"])

affy_roc_fur = calculate_roc_scores("FUR_MEAN", data_affy_hk)
 affy_auc_fur = trapz(affy_roc_fur[,"specificity"],affy_roc_fur[,"sensitivity"])

#Create ROC graphs with the data generated above

#Mean Data
#Create names for each experiment that include the AUC value
name_n_MIP = paste("Nimble_MIP", "(AUC =", round(nimble_auc_mip, digit=4), ")")
name_n_FUR = paste("Nimble_FUR", "(AUC =", round(nimble_auc_fur, digit=4), ")")
name_a_MIP = paste("Affy_MIP", "(AUC =", round(affy_auc_mip, digit=4), ")")
name_a_FUR = paste("Affy_FUR", "(AUC =", round(affy_auc_fur, digit=4), ")")

plot(1-(nimble_roc_mip[,"specificity"]), nimble_roc_mip[,"sensitivity"], ylim=c(0,1), xlim=c(0,1), col="blue", type="l",
     xlab="1 - Specificity", ylab="Sensitivity", main="ROC Curves - Affy vs. NimbleGen - Mean Data", lty=1, lwd=2)
lines(1-(nimble_roc_fur[,"specificity"]), nimble_roc_fur[,"sensitivity"], col="blue", type="l", lty=2, lwd=2)
lines(1-(affy_roc_mip[,"specificity"]), affy_roc_mip[,"sensitivity"], col="red", type="l", lty=1, lwd=2)
lines(1-(affy_roc_fur[,"specificity"]), affy_roc_fur[,"sensitivity"], col="red", type="l", lty=2, lwd=2)
legend(locator(n=1,type="n"), c(name_n_MIP,name_n_FUR,name_a_MIP,name_a_FUR), lty=c(1,2,1,2), col=c("blue","blue","red","red"), lwd=c(2,2,2,2))














