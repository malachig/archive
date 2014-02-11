#Written by Malachi Griffith
#Identify a list of candidate significant events (exons, junctions, boundaries)
#This will involve:
#1.) first filtering events to eliminate those that have poor evidence for expression in either condition
#2.) then filtering to eliminate those with only a minor change in expression between conditions
#3.) then calculating a p-value for this observed change between conditions
#4.) then correcting this p-value for multiple testing
#5.) then filtering out those events with a p-value > 0.05
#6.) and finally printing out the final list of candidate events

#Load neccessary libraries
library(multtest)
library(nortest)

#1.) Import DATA
#Load exon expression values for each probeset
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)

exon_exp_data = read.table(file="Standard_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])

#2.) Testing the assumption of normality
#Before deciding on whether to use t-tests, examine the assumption of normality
#Use QQ-Plots and statistical tests of normality such as the following from the 'nortest' package:
#Lilliefors (Kolmogorov-Smirnov) test for the composite hypothesis of normality, see e.g. Thode (2002, Sec. 5.1.1)
#Cramer-von Mises test for the composite hypothesis of normality, see e.g. Thode (2002, Sec. 5.1.3).
#Anderson-Darling test for the composite hypothesis of normality, see e.g. Thode (2002, Sec. 5.1.4).
#This test should be conducted on each of the populations being compared (sensitive and resistant) for each 
#test conducted (each exon)??
#If the p-value is sufficiently small it suggests that the populations being compared are NOT normally distributed
#Examples of normality tests follow

x = rnorm(100, mean=5, sd=3)
y = runif(100, min=2, max=4)
qqnorm(x)
qqline(x, datax = FALSE, col="red", lwd=2)
ad.test(x)
cvm.test(x)
lillie.test(x)
sf.test(x)

qqnorm(y)
qqline(y, datax = FALSE, col="red", lwd=2)
ad.test(y)
cvm.test(y)
lillie.test(y)
sf.test(y)

#For a sample size similar to those I will be testing below (n=9 say) how often does this test correctly identify a sample
#population as normal if it was randomly drawn from a random normal population with similar characteristics to my expression data
x = rnorm(9, mean=7.66, sd=2.80)
y = runif(9, min=4.0000, max=15.9313)
ad.test(x)
ad.test(y)

#Problem.  With sample sizes as small as n=9, the AD test will not indicate a lack of normality simply because it lacks power
#It also fails to detect the non-normality of samples drawn from a UNIFORM distribution.
#For this reason it is not safe to use a t-test or other parametric test because we can't convince ourselves that the samples are normal

#3.) Before identifying significant events first filter to eliminate those that have poor evidence for expression in either condition

#What is a good expression cutoff?  Base this on the mean expression values of negative control probes.
nc_exon_exp_data = exon_exp_data[which(exon_exp_data[,"Probe_Type"] == "Control-Negative"),]
length(nc_exon_exp_data[,1])
quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=seq(0, 1, 0.05))
quantile(nc_exon_exp_data[,"fur_LOG2_I_U"], probs=seq(0, 1, 0.05))
expression_cutoff = quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=0.975)

#Now identify the probesets that meet the cutoff
expressed_probesets = exon_exp_data[which(exon_exp_data[,"mip_LOG2_I_U"] >= expression_cutoff | exon_exp_data[,"fur_LOG2_I_U"] >= expression_cutoff),]
length(expressed_probesets[,1])

#4.) Eliminate events that are not variant between the conditions (and therefore not of interest 
#Now filter out those without an absolute DE value of > 1 (2-fold)

de_cutoff = log2(2)
de_probesets = exon_exp_data[which(abs(exon_exp_data[,"mip_v_fur_DE_U"]) > de_cutoff),]
length(de_probesets[,1])
de_probesets_list = de_probesets[,c("ProbeSet_ID","Probe_Type")]

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_excluded"
dir(datadir)
setwd(datadir)
exon_data_mtp = read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(exon_data_mtp[,1])

#Load H19 data seperately
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_only"
dir(datadir)
setwd(datadir)
exon_data_mtp_H19 = read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(exon_data_mtp_H19[,1])

#5.) Now filter out probesets that appear to be invariant or 'non-responsive' across conditions
#Try calculating the 'Coefficient of Variation (CV)' as a means of identifying non-responding probes to be filtered out
#This is a metric used by D-CHIP for this purpose.  CV = SD/mean
#It should not be calculated on log2 values because the log transformation results in stabilizing of variance  
#NOTE: CV doesnt seem to work well for my situation.  Even events with large changes have very low CV because the means are so high

#calculate_sd = function(dataset){
#  start_data_column = 6
#  end_data_column = length(dataset)-1
#  variable_data=as.numeric((dataset[start_data_column:end_data_column]))
#  mean_value = mean(as.numeric(variable_data), na.rm = TRUE)
#  sd_value = sd(as.numeric(variable_data), na.rm = TRUE)
#  test_value = sd_value/mean_value
#  observations = length(which(!is.na(variable_data)))
#  return(sd_value)
#}
#sd_values = apply(exon_data_mtp[,], 1, calculate_sd) 
#temp_data = exon_data_mtp[,c("ProbeSet_ID","AlexaGene_ID")]
#temp_data[,"SD"] = sd_values
#sd_values_H19 = apply(exon_data_mtp_H19[,], 1, calculate_sd)
#temp_data_H19 = exon_data_mtp_H19[,c("ProbeSet_ID","AlexaGene_ID")]
#temp_data_H19[,"SD"] = sd_values_H19
#sd_cutoff = 1
#variant_probesets = temp_data[which(temp_data[,"SD"] > sd_cutoff),]
#length(variant_probesets[,1])
#variant_probesets_H19 = temp_data_H19[which(temp_data_H19[,"SD"] > sd_cutoff),]
#length(variant_probesets_H19[,1])

#Using the lists of probesets determined above, grab the subset of probesets which passed the filters
exon_data_mtp_filtered_1 = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"] %in% expressed_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_1[,1])
#exon_data_mtp_filtered = exon_data_mtp_filtered_1[which(exon_data_mtp_filtered_1[,"ProbeSet_ID"] %in% variant_probesets[,"ProbeSet_ID"]),]
#length(exon_data_mtp_filtered[,1])
exon_data_mtp_filtered = exon_data_mtp_filtered_1[which(exon_data_mtp_filtered_1[,"ProbeSet_ID"] %in% de_probesets_list[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered[,1])

exon_data_mtp_H19_filtered_1 = exon_data_mtp_H19[which(exon_data_mtp_H19[,"ProbeSet_ID"] %in% expressed_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_1[,1])
#exon_data_mtp_H19_filtered = exon_data_mtp_H19_filtered_1[which(exon_data_mtp_H19_filtered_1[,"ProbeSet_ID"] %in% variant_probesets_H19[,"ProbeSet_ID"]),]
#length(exon_data_mtp_H19_filtered[,1])
exon_data_mtp_H19_filtered = exon_data_mtp_H19_filtered_1[which(exon_data_mtp_H19_filtered_1[,"ProbeSet_ID"] %in% de_probesets_list[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered[,1])

#6.) Now calculate a p-value for this observed change between conditions
# For this we need to load the MTP formatted data file

#Start with exon level log2 values for 3 replicates and 2 conditions (sensitive versus resistant)
#Each row contains all of the probe intensity observations for a probeset (which corresponds to an exon) for all replicates 
#Make sure the file is tab-delimited and missing values are set to NA 

#Use a simple Mann-Whitney test (aka Wilcoxen Rank Test)="t.twosamp.equalvar" for the reasons described above  

#A.) Create a function that will get a p-value for a single row of values consisting of two populations

#First 5 columns are not data

#wilcox.test
compute_pvalues = function(dataset){
  start_data_column = 6
  end_data_column = length(dataset)-1
  variable_data=dataset[start_data_column:end_data_column]
  half_data = length(variable_data)/2
  x=as.numeric(variable_data[1:half_data])
  y=as.numeric(variable_data[(half_data+1):(length(variable_data))])

  result = wilcox.test(x, y, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
  return(result$p.value)
}

#Try testing on different populations of probes
nc_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Control-Negative"),]
intron_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Intron"),]
exon_boundary_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Intron-Exon" | exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Intron"),]
exon_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon"),]
canonical_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Exon" & exon_data_mtp_filtered[,"Exons_Skipped"] == 0),]
skipping_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Exon" & exon_data_mtp_filtered[,"Exons_Skipped"] >= 1),]

#further breakdowns
s1_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Exon" & exon_data_mtp_filtered[,"Exons_Skipped"] == 1),]
s2_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Exon" & exon_data_mtp_filtered[,"Exons_Skipped"] == 2),]
s3_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Exon" & exon_data_mtp_filtered[,"Exons_Skipped"] >= 3),]
ei_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon-Intron"),]
ie_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Intron-Exon"),]

nc_pvals = apply(nc_data, 1, compute_pvalues)
intron_pvals = apply(intron_data, 1, compute_pvalues)
exon_boundary_pvals = apply(exon_boundary_data, 1, compute_pvalues)
exon_pvals = apply(exon_data, 1, compute_pvalues)
canonical_pvals = apply(canonical_data, 1, compute_pvalues)
skipping_pvals = apply(skipping_data, 1, compute_pvalues)
s1_pvals = apply(s1_data, 1, compute_pvalues)
s2_pvals = apply(s2_data, 1, compute_pvalues)
s3_pvals = apply(s3_data, 1, compute_pvalues)
ei_pvals = apply(ei_data, 1, compute_pvalues)
ie_pvals = apply(ie_data, 1, compute_pvalues)
H19_pvals = apply(exon_data_mtp_H19_filtered, 1, compute_pvalues)

#summarise the distribution of P-values for each of the groups - BEFORE MTP correction
summary (nc_pvals)
summary (intron_pvals)
summary (exon_boundary_pvals)
summary (exon_pvals)
summary (canonical_pvals)
summary (skipping_pvals)
summary (s1_pvals)
summary (s2_pvals)
summary (s3_pvals)
summary (ei_pvals)
summary (ie_pvals)
summary (H19_pvals)

#7.) P-value correction for multiple testing
#The mt.rawp2adjp Function
#Used to calculate adjusted p-values for simple multiple testing procedures from a vector of 
#unadjusted p-values.  7 procedures for getting adjusted p-values are used:
#Two groups (A) Strong control of family-wise Type I error rate (FWER)
#           (B) Strong control of false discovery rate (FDR)
#1.) Bonferoni FWER
#2.) Holm FWER
#3.) Hochberg FWER
#4.) Sidak-SS FWER
#5.) Sidak-SD FWER
#6.) Benjamini and Hochberg FDR
#7.) Benjamini and Yekutieli FDR

procs = c("Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
res_nc = mt.rawp2adjp(nc_pvals, procs)
res_intron = mt.rawp2adjp(intron_pvals, procs)
res_exon_boundary = mt.rawp2adjp(exon_boundary_pvals, procs)
res_exon = mt.rawp2adjp(exon_pvals, procs)
res_canonical = mt.rawp2adjp(canonical_pvals, procs)
res_skipping = mt.rawp2adjp(skipping_pvals, procs)
res_s1 = mt.rawp2adjp(s1_pvals, procs)
res_s2 = mt.rawp2adjp(s2_pvals, procs)
res_s3 = mt.rawp2adjp(s3_pvals, procs)
res_ei = mt.rawp2adjp(ei_pvals, procs)
res_ie = mt.rawp2adjp(ie_pvals, procs)
res_H19 = mt.rawp2adjp(H19_pvals, procs)

#Get the adjusted p-values and store in the original order
adjp_nc = res_nc$adjp[order(res_nc$index),"BH"]
adjp_intron = res_intron$adjp[order(res_intron$index),"BH"]
adjp_exon_boundary = res_exon_boundary$adjp[order(res_exon_boundary$index),"BH"]
adjp_exon = res_exon$adjp[order(res_exon$index),"BH"]
adjp_canonical = res_canonical$adjp[order(res_canonical$index),"BH"]
adjp_skipping = res_skipping$adjp[order(res_skipping$index),"BH"]
adjp_s1 = res_s1$adjp[order(res_s1$index),"BH"]
adjp_s2 = res_s2$adjp[order(res_s2$index),"BH"]
adjp_s3 = res_s3$adjp[order(res_s3$index),"BH"]
adjp_ei = res_ei$adjp[order(res_ei$index),"BH"]
adjp_ie = res_ie$adjp[order(res_ie$index),"BH"]
adjp_H19 = res_H19$adjp[order(res_H19$index),"BH"]

nc_data[,"adjp_BH"] = adjp_nc
intron_data[,"adjp_BH"] = adjp_intron
exon_boundary_data[,"adjp_BH"] = adjp_exon_boundary
exon_data[,"adjp_BH"] = adjp_exon
canonical_data[,"adjp_BH"] = adjp_canonical
skipping_data[,"adjp_BH"] = adjp_skipping
s1_data[,"adjp_BH"] = adjp_s1
s2_data[,"adjp_BH"] = adjp_s2
s3_data[,"adjp_BH"] = adjp_s3
ei_data[,"adjp_BH"] = adjp_ei
ie_data[,"adjp_BH"] = adjp_ie
exon_data_mtp_H19_filtered[,"adjp_BH"] = adjp_H19

#Summarize the adjusted p-values for each of the probe categories or sub-categories
summary(adjp_nc)
summary(adjp_intron)
summary(adjp_exon_boundary)
summary(adjp_exon)
summary(adjp_canonical)
summary(adjp_skipping)
summary(adjp_s1)
summary(adjp_s2)
summary(adjp_s3)
summary(adjp_ei)
summary(adjp_ie)
summary(adjp_H19)

length(which(adjp_nc < 0.05))
length(which(adjp_intron < 0.05))
length(which(adjp_exon_boundary < 0.05))
length(which(adjp_exon < 0.05))
length(which(adjp_canonical < 0.05))
length(which(adjp_skipping < 0.05))
length(which(adjp_s1 < 0.05))
length(which(adjp_s2 < 0.05))
length(which(adjp_s3 < 0.05))
length(which(adjp_ei < 0.05))
length(which(adjp_ie < 0.05))
length(which(adjp_H19 < 0.05))

total_alt_events = length(which(adjp_s1 < 0.05)) + length(which(adjp_s2 < 0.05)) + length(which(adjp_s3 < 0.05)) + length(which(adjp_ei < 0.05)) + length(which(adjp_ie < 0.05))
total_alt_events

#8.) Get the gene-level expression changes for these events
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)
gene_exp_data = read.table(file="Standard_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
names(gene_exp_data)[14] = "gene_mip_v_fur_DE_U"

#Gather actual mean DE values to be printed out with p-values
nc_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% nc_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
intron_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% intron_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
exon_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% exon_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
canonical_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% canonical_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
skipping_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% skipping_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
#s1_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% s1_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
#s2_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% s2_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
#s3_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% s3_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
ei_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% ei_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
ie_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% ie_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
exon_data_mtp_H19_filtered[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% exon_data_mtp_H19_filtered[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]

#Add gene-level data
nc_data = merge(x=nc_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
intron_data = merge(x=intron_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_data = merge(x=exon_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
canonical_data = merge(x=canonical_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
skipping_data = merge(x=skipping_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
#s1_data = merge(x=s1_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
#s2_data = merge(x=s2_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
#s3_data = merge(x=s3_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
ei_data = merge(x=ei_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
ie_data = merge(x=ie_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_data_mtp_H19_filtered = merge(x=exon_data_mtp_H19_filtered, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)


#9.) Print out the probeset IDs for these events
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/"
dir(datadir)
setwd(datadir)
write.table(nc_data[which(nc_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(intron_data[which(intron_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(exon_data[which(exon_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(canonical_data[which(canonical_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(skipping_data[which(skipping_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
#write.table(s1_data[which(s1_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
#write.table(s2_data[which(s2_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
#write.table(s3_data[which(s3_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(ei_data[which(ei_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(ie_data[which(ie_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(exon_data_mtp_H19_filtered[which(exon_data_mtp_H19_filtered[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))


#10.) Similar to the analysis for events above, identify GENES that are significantly DE between Sensitive and Resistant  

#10a) Load gene expression values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)
gene_exp_data = read.table(file="Standard_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(gene_exp_data[,1])

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/"
dir(datadir)
setwd(datadir)
gene_data_mtp = read.table(file="Standard_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(gene_data_mtp[,1])

#10b) Use NC probeset expression values to apply a simple expression cutoff 
#Identify genes that are not expressed in either condition
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)
exon_exp_data = read.table(file="Standard_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])
nc_exon_exp_data = exon_exp_data[which(exon_exp_data[,"Probe_Type"] == "Control-Negative"),]
length(nc_exon_exp_data[,1])
quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=seq(0, 1, 0.05))
quantile(nc_exon_exp_data[,"fur_LOG2_I_U"], probs=seq(0, 1, 0.05))
expression_cutoff = quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=0.75)
expressed_genes = gene_exp_data[which(gene_exp_data[,"mip_LOG2_I_U"] >= expression_cutoff | gene_exp_data[,"fur_LOG2_I_U"] >= expression_cutoff),]
length(expressed_genes[,1])

#10c) Now identify those without an absolute DE value of > 1 (2-fold)
de_cutoff = log2(2)
de_genes = gene_exp_data[which(abs(gene_exp_data[,"mip_v_fur_DE_U"]) > de_cutoff),]
length(de_genes[,1])

#10d) Calculate P-values for DE and MTP correct
gene_data_mtp_filtered = gene_data_mtp[which(gene_data_mtp[,"AlexaGene_ID"] %in% expressed_genes[,"AlexaGene_ID"]),]
gene_pvals = apply(gene_data_mtp_filtered, 1, compute_pvalues)
summary (gene_pvals)
procs = c("Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
res_gene = mt.rawp2adjp(gene_pvals, procs)
adjp_gene = res_gene$adjp[order(res_gene$index),"BH"]
gene_data_mtp_filtered[,"adjp_BH"] = adjp_gene
summary(adjp_gene)
length(which(adjp_gene < 0.05))

#10e) Now filter out those genes that had low absolute DE values
gene_data_mtp_de = gene_data_mtp_filtered[which(gene_data_mtp_filtered[,"AlexaGene_ID"] %in% de_genes[,"AlexaGene_ID"]),]

#The final number of DE genes that will be considered 'significant'
length(which(gene_data_mtp_de[,"adjp_BH"] < 0.05))





#11a.) Specifically examine the expression values for probesets corresponding to isoforms of UMPS!
#Isoform A probe data 
E2_data = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"]==2381692),6:41]
E1_E2_data = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"]==602470),6:41]
E2_E3_data = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"]==602475),6:41]

#Isoform B probe data
E1_E3_a_data = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"]==602471),6:41]
#E1_E3_b_data = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"]==602466),6:41]

Sensitive_isoformA = c(as.numeric(E1_E2_data[1:18]), as.numeric(E2_data[1:18]), as.numeric(E2_E3_data[1:18]))
Resistant_isoformA = c(as.numeric(E1_E2_data[19:36]), as.numeric(E2_data[19:36]), as.numeric(E2_E3_data[19:36]))
Sensitive_isoformB = c(as.numeric(E1_E3_a_data[1:18]))
Resistant_isoformB = c(as.numeric(E1_E3_a_data[19:36]))
#Sensitive_isoformB = c(as.numeric(E1_E3_a_data[1:18]), as.numeric(E1_E3_b_data[1:18]))
#Resistant_isoformB = c(as.numeric(E1_E3_a_data[19:36]), as.numeric(E1_E3_b_data[19:36]))


#Plot boxplots of these individual probe expression values
#Compare these to mean expression of all canonical probes and negative control probesets

canonical_probesets = exon_exp_data[which((exon_exp_data[,"Probe_Type"] == "Exon") | (exon_exp_data[,"Probe_Type"] == "Exon-Exon" & exon_exp_data[,"Exons_Skipped"] == 0)),]
canonical_exp_values = c(canonical_probesets[,"mip_LOG2_I_U"], canonical_probesets[,"fur_LOG2_I_U"])
negative_probesets = exon_exp_data[which(exon_exp_data[,"Probe_Type"] == "Control-Negative"),]
negative_exp_values = c(negative_probesets[,"mip_LOG2_I_U"], negative_probesets[,"fur_LOG2_I_U"])
 
#x = list(Sensitive_isoformA, Resistant_isoformA, Sensitive_isoformB, Resistant_isoformB, canonical_exp_values, negative_exp_values)
x = list(Sensitive_isoformA, Sensitive_isoformB,Resistant_isoformA, Resistant_isoformB, canonical_exp_values, negative_exp_values)


boxplot(x, main="Expression of UMPS Isoforms in Sensitive and Resistant Cells", ylab="Log2 Expression value from ALEXA platform", 
	  names=c("Isoform A","Isoform B","Isoform A","Isoform B","All exons","NC"), 
	  col=c("#92C2E8","#FFF100","#92C2E8","#FFF100","Purple","Red"),
	  border = c("Black","Black","Black","Black","Black","Black"))

wilcox.test(Sensitive_isoformA, Resistant_isoformA, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
wilcox.test(Sensitive_isoformB, Resistant_isoformB, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)

#11b.) Do the same thing except using the probe-level Affy data for the variant exon of UMPS (Affy probeset: 2639881)
# - Use the AO affy data (limited to only tose genese covered on alexa array and processed the same way)
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/"
dir(datadir)
setwd(datadir)
data_control_affy = read.table("AO_affy_probe_data.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4), na.strings='na')
data_control_affy = data_control_affy[which(data_control_affy[,"Probe_Type"]=="Control-Negative"),]
negative_exp_values = c(data_control_affy[,"mip_LOG2_I_U"], data_control_affy[,"fur_LOG2_I_U"])


#Load data for all exon and intron probesets that map to genes targeted by the ALEXA design
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/"
dir(datadir)
setwd(datadir)
exon_exp_data = read.table(file="AO_affy_probe_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])

canonical_probesets = exon_exp_data[which((exon_exp_data[,"Probe_Type"] == "Exon")),]
length(canonical_probesets[,1])
canonical_exp_values = c(canonical_probesets[,"mip_LOG2_I_U"], canonical_probesets[,"fur_LOG2_I_U"])

Sensitive_isoformA = c(7.8747,9.5199,9.9979,11.5256,7.1135,9.6998,10.1315,11.6842,7.5388,9.9931,9.9345,11.4948)
Resistant_isoformA = c(5.7746,7.3475,7.7847,9.2298,5.3266,6.9918,7.2592,9.5943,5.4473,6.8650,7.2707,9.2710)	
x = list(Sensitive_isoformA,Resistant_isoformA,canonical_exp_values,negative_exp_values)

boxplot(x, main="Expression of UMPS Isoforms in Sensitive and Resistant Cells", ylab="Log2 Expression value from Affymetrix platform", 
	  names=c("Isoform A","Isoform A","All exons","NC"), 
	  col=c("#92C2E8","#92C2E8","Purple","Red"),
	  border = c("Black","Black","Black","Black"))

wilcox.test(Sensitive_isoformA, Resistant_isoformA, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)



