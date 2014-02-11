#Written by Malachi Griffith
#Consider the Affymetrix exon array data for only those probes that occur within an EnsEMBL gene targeted on the ALEXA design
#All Affymetrix probesets are left as is on the array, all those within an EnsEMBL exon are called 'Exon' probes and all others are 'Intron'
#Identify a list of candidate significant events (exons or introns)
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

#1.) Import Affymetrix 'AO' and Negative control DATA
#Load exon expression values for each probeset

#Load Affymetrix Negative control data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel"
dir(datadir)
setwd(datadir)
data_hk_control_affy = read.table("HuEX-1_001-006_RawProbe_HK+NC_ControlsOnly.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,12,15), na.strings='na')
data_hk_control_affy = data_hk_control_affy[,c("Probe_Type","probe_id","probeset_id","MIP_C","FUR_C","MIP_GH_A","FUR_GH_A","MIP_EF_B","FUR_EF_B","Probe_Tm")]
data_control_affy = data_hk_control_affy[which(data_hk_control_affy[,"Probe_Type"]=="Control-Negative"),]

#Load data for all exon and intron probesets that map to genes targeted by the ALEXA design
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/"
dir(datadir)
setwd(datadir)
exon_exp_data = read.table(file="AO_affy_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])
table(exon_data_mtp[,"Probe_Type"])

#Load data for all exon and intron probes in an MTP suitable format
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/MTP_formatted"
dir(datadir)
setwd(datadir)
exon_data_mtp = read.table(file="AO_affy_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(exon_data_mtp[,1])


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

#5.) Now filter out probesets that appear to be invariant or 'non-responsive' across conditions

#Using the lists of probesets determined above, grab the subset of probesets which passed the filters
exon_data_mtp_filtered_1 = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"] %in% expressed_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_1[,1])
exon_data_mtp_filtered = exon_data_mtp_filtered_1[which(exon_data_mtp_filtered_1[,"ProbeSet_ID"] %in% de_probesets_list[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered[,1])


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
intron_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Intron"),]
exon_data = exon_data_mtp_filtered[which(exon_data_mtp_filtered[,"Probe_Type"]=="Exon"),]

intron_pvals = apply(intron_data, 1, compute_pvalues)
exon_pvals = apply(exon_data, 1, compute_pvalues)

#summarise the distribution of P-values for each of the groups - BEFORE MTP correction
summary (intron_pvals)
summary (exon_pvals)


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
res_intron = mt.rawp2adjp(intron_pvals, procs)
res_exon = mt.rawp2adjp(exon_pvals, procs)

#Get the adjusted p-values and store in the original order
adjp_intron = res_intron$adjp[order(res_intron$index),"BH"]
adjp_exon = res_exon$adjp[order(res_exon$index),"BH"]

intron_data[,"adjp_BH"] = adjp_intron
exon_data[,"adjp_BH"] = adjp_exon

#Summarize the adjusted p-values for each of the probe categories or sub-categories
summary(adjp_intron)
summary(adjp_exon)

length(which(adjp_intron < 0.05))
length(which(adjp_exon < 0.05))


#8.) Get the gene-level expression changes for these events
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/"
dir(datadir)
setwd(datadir)
gene_exp_data = read.table(file="AO_affy_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
names(gene_exp_data)[14] = "gene_mip_v_fur_DE_U"

#Gather actual mean DE values to be printed out with p-values
intron_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% intron_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]
exon_data[,"mip_v_fur_DE_U"] = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% exon_data[,"ProbeSet_ID"]),c("mip_v_fur_DE_U")]

#Add gene-level data
intron_data = merge(x=intron_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_data = merge(x=exon_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)

#9.) Print out the probeset IDs for these events
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/MTP_formatted"
dir(datadir)
setwd(datadir)
write.table(intron_data[which(intron_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(exon_data[which(exon_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_DE_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_DE.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))

#10.) Similar to the analysis for events above, identify GENES that are significantly DE between Sensitive and Resistant  

#10a) Load gene expression values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets"
dir(datadir)
setwd(datadir)
gene_exp_data = read.table(file="AO_affy_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(gene_exp_data[,1])

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets/MTP_formatted/"
dir(datadir)
setwd(datadir)
gene_data_mtp = read.table(file="AO_affy_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(gene_data_mtp[,1])

#10b) Use NC probeset expression values to apply a simple expression cutoff 
#Identify genes that are not expressed in either condition
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AO_probesets"
dir(datadir)
setwd(datadir)
exon_exp_data = read.table(file="AO_affy_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
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







