#Written by Malachi Griffith
#Identify a list of candidate significant events (exons, junctions, boundaries) that ARE NOT the result of DE at the gene level

#This will involve:
#1.) first filtering events to eliminate those that have poor evidence for expression in either condition
#2.) then filtering to eliminate those with only a minor change in expression between conditions
#3.) also filter those where the gene had a large DE value
#3.) then calculating a p-value for this observed change between conditions
#4.) then correcting this p-value for multiple testing
#5.) then filtering out those events with a p-value > 0.05
#6.) and finally printing out the final list of candidate events

#Load neccessary libraries
library(multtest)

#Load exon expression, DE, and SI values for each probeset and gene

#GENE values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)
gene_exp_data = read.table(file="Standard_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
names(gene_exp_data)[14] = "gene_mip_v_fur_DE_U"
length(gene_exp_data[,1])

#Probeset I and DE values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)
exon_exp_data = read.table(file="Standard_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])

#Probeset GNI and SI values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values"
dir(datadir)
setwd(datadir)
exon_si_data = read.table(file="Standard_probeset_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))

#1.) first filtering events to eliminate those that have poor evidence for expression in either condition
#What is a good expression cutoff?  Base this on the mean expression values of negative control probes.
nc_exon_exp_data = exon_exp_data[which(exon_exp_data[,"Probe_Type"] == "Control-Negative"),]
length(nc_exon_exp_data[,1])
quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=seq(0, 1, 0.05))
quantile(nc_exon_exp_data[,"fur_LOG2_I_U"], probs=seq(0, 1, 0.05))
expression_cutoff = quantile(nc_exon_exp_data[,"mip_LOG2_I_U"], probs=0.975) #0.975

#Now identify the probesets that meet the cutoff
expressed_probesets = exon_exp_data[which(exon_exp_data[,"mip_LOG2_I_U"] >= expression_cutoff | exon_exp_data[,"fur_LOG2_I_U"] >= expression_cutoff),]
length(expressed_probesets[,1])

#2) Determine those with an absolute SI value of > 1 (2-fold)
si_cutoff = log2(2)
si_probesets = exon_si_data[which(abs(exon_si_data[,"mip_v_fur_SI_U"]) > si_cutoff),]
length(si_probesets[,1])

#3.) Determine those with an absolute DE value of > 1 (2-fold)
#    - Only do this for alternative events!  For exons and canonical junctions we want to be able to identify when a few exons do not 
#      change and the rest of the gene does.
de_cutoff = log2(2)
de_probesets = exon_exp_data[which((abs(exon_exp_data[,"mip_v_fur_DE_U"]) >= de_cutoff) | (exon_exp_data[,"Probe_Type"] == "Exon") | (exon_exp_data[,"Probe_Type"] == "Exon-Exon" & exon_exp_data[,"Exons_Skipped"] == 0)),]
length(de_probesets[,1])

#4.) Determine those with a abs(SI_event - DE_gene) > log2(3)
event_si_vs_gene_de_cutoff = log2(3)
exon_si_data = merge(x=exon_si_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_si_data[,"event_SI_vs_gene_DE"] = exon_si_data[,"mip_v_fur_SI_U"] - exon_si_data[,"gene_mip_v_fur_DE_U"]
event_si_vs_gene_de_probesets = exon_si_data[which(abs(exon_si_data[,"event_SI_vs_gene_DE"]) > event_si_vs_gene_de_cutoff),]
length(event_si_vs_gene_de_probesets[,1])

#5.) Determine those with an absolute GENE DE value less than (8-fold)
gene_de_cutoff = log2(8)
non_de_gene_probesets = exon_si_data[which(abs(exon_si_data[,"gene_mip_v_fur_DE_U"]) < gene_de_cutoff),]
length(non_de_gene_probesets[,1])

#Actually load the MTP formatted data
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/H19_excluded/"
dir(datadir)
setwd(datadir)
exon_data_mtp = read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(exon_data_mtp[,1])

#Load H19 data seperately
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/H19_only/"
dir(datadir)
setwd(datadir)
exon_data_mtp_H19 = read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
length(exon_data_mtp_H19[,1])

#6.) Now filter out probesets according to the filters defined above

#Using the lists of probesets determined above, grab the subset of probesets which passed the filters
exon_data_mtp_filtered_1 = exon_data_mtp[which(exon_data_mtp[,"ProbeSet_ID"] %in% expressed_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_1[,1])
exon_data_mtp_filtered_2 = exon_data_mtp_filtered_1[which(exon_data_mtp_filtered_1[,"ProbeSet_ID"] %in% si_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_2[,1])
exon_data_mtp_filtered_3 = exon_data_mtp_filtered_2[which(exon_data_mtp_filtered_2[,"ProbeSet_ID"] %in% de_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_3[,1])
exon_data_mtp_filtered_4 = exon_data_mtp_filtered_3[which(exon_data_mtp_filtered_3[,"ProbeSet_ID"] %in% event_si_vs_gene_de_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_4[,1])
exon_data_mtp_filtered_5 = exon_data_mtp_filtered_4[which(exon_data_mtp_filtered_4[,"ProbeSet_ID"] %in% non_de_gene_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_filtered_5[,1])
exon_data_mtp_filtered = exon_data_mtp_filtered_5

exon_data_mtp_H19_filtered_1 = exon_data_mtp_H19[which(exon_data_mtp_H19[,"ProbeSet_ID"] %in% expressed_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_1[,1])
exon_data_mtp_H19_filtered_2 = exon_data_mtp_H19_filtered_1[which(exon_data_mtp_H19_filtered_1[,"ProbeSet_ID"] %in% si_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_2[,1])
exon_data_mtp_H19_filtered_3 = exon_data_mtp_H19_filtered_2[which(exon_data_mtp_H19_filtered_2[,"ProbeSet_ID"] %in% de_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_3[,1])
exon_data_mtp_H19_filtered_4 = exon_data_mtp_H19_filtered_3[which(exon_data_mtp_H19_filtered_3[,"ProbeSet_ID"] %in% event_si_vs_gene_de_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_4[,1])
exon_data_mtp_H19_filtered_5 = exon_data_mtp_H19_filtered_4[which(exon_data_mtp_H19_filtered_4[,"ProbeSet_ID"] %in% non_de_gene_probesets[,"ProbeSet_ID"]),]
length(exon_data_mtp_H19_filtered_5[,1])
exon_data_mtp_H19_filtered = exon_data_mtp_H19_filtered_5


#7.) Now calculate a p-value for the observed change between conditions
# For this we need to load the MTP formatted data file

#Start with exon level log2 values for 3 replicates and 2 conditions (sensitive versus resistant)
#Each row contains all of the probe intensity observations for a probeset (which corresponds to an exon) for all replicates 
#Make sure the file is tab-delimited and missing values are set to NA 

#Try a simple Mann-Whitney test (aka Wilcoxen Rank Test)="t.twosamp.equalvar".  

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

#P-value correction for multiple testing
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

nc_data[,"adjp_BH"] = nc_pvals  #**Only one p-value so MTP Correction meaningless 
intron_data[,"adjp_BH"] = intron_pvals  #**Only one p-value so MTP Correction meaningless 
exon_boundary_data[,"adjp_BH"] = adjp_exon_boundary
exon_data[,"adjp_BH"] = adjp_exon
canonical_data[,"adjp_BH"] = adjp_canonical
skipping_data[,"adjp_BH"] = adjp_skipping
s1_data[,"adjp_BH"] = adjp_s1
s2_data[,"adjp_BH"] = adjp_s2
s3_data[,"adjp_BH"] = adjp_s3
ei_data[,"adjp_BH"] = adjp_ei
ie_data[,"adjp_BH"] = adjp_ie
exon_data_mtp_H19_filtered[,"adjp_BH"] = H19_pvals  #**Only one p-value so MTP Correction meaningless

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

#Gather actual mean SI values to be printed out with p-values
nc_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% nc_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
intron_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% intron_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
exon_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% exon_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
canonical_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% canonical_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
s1_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% s1_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
s2_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% s2_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
s3_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% s3_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
ei_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% ei_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
ie_data[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% ie_data[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]
exon_data_mtp_H19_filtered[,"mip_v_fur_SI_U"] = exon_si_data[which(exon_si_data[,"ProbeSet_ID"] %in% exon_data_mtp_H19_filtered[,"ProbeSet_ID"]),c("mip_v_fur_SI_U")]

#Add gene-level data
nc_data = merge(x=nc_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
intron_data = merge(x=intron_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_data = merge(x=exon_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
canonical_data = merge(x=canonical_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
s1_data = merge(x=s1_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
s2_data = merge(x=s2_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
s3_data = merge(x=s3_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
ei_data = merge(x=ei_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
ie_data = merge(x=ie_data, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)
exon_data_mtp_H19_filtered = merge(x=exon_data_mtp_H19_filtered, y=gene_exp_data[,c("AlexaGene_ID","gene_mip_v_fur_DE_U")], by.x=2, by.y=1)


#4.) Print out the probeset IDs for these events
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted"
dir(datadir)
setwd(datadir)
write.table(nc_data[which(nc_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(intron_data[which(intron_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
write.table(exon_data[which(exon_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(canonical_data[which(canonical_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(s1_data[which(s1_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(s2_data[which(s2_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(s3_data[which(s3_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(ei_data[which(ei_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(ie_data[which(ie_data[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
write.table(exon_data_mtp_H19_filtered[which(exon_data_mtp_H19_filtered[,"adjp_BH"] < 0.05),c("ProbeSet_ID","AlexaGene_ID","Probe_Type","Exons_Skipped","mip_v_fur_SI_U","adjp_BH","gene_mip_v_fur_DE_U")], file = "SigEvents_Log2Exp_FilteredOnExp_SI.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))

