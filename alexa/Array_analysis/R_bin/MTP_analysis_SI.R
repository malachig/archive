#Written by Malachi Griffith
#Process data with the MTP package to identify significant DE exons and genes

#Start with exon level log2 values for 3 replicates and 2 conditions (sensitive versus resistant)
#Each row contains all of the probe intensity observations for a probeset (which corresponds to an exon) for all replicates 
#Make sure the file is tab-delimited and missing values are set to NA 

library(multtest)
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/H19_excluded"
dir(datadir)
setwd(datadir)
rawdata=read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4), comment.char = "")

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/H19_only"
dir(datadir)
setwd(datadir)
rawdata_H19=read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4), comment.char = "")


#Make note of the total number of statistical tests that are going to be conducted here!!
#This will either be the number of probesets (exons, introns or junctions) or the number of gene
length(rawdata[,"Probe_Type"])
length(rawdata_H19[,"Probe_Type"])

#A.) Create a function that will get a p-value for a single row of values consisting of two populations
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
nc_data = rawdata[which(rawdata[,"Probe_Type"]=="Control-Negative"),]
intron_data = rawdata[which(rawdata[,"Probe_Type"]=="Intron"),]
exon_boundary_data = rawdata[which(rawdata[,"Probe_Type"]=="Intron-Exon" | rawdata[,"Probe_Type"]=="Exon-Intron"),]
exon_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon"),]
canonical_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 0),]
skipping_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] >= 1),]

#further breakdowns
s1_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 1),]
s2_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 2),]
s3_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 3),]
ei_data = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Intron"),]
ie_data = rawdata[which(rawdata[,"Probe_Type"]=="Intron-Exon"),]

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

#Get p-values for the whole dataset and write the raw p-values to a new file
all_pvals = apply(rawdata, 1, compute_pvalues)
h19_pvals = apply(rawdata_H19, 1, compute_pvalues)

rawdata[,"raw_ttest_pval"] = all_pvals
rawdata_H19[,"raw_ttest_pval"] = h19_pvals

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/"
dir(datadir)
setwd(datadir)

write.table(rawdata[,c("ProbeSet_ID","AlexaGene_ID","Probe_Count","Probe_Type","Exons_Skipped","raw_ttest_pval")], 
		file = "exon_log2_data_pvals.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "na", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))

write.table(rawdata_H19[,c("ProbeSet_ID","AlexaGene_ID","Probe_Count","Probe_Type","Exons_Skipped","raw_ttest_pval")], 
		file = "exon_log2_data_pvals.txt", append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "na", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))


