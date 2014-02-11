#Written by Malachi Griffith
#Process data with the MTP package to identify significant DE exons and genes

#Start with exon level log2 values for 3 replicates and 2 conditions (sensitive versus resistant)
#Each row contains all of the probe intensity observations for a probeset (which corresponds to an exon) for all replicates 
#Make sure the file is tab-delimited and missing values are set to NA 

library(multtest)
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_excluded"
dir(datadir)
setwd(datadir)
rawdata=read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4), comment.char = "")

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_only"
dir(datadir)
setwd(datadir)
rawdata_H19=read.table(file="Standard_exon_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4), comment.char = "")


#Make note of the total number of statistical tests that are going to be conducted here!!
#This will either be the number of probesets (exons, introns or junctions) or the number of gene
length(rawdata[,"Probe_Type"])
length(rawdata_H19[,"Probe_Type"])

#Try a simple t.test

#A.) Create a function that will get a p-value for a single row of values consisting of two populations
compute_pvalues = function(dataset){
  start_data_column = 6
  end_data_column = length(dataset)-1
  variable_data=dataset[start_data_column:end_data_column]
  half_data = length(variable_data)/2
  result = t.test(x=as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]),
                  alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
  return(result$p.value)
}

#compute_wilcox_values = function(dataset){
#  variable_data=dataset[6:41]
#  half_data = length(variable_data)/2
#  result = wilcox.test(as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]) , alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
#  return(result$p.value)
#}

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

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/"
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

#Calculate adjusted p-values using the 7 alternate procedures listed above
procs = c("Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
res = mt.rawp2adjp(all_pvals, procs)

#Get the adjusted p-values and store in the original order
adjp = res$adjp[order(res$index), ]
rawdata[,"adjp_BH"] = adjp[,"BH"]
summary(rawdata[,"adjp_BH"])

#Display distributions of P-values
hist(rawdata[,"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - All Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Control-Negative"),"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - NC Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Intron"),"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - Intron Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Intron-Exon" | rawdata[,"Probe_Type"]=="Exon-Intron"),"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - Exon Boundary Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Exon"),"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - Exon Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 0),"adjp_BH"], breaks=20, col="blue", main="Distribution of corrected p-values - Canonical junction Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")
hist(rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] >= 1),"adjp_BH"], breaks=35, col="blue", main="Distribution of corrected p-values - Skipping junction Probesets", ylab="Probeset Count", xlab="Adjusted P-value (BH method)")

#How many probes are actually significant at 0.05 and what percent is this?
x = rawdata[,"adjp_BH"]
summary(rawdata[,"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Control-Negative"),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Control-Negative"),"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Intron"),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Intron"),"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Intron-Exon" | rawdata[,"Probe_Type"]=="Exon-Intron"),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Intron-Exon" | rawdata[,"Probe_Type"]=="Exon-Intron"),"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Exon"),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Exon"),"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 0),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] == 0),"adjp_BH"])
length(which(x < 0.1))

x = rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] >= 1),"adjp_BH"]
summary(rawdata[which(rawdata[,"Probe_Type"]=="Exon-Exon" & rawdata[,"Exons_Skipped"] >= 1),"adjp_BH"])
length(which(x < 0.1))


##THE ONLY PROBES THAT COME OUT AS SIGNIFICANT AFTER MTP CORRECTION ARE EXON PROBES
#THIS IS PERHAPS BECAUSE THE JUNCTION PROBES TEND TO HAVE FEWER OBSERVATIONS
#TRY CORRECTING FOR MTP FOR EACH PROBESET SEPERATELY TO SEE IF THIS HELPS THE SITUATION
#EACH OF THE PROBESETS HAS A DIFFERENT EXPECTATION AND POTENTIALLY VARIANCE AS WELL
#THE ARGUMENT COULD BE MADE THAT THEY CAN BE FAIRLY PROCESSED SEPERATELY
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

length(which(adjp_exon < 0.05))
length(which(adjp_canonical < 0.05))

x = exon_data[which(adjp_exon < 0.05),"AlexaGene_ID"]
length(unique(x))

x = canonical_data[which(adjp_canonical < 0.05), "AlexaGene_ID"]
length(unique(x))



#B.) Try the same thing except this time use a paired test (sensitive versus resistant but WITHIN sample pairs)
#Create a function that will get a p-value for a single row of values consisting of two populations
#compute_wilcox_values_paired = function(dataset){
#  variable_data=dataset[6:41]
#  half_data = length(variable_data)/2
#  result = wilcox.test(as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]) , alternative = "two.sided", mu = 0, paired = TRUE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
#  return(result$p.value)
#}

#nc_pvals_paired = apply(nc_data, 1, compute_wilcox_values_paired)
#intron_pvals_paired = apply(intron_data, 1, compute_wilcox_values_paired)
#exon_boundary_pvals_paired = apply(exon_boundary_data, 1, compute_wilcox_values_paired)
#exon_pvals_paired = apply(exon_data, 1, compute_wilcox_values_paired)
#canonical_pvals_paired = apply(canonical_data, 1, compute_wilcox_values_paired)
#skipping_pvals_paired = apply(skipping_data, 1, compute_wilcox_values_paired)

#summarise the distribution of P-values for each of the groups
#summary (nc_pvals_paired)
#summary (intron_pvals_paired)
#summary (exon_boundary_pvals_paired)
#summary (exon_pvals_paired)
#summary (canonical_pvals_paired)
#summary (skipping_pvals_paired)

#This method doesnt seem to work well.  Gives lots of significant p-values for nc and intron probes, most likely due to 
#differences in background between sample pairs which were not adequately corrected ...


#C.) GENE-LEVEL ANALYSIS

#Each row contains all of the probe intensity observations for a  gene (exon + canonical probes) for all replicates 
#Make sure the file is tab-delimited and missing values are set to NA 

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_excluded"
dir(datadir)
setwd(datadir)
rawdata=read.table(file="Standard_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/H19_only"
dir(datadir)
setwd(datadir)
rawdata_H19=read.table(file="Standard_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))


#Make note of the total number of statistical tests that are going to be conducted here!!
#This will either be the number of probesets (exons, introns or junctions) or the number of gene
#Note that H19 is processed seperately in this analysis
length(rawdata[,"AlexaGene_ID"])
length(rawdata_H19[,"AlexaGene_ID"])

#Note that there are many columns in this data.  The first 2 are not data points!
compute_pvalues = function(dataset){
  if (dataset[2] < 1){
    return(NA)
  }

  start_data_column = 3
  end_data_column = length(dataset)-1
  variable_data=dataset[start_data_column:end_data_column]
  half_data = length(variable_data)/2
  result = t.test(x=as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]),
                  alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
  return(result$p.value)
}

#compute_wilcox_values = function(dataset){
#  variable_data=dataset[3:2108]
#  #Make sure the gene has at least some suitable probes!!
#  if (dataset[2] < 1){
#    return(NA)
#  }
#  half_data = length(variable_data)/2
#  result = wilcox.test(as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]) , alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
#  return(result$p.value)
#}

pvals = apply(rawdata[,], 1, compute_pvalues)
pvals_H19 = apply(rawdata_H19[,], 1, compute_pvalues)
pvals_combined = c(pvals,pvals_H19)

procs = c("Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
res = mt.rawp2adjp(pvals_combined, procs)
adjp = res$adjp[order(res$index),"BH"]

rawdata[,"raw_ttest_pval"] = pvals
rawdata_H19[,"raw_ttest_pval"] = pvals_H19

rawdata[,"adjp_BH"] = adjp[1:2508]
rawdata_H19[,"adjp_BH"] = adjp[2509]

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted/"
dir(datadir)
setwd(datadir)

write.table(rawdata[,c("AlexaGene_ID","Probe_Count","raw_ttest_pval","adjp_BH")], 
		file = "gene_log2_data_pvals.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "na", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))

write.table(rawdata_H19[,c("AlexaGene_ID","Probe_Count","raw_ttest_pval","adjp_BH")], 
		file = "gene_log2_data_pvals.txt", append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "na", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))











#MTP ITERATION METHOD

#Try processing different probe types individually
table(rawdata[,"Probe_Type"])

#NOTE: The following test fails for some probesets if the number of observations is too small relative to the number of NA's
#This might be due to the iterations resulting in populations repeatedly consisting entirely of NA's by chance 
#It is simpler to just use the one-step MTP method provided with multtest

conduct_mtp_test = function(dataset){
  #Create a matrix with only the variable data (i.e. exclude the first five columns of info)
  #Also watch out for a dead column at the end which contains no data!

  #The MTP function expects one test for each row of data (as is conventional in gene x exp matrices).
  variable_data=dataset[,6:41]

  #Create a vector for the class labels (sensitive vs resistant).  This determines what values are tested against what
  half_data = length(variable_data[1,])/2
  class_labels = array(data = "test", dim = length(variable_data[1,]), dimnames = NULL)
  class_labels[1:half_data] = "sensitive"
  class_labels[(half_data+1):(length(variable_data[1,]))] = "resistant"

  print (length(class_labels))
  print (length(variable_data[1,]))

  #Create a vector of variable names for output later
  col_names=colnames(variable_data)
  full_col_names=colnames(dataset)

  row_names = dataset[,1]

  #Run the MTP procedure.
  #The Mann-Whitney test (or Wilcoxen Rank Test)="t.twosamp.equalvar".  
  #You might also consider using the Welch T-test (t.twosamp.unequalvar).  Are variances equal?
  #Also, note that B=1000 (the number of permutations) or higher is recommended.  But, for testing 100 is convenient.
  #If the sample size is small, the authors suggest setting standardize=FALSE
  #If you are getting rawp values equal to 0 try setting smooth.null=TRUE
  #For robust version of tests you can set robust=TRUE.  I'm not sure what this really means.

  MTP_results=MTP(X=variable_data, Y=class_labels, na.rm=TRUE, test="t.twosamp.equalvar", alternative="two.sided", robust=FALSE, typeone="fdr", fdr.method="conservative", alpha=0.05, B=100, method="sd.minP", smooth.null=FALSE, standardize=TRUE)

  #Output results.
  MTP_summary=cbind(row_names,MTP_results@statistic,MTP_results@estimate,MTP_results@rawp, MTP_results@adjp, MTP_results@reject)
  colnames(MTP_summary)=c("variable", "statistic", "estimate", "rawp", "adjp", "reject")

  return(MTP_summary)
}
test_data = rawdata[which(rawdata[,"Probe_Type"] == "Intron"),]

conduct_mtp_test(test_data[1:100,])

