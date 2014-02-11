#The purpose of this script is to compare the expression values (processed data - BGC, quantiles) for common probesets
#All of the data compared below corresponds to probesets that map to a common exon region in both platforms
#For both platforms: 
#- AA_probesets were defined such that all of the probes from both platforms that cover a particular exon may be compared directly

#The comparison of data will be done at the level of gene values
#A MTP correction will also be done and the number of genes emerging as significantly DE in one or both platforms will be noted

#1.) COMPARISON OF GENE EXPRESSION AND DIFFERENTIAL EXPRESSION VALUES
library(geneplotter)
library(RColorBrewer)
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AA_probesets"
dir(datadir)
setwd(datadir)

#Import the gene-level (probeset) data for common probesets for both platforms

#Alexa probeset data - These are values that have been summarized to the probeset level (combining all individual probes)
alexa_gene_data=read.table(file="AA_alexa_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

#Affy probeset data - These are values that have been summarized to the probeset level (combining all individual probes)
affy_gene_data=read.table(file="AA_affy_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

#Test the length of these datasets - Since these are common probesets they should be the same length!
length(alexa_gene_data[,"Probe_Count"])
length(affy_gene_data[,"Probe_Count"])

#Summarize the distribution of probeset sizes (number of probes making up each probeset) on both platforms
table(alexa_gene_data[,"Probe_Count"])
table(affy_gene_data[,"Probe_Count"])

#INTENSITY DATA CORRELATIONS
#Calculate correlations between replicates within each platform
alexa_mip_I_names = c("MIP_EF_A_LOG2_I","MIP_EF_B_LOG2_I","MIP_GH_A_LOG2_I")
alexa_fur_I_names = c("FUR_EF_A_LOG2_I","FUR_EF_B_LOG2_I","FUR_GH_A_LOG2_I")
alexa_mip_cor = cor(alexa_gene_data[,alexa_mip_I_names], method = "spearman")
alexa_fur_cor = cor(alexa_gene_data[,alexa_fur_I_names], method = "spearman")
alexa_de_names = c("ef_a_LOG2_DE","ef_b_LOG2_DE","gh_a_LOG2_DE")
alexa_de_cor = cor(alexa_gene_data[,alexa_de_names], method = "pearson")

affy_mip_I_names = c("MIP_C_LOG2_I","MIP_EF_B_LOG2_I","MIP_GH_A_LOG2_I")
affy_fur_I_names = c("FUR_C_LOG2_I","FUR_EF_B_LOG2_I","FUR_GH_A_LOG2_I")
affy_mip_cor = cor(affy_gene_data[,affy_mip_I_names], method = "spearman")
affy_fur_cor = cor(affy_gene_data[,affy_fur_I_names], method = "spearman")
affy_de_names = c("c_LOG2_DE","ef_b_LOG2_DE","gh_a_LOG2_DE")
affy_de_cor = cor(affy_gene_data[,affy_de_names], method = "pearson")

#ALEXA
alexa_mip_cor
alexa_fur_cor
alexa_de_cor

#AFFY
affy_mip_cor
affy_fur_cor
affy_de_cor

#Now calculate correlations between replicates across the two platforms (try replicate and mean values)
alexa_vs_affy_mip_cor = cor(x=alexa_gene_data[,alexa_mip_I_names], y=affy_gene_data[,affy_mip_I_names], method="spearman")
alexa_vs_affy_fur_cor = cor(x=alexa_gene_data[,alexa_fur_I_names], y=affy_gene_data[,affy_fur_I_names], method="spearman")
alexa_vs_affy_de_cor = cor(x=alexa_gene_data[,alexa_de_names], y=affy_gene_data[,affy_de_names], method="pearson")

alexa_vs_affy_mean_mip_cor = cor(x = alexa_gene_data[,"mip_LOG2_I_U"], y = affy_gene_data[,"mip_LOG2_I_U"], method = "spearman") 
alexa_vs_affy_mean_fur_cor = cor(x = alexa_gene_data[,"fur_LOG2_I_U"], y = affy_gene_data[,"fur_LOG2_I_U"], method = "spearman") 
alexa_vs_affy_mean_de_cor = cor(x = alexa_gene_data[,"mip_v_fur_DE_U"], y = affy_gene_data[,"mip_v_fur_DE_U"], method = "pearson") 

#Summary
alexa_vs_affy_mip_cor
alexa_vs_affy_fur_cor
alexa_vs_affy_de_cor

alexa_vs_affy_mean_mip_cor
alexa_vs_affy_mean_fur_cor
alexa_vs_affy_mean_de_cor

#Plot the Gene expression values from each platform against each other
#Actual log2 expression values for MIP samples
x_alexa_mip = alexa_gene_data[,"mip_LOG2_I_U"]
y_affy_mip = affy_gene_data[,"mip_LOG2_I_U"]
plot (x_alexa_mip, y_affy_mip, xlab="Alexa Log2 Gene Expression (means)", ylab="Affymetrix Log2 Gene Expression (means)", main="ALEXA and Affymetrix Expression Values (Sensitive data)", 
      col="blue", xlim=c(4,16), ylim=c(4,16))
test = cbind (x_alexa_mip,y_affy_mip)
lines(lowess(test), col="red", lwd=3, lty=2)
abline(0,1, col="black", lwd=2)

#Try smoothScatter plot
x_alexa_mip = alexa_gene_data[,"mip_LOG2_I_U"]
y_affy_mip = affy_gene_data[,"mip_LOG2_I_U"]
test = cbind (x_alexa_mip,y_affy_mip)
smoothScatter (x_alexa_mip, y_affy_mip, xlab="Alexa Log2 Gene Expression (means)", ylab="Affymetrix Log2 Gene Expression (means)", main="ALEXA and Affymetrix Expression Values (Sensitive data)", 
               xlim=c(4,16), ylim=c(4,16))
lines(lowess(test), col="red", lwd=3, lty=2)
abline(0,1, col="black", lwd=2)

#Plot ranks
plot (rank(x_alexa_mip), rank(y_affy_mip), xlab="Alexa Log2 Gene Expression (rank)", ylab="Affymetrix Log2 Gene Expression (rank)", main="Correlation of ALEXA and Affymetrix Expression Ranks (MIP101)", col="blue")
abline(0,1, col="red", lwd=2)

#Actual log2 expression values for 5FUR samples
x_alexa_fur = alexa_gene_data[,"fur_LOG2_I_U"]
y_affy_fur = affy_gene_data[,"fur_LOG2_I_U"]
plot (x_alexa_fur, y_affy_fur, xlab="Alexa Log2 Gene Expression (means)", ylab="Affymetrix Log2 Gene Expression (means)", main="ALEXA and Affymetrix Expression Values (Resistant data)", 
      col="blue", xlim=c(4,16), ylim=c(4,16))
test = cbind (x_alexa_fur,y_affy_fur)
lines(lowess(test), col="red", lwd=3, lty=2)
abline(0,1, col="black", lwd=2)

#Try smoothScatter plot
x_alexa_fur = alexa_gene_data[,"fur_LOG2_I_U"]
y_affy_fur = affy_gene_data[,"fur_LOG2_I_U"]
smoothScatter (x_alexa_fur, y_affy_fur, xlab="Alexa Log2 Gene Expression (means)", ylab="Affymetrix Log2 Gene Expression (means)", main="ALEXA and Affymetrix Expression Values (Resistant data)", 
               xlim=c(4,16), ylim=c(4,16))
test = cbind (x_alexa_fur,y_affy_fur)
lines(lowess(test), col="red", lwd=3, lty=2)
abline(0,1, col="black", lwd=2)


#Plot ranks
plot (rank(x_alexa_fur), rank(y_affy_fur), xlab="Alexa Log2 Gene Expression (rank)", ylab="Affymetrix Log2 Gene Expression (rank)", main="Correlation of ALEXA and Affymetrix Expression Ranks (MIP-5FUR)", col="blue")
abline(0,1, col="red", lwd=2)


#Plot the gene DE values from each platform against each other
x_alexa = alexa_gene_data[,"mip_v_fur_DE_U"]
y_affy = affy_gene_data[,"mip_v_fur_DE_U"]
plot (x_alexa, y_affy, xlab="Alexa Log2 Differential Expression", ylab="Affymetrix Log2 Differential Expression", main="Correlation of ALEXA and Affymetrix Gene DE Values", col="blue")
lm_fitted = predict(lm(y_affy~x_alexa))
lines(x_alexa, lm_fitted, col="red", lwd=2)

#Try smoothScatter plot
x_alexa = alexa_gene_data[,"mip_v_fur_DE_U"]
y_affy = affy_gene_data[,"mip_v_fur_DE_U"]
smoothScatter (x_alexa, y_affy, xlab="Alexa Log2 Differential Expression", ylab="Affymetrix Log2 Differential Expression", main="Correlation of ALEXA and Affymetrix Gene DE Expression Values",
 		   colramp=colorRampPalette(c("white", "blue", "dark blue")))
lm_fitted = predict(lm(y_affy~x_alexa))
lines(x_alexa, lm_fitted, col="red", lwd=2)

display.brewer.all()
colors()
col2rgb("yellow")

x_alexa = alexa_gene_data[,"mip_v_fur_DE_U"]
y_affy = affy_gene_data[,"mip_v_fur_DE_U"]
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#colors = colorRampPalette(c("white", "#007FFF", "blue", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#colors = colorRampPalette(c("white", brewer.pal(9, "Blues")))
#colors = colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))
#colors = colorRampPalette(c("white", brewer.pal(9, "RdYlGn")))

smoothScatter (x_alexa, y_affy, xlab="Alexa Log2 Differential Expression", ylab="Affymetrix Log2 Differential Expression", main="Correlation of ALEXA and Affymetrix Gene DE Expression Values",
 		   colramp=colors, nbin=275) # nbin=225 
lm_fitted = predict(lm(y_affy~x_alexa))
lines(x_alexa, lm_fitted, col="red", lwd=2)



#2.)IDENTIFY SIGNIFICANT DE EXONS ON BOTH PLATFORMS, CORRECT FOR MTP AND DESCRIBE OVERLAP OF THE TWO DATASETS
library(multtest)
library(nortest)

#Open MTP formatted files for both platforms
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AA_probesets/MTP_formatted"
dir(datadir)
setwd(datadir)

#Alexa gene data - These are values that have been summarized to the gene level (combining all individual probes)
alexa_gene_data=read.table(file="AA_alexa_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

#Affy gene data - These are values that have been summarized to the gene level (combining all individual probes)
affy_gene_data=read.table(file="AA_affy_gene_log2_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

#Check length of these two datasets
length(alexa_gene_data[,"Probe_Count"])
length(affy_gene_data[,"Probe_Count"])


compute_wilcox_values = function(dataset){

  #Get the number of columns
  column_count = length(dataset)

  #NOTE: First 2 columns contain extra info (not data), last column is empty
  variable_data = dataset[3:(column_count-1)]

  half_data = length(variable_data)/2
  result = wilcox.test(as.numeric(variable_data[1:half_data]), y = as.numeric(variable_data[(half_data+1):(length(variable_data))]) , alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
  return(result$p.value)
}

alexa_pvals = apply(alexa_gene_data, 1, compute_wilcox_values)
affy_pvals = apply(affy_gene_data, 1, compute_wilcox_values)

alexa_gene_data[,"raw_wilcox_pvals"] = alexa_pvals
affy_gene_data[,"raw_wilcox_pvals"] = affy_pvals

alexa_gene_data[,"affy_raw_wilcox_pvals"] = affy_pvals
affy_gene_data[,"alexa_raw_wilcox_pvals"] = alexa_pvals

alexa_sig = which((alexa_gene_data[,"raw_wilcox_pvals"] < 0.05) & (alexa_gene_data[,"affy_raw_wilcox_pvals"] >= 0.05))
affy_sig = which((alexa_gene_data[,"affy_raw_wilcox_pvals"] < 0.05) & (alexa_gene_data[,"raw_wilcox_pvals"] >= 0.05))
both_sig = which((alexa_gene_data[,"affy_raw_wilcox_pvals"] < 0.05) & (alexa_gene_data[,"raw_wilcox_pvals"] < 0.05))

#Summarize overlap of raw p-values
length(alexa_sig)
length(affy_sig)
length(both_sig)

#Now do the multiple testing correction
procs = c("Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
res = mt.rawp2adjp(alexa_pvals, procs)
adjp = res$adjp[order(res$index), ]
alexa_gene_data[,"adjp_BH"] = adjp[,"BH"]

res = mt.rawp2adjp(affy_pvals, procs)
adjp = res$adjp[order(res$index), ]
affy_gene_data[,"adjp_BH"] = adjp[,"BH"]

alexa_gene_data[,"affy_adjp_BH"] = affy_gene_data[,"adjp_BH"]
affy_gene_data[,"alexa_adjp_BH"] = alexa_gene_data[,"adjp_BH"]

alexa_sig = which((alexa_gene_data[,"adjp_BH"] < 0.05) & ((alexa_gene_data[,"affy_adjp_BH"] >= 0.05) | is.na(alexa_gene_data[,"affy_adjp_BH"])))
affy_sig = which((alexa_gene_data[,"affy_adjp_BH"] < 0.05) & ((alexa_gene_data[,"adjp_BH"] >= 0.05) | is.na(alexa_gene_data[,"adjp_BH"])))
both_sig = which((alexa_gene_data[,"affy_adjp_BH"] < 0.05) & (alexa_gene_data[,"adjp_BH"] < 0.05))
either_sig = which((alexa_gene_data[,"affy_adjp_BH"] < 0.05) | (alexa_gene_data[,"adjp_BH"] < 0.05))

#Summarize overlap of raw p-values
length(either_sig)
length(alexa_sig)
length(affy_sig)
length(both_sig)
significant_observations = length(both_sig)

#Define a class in each dataframe based on p-values 
alexa_gene_data[,"Pvalue_Class"] = "Neither"
alexa_gene_data[alexa_sig,"Pvalue_Class"] = "ALEXA"
alexa_gene_data[affy_sig,"Pvalue_Class"] = "AFFY"
alexa_gene_data[both_sig,"Pvalue_Class"] = "BOTH"

affy_gene_data[,"Pvalue_Class"] = "Neither"
affy_gene_data[alexa_sig,"Pvalue_Class"] = "ALEXA"
affy_gene_data[affy_sig,"Pvalue_Class"] = "AFFY"
affy_gene_data[both_sig,"Pvalue_Class"] = "BOTH"

#What is the distribution of observed intensities for each of these classes?
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AA_probesets"
dir(datadir)
setwd(datadir)

#Alexa probeset data - These are values that have been summarized to the probeset level (combining all individual probes)
alexa_summary_data=read.table(file="AA_alexa_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))
affy_summary_data=read.table(file="AA_affy_gene_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"))

#Using Alexa expression values
alexa_sig_genes = c(alexa_summary_data[alexa_sig,"mip_LOG2_I_U"],alexa_summary_data[alexa_sig,"fur_LOG2_I_U"]) 
both_sig_genes = c(alexa_summary_data[both_sig,"mip_LOG2_I_U"],alexa_summary_data[both_sig,"fur_LOG2_I_U"])
affy_sig_genes = c(alexa_summary_data[affy_sig,"mip_LOG2_I_U"],alexa_summary_data[affy_sig,"fur_LOG2_I_U"])
x = list(alexa_sig_genes,both_sig_genes,affy_sig_genes)
boxplot(x, main="Distribution of expression for each class of significant DE genes", ylab="Gene expression values from ALEXA platform", names=c("ALEXA only","BOTH","AFFY only"), col=c("Light Blue","Blue","Dark Green"))

#Normalcy test (n, qq-plot, p-value for whether we can prove that the distribution in NOT normal)
length(alexa_sig_genes)
qqnorm(alexa_sig_genes)
qqline(alexa_sig_genes, col="red", lwd=2)
ad.test(alexa_sig_genes)
length(affy_sig_genes)
qqnorm(affy_sig_genes)
qqline(affy_sig_genes, col="red", lwd=2)
ad.test(affy_sig_genes)
#t.test(x = alexa_sig_genes, y = affy_sig_genes, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
wilcox.test(alexa_sig_genes, affy_sig_genes, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)



#Using Affy expression values
alexa_sig_genes = c(affy_summary_data[alexa_sig,"mip_LOG2_I_U"],affy_summary_data[alexa_sig,"fur_LOG2_I_U"]) 
both_sig_genes = c(affy_summary_data[both_sig,"mip_LOG2_I_U"],affy_summary_data[both_sig,"fur_LOG2_I_U"])
affy_sig_genes = c(affy_summary_data[affy_sig,"mip_LOG2_I_U"],affy_summary_data[affy_sig,"fur_LOG2_I_U"])
x = list(alexa_sig_genes,both_sig_genes,affy_sig_genes)
boxplot(x, main="Distribution of expression for each class of significant DE genes", ylab="Gene expression values from Affymetrix platform", names=c("ALEXA only","BOTH","AFFY only"), col=c("Light Blue","Blue","Dark Green"))

#Normalcy test (n, qq-plot, p-value for whether we can prove that the distribution in NOT normal)
length(alexa_sig_genes)
qqnorm(alexa_sig_genes)
qqline(alexa_sig_genes, col="red", lwd=2)
ad.test(alexa_sig_genes)
length(affy_sig_genes)
qqnorm(affy_sig_genes)
qqline(affy_sig_genes, col="red", lwd=2)
ad.test(affy_sig_genes)

#t.test(x = alexa_sig_genes, y = affy_sig_genes, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
wilcox.test(alexa_sig_genes, affy_sig_genes, alternative = "two.sided", mu = 0, paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)


#Write out appended data files containing the desired P-values
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AA_probesets/MTP_formatted"
dir(datadir)
setwd(datadir)

alexa_data_out = alexa_gene_data[,c("AlexaGene_ID","Probe_Count","raw_wilcox_pvals","adjp_BH","affy_raw_wilcox_pvals","affy_adjp_BH","Pvalue_Class")]
alexa_data_out[,"raw_wilcox_pvals"] = format(alexa_data_out[,"raw_wilcox_pvals"], scientific=FALSE)
alexa_data_out[,"adjp_BH"] = format(alexa_data_out[,"adjp_BH"], scientific=FALSE)
alexa_data_out[,"affy_raw_wilcox_pvals"] = format(alexa_data_out[,"affy_raw_wilcox_pvals"], scientific=FALSE)
alexa_data_out[,"affy_adjp_BH"] = format(alexa_data_out[,"affy_adjp_BH"], scientific=FALSE)

write.table(alexa_data_out, file = "AA_ALEXA_AFFY_gene_pvals.txt", append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "na", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = c("escape", "double"))


#Determine the significance of the observed overlap (511) between methods using a randomization test
gene_count = length(alexa_gene_data[,1])
indexes = 1:gene_count

#Get a random sample of the p-values observed for each dataset
rand_alexa_index = sample(indexes, gene_count, replace = FALSE)
rand_affy_index = sample(indexes, gene_count, replace = FALSE)
length(rand_alexa_index)

#Create new dataframes using these randomized indexes
rand_alexa_pvals = alexa_gene_data[rand_alexa_index,"adjp_BH"]
rand_affy_pvals = affy_gene_data[rand_affy_index,"adjp_BH"]

alexa_sig = which((rand_alexa_pvals < 0.05) & ((rand_affy_pvals >= 0.05) | is.na(rand_affy_pvals)))
affy_sig = which((rand_affy_pvals < 0.05) & ((rand_alexa_pvals >= 0.05) | is.na(rand_alexa_pvals)))
both_sig = which((rand_alexa_pvals < 0.05) & (rand_affy_pvals < 0.05))
either_sig = which((rand_alexa_pvals < 0.05) | (rand_affy_pvals < 0.05))

#Note: (alexa_sig + both_sig) is the same for the original dataset and every random dataset  
#Similary (affy_sig + both_sig) is the same for the original dataset and every random dataset  
length(either_sig)
length(alexa_sig)
length(affy_sig)
length(both_sig)

#Now try this many times and determine how ofter 'both_sig' is >= 511 by chance!
iter = 1000
random_test = matrix(data = NA, nrow = iter, ncol = 4, dimnames = list(1:iter, c("total", "alexa", "affy", "both")))
exon_count = length(alexa_gene_data[,1])
indexes = 1:gene_count
for (i in 1:iter){
  #Get a random sample of the p-values observed for each dataset
  rand_alexa_index = sample(indexes, exon_count, replace = FALSE)
  rand_affy_index = sample(indexes, exon_count, replace = FALSE)

  #Create new dataframes using these randomized indexes
  rand_alexa_pvals = alexa_gene_data[rand_alexa_index,"adjp_BH"]
  rand_affy_pvals = affy_gene_data[rand_affy_index,"adjp_BH"]

  alexa_sig = which((rand_alexa_pvals < 0.05) & ((rand_affy_pvals >= 0.05) | is.na(rand_affy_pvals)))
  affy_sig = which((rand_affy_pvals < 0.05) & ((rand_alexa_pvals >= 0.05) | is.na(rand_alexa_pvals)))
  both_sig = which((rand_alexa_pvals < 0.05) & (rand_affy_pvals < 0.05))

  random_test[i, "total"] = exon_count
  random_test[i, "alexa"] = length(alexa_sig) + length(both_sig)
  random_test[i, "affy"] = length(affy_sig) + length(both_sig)
  random_test[i, "both"] = length(both_sig)

  #Calculate apparent p-value thus far
  cum_pvalue = (length(which(random_test[,"both"] >= significant_observations))) / i 

  message1 = paste(i, "Apparent Pvalue = ", cum_pvalue)
  print (message1)
}

#Check out the results to see if they make sense
summary(random_test[, "both"])

#Write out the results and save in case you want to add multiple runs together
#If the file is aleady there, append your results

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/AA_probesets/MTP_formatted"
dir(datadir)
setwd(datadir)

write.table(as.data.frame(random_test), file = "pvalue_permutations_genes.txt", append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "na", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = c("escape", "double"))

























