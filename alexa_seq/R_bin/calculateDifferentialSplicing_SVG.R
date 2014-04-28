#!/usr/bin/env Rscript
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Load the multtest library
#Note that this requires an installation of Bioconductor to your existing R installation
#http://www.bioconductor.org/docs/install/
library(multtest)
library(Cairo)

#Options:
#[1] dataname
#[2] datafile
#[3] outdir
#[4] shortname

#NOTES:
#The input file must contain: a unique ID (gene, exon, junction), Gene ID, three additional descriptive columns and four columns of expression data (Gene and Seq expression data for library A and B)
#These expression values should be raw expression values (normalization to the smaller library will be done by this script)
#The output file will contain all input columns, normalized data values, fold change, the splicing Index value, pvalue, and corrected pvalues (qvalue)
#Not that a basic expression filter was applied when creating the input files for this script.  The Gene AND Feature must be expressed in at least one of the two conditions
args = (commandArgs(TRUE))
dataname = args[1];
datafile = args[2];
outdir = args[3];
shortname = args[4]

bg = "white"

#DEBUG
#dataname = "MIP101_vs_MIP5FUR"
#datafile = "/projects/malachig/solexa/read_records/SI/Junctions_v53/temp_data.txt";
#outdir = "/projects/malachig/solexa/read_records/SI/Junctions_v53/";
#shortname = "Junctions"

#Define a fold-change cutoff level (not on a log scale)

#Also Select cutoff level and type for signficance p-values or q-values
#Type could be BH, Bonferroni, pvalue
#For no cutoff, specify 1

#5-FU DATA
#fc_cutoff= 2.25
#cutoff = 0.05
#cutoff_type = "BH"

#IRESSA RESPONSE DATA
fc_cutoff= 2.0
cutoff = 0.05
cutoff_type = "pvalue"

#Files to be created
text_si_values = paste(gsub(" +", "", dataname), "_SI_Values.txt", sep ="")
text_si_values_sorted = paste(gsub(" +", "", dataname), "_SI_Values_Sorted.txt", sep ="")
text_si_values_sorted_cutoff = paste(gsub(" +", "", dataname), "_SI_Values_Sorted_Cutoff.txt", sep ="")
text_si_values_sorted_cutoff_recip = paste(gsub(" +", "", dataname), "_SI_Values_Sorted_Cutoff_Recip.txt", sep ="")
text_si_values_sorted_cutoff_percent_SEQ_DE = paste(gsub(" +", "", dataname), "_SI_Values_Sorted_Cutoff_Percent_SEQ_DE.txt", sep ="")

svg_si_values_hist = paste(gsub(" +", "", dataname), "_SI_Values_Hist.svg", sep ="")
svg_si_values_cutoff_hist = paste(gsub(" +", "", dataname), "_SI_Values_Cutoff_Hist.svg", sep ="")

#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
setwd(outdir)

#Create normalized versions of each data column - Adjust scale of values in the larger library to that of the smaller library 
print ("Normalizing data to the smaller library")
libraryA_gene_size = sum(data[,"A_GENE_RAW"], na.rm=TRUE) 
libraryB_gene_size = sum(data[,"B_GENE_RAW"], na.rm=TRUE)
libraryA_seq_size = sum(data[,"A_SEQ_RAW"], na.rm=TRUE) 
libraryB_seq_size = sum(data[,"B_SEQ_RAW"], na.rm=TRUE)

larger_library_gene_size = max(libraryA_gene_size, libraryB_gene_size, na.rm=TRUE)
smaller_library_gene_size = min(libraryA_gene_size, libraryB_gene_size, na.rm=TRUE)
larger_library_seq_size = max(libraryA_seq_size, libraryB_seq_size, na.rm=TRUE)
smaller_library_seq_size = min(libraryA_seq_size, libraryB_seq_size, na.rm=TRUE)

data[,"A_GENE_Norm"] = data[,"A_GENE_RAW"]*(smaller_library_gene_size/libraryA_gene_size)
data[,"B_GENE_Norm"] = data[,"B_GENE_RAW"]*(smaller_library_gene_size/libraryB_gene_size)
data[,"A_SEQ_Norm"] = data[,"A_SEQ_RAW"]*(smaller_library_seq_size/libraryA_seq_size)
data[,"B_SEQ_Norm"] = data[,"B_SEQ_RAW"]*(smaller_library_seq_size/libraryB_seq_size)

z_gene = data[,c("A_GENE_Norm","B_GENE_Norm")]
z_seq = data[,c("A_SEQ_Norm","B_SEQ_Norm")]

#Calculate fold change and log2 difference for both GENES and SEQS
#Add '1' to all values to stabilize variance and prevent divion by 0
#Define function to calculate fold change - Fold changes will always be >= 1 (no change) and +ve versus -ve indicates the direction of the change
#Note that this will be library B - library A.  (e.g. MIP/5FU - MIP101)
fold.change=function(x){
  dataA = x[1]
  dataB = x[2]

  #Watch out for NA values
  if (is.na(dataA) | is.na(dataB)){
    result = NA;
  }else{
    #Otherwise calculate fold change value
    if (dataA >= dataB){
      result = (dataA+1)/(dataB+1)
    }else{
      result = ((dataB+1)/(dataA+1))*-1
    }
  }
  return(result)
}
print ("Calculating fold change values and log2 differences")
data[,"GENE_Fold_Change"] = apply(z_gene, 1, fold.change)
data[,"GENE_Log2_Diff"] = log2(data[,"A_GENE_Norm"]+1) - log2(data[,"B_GENE_Norm"]+1)
data[,"SEQ_Fold_Change"] = apply(z_seq, 1, fold.change)
data[,"SEQ_Log2_Diff"] = log2(data[,"A_SEQ_Norm"]+1) - log2(data[,"B_SEQ_Norm"]+1)

#Calculate an SI value for each expressed SEQ...

#Possible SI formula's (all produce equivalent values)
#SI = (log2(exon1_fur)-log2(gene1_fur)) - (log2(exon1_mip) - log2(gene1_mip))
#SI = log2(exon1_fur/gene1_fur) - log2(exon1_mip/gene1_mip)
#SI = log2((exon1_fur/gene1_fur) / (exon1_mip/gene1_mip))

#data[,"SI_1"] = (log2(data[,"B_SEQ_Norm"]+1) - log2(data[,"B_GENE_Norm"]+1)) - (log2(data[,"A_SEQ_Norm"]+1) - log2(data[,"A_GENE_Norm"]+1))
#data[,"SI_2"] = (log2((data[,"B_SEQ_Norm"]+1)/(data[,"B_GENE_Norm"]+1))) - (log2((data[,"A_SEQ_Norm"]+1)/(data[,"A_GENE_Norm"]+1)))
#data[,"SI_3"] = log2(((data[,"B_SEQ_Norm"]+1)/(data[,"B_GENE_Norm"]+1)) / ((data[,"A_SEQ_Norm"]+1)/(data[,"A_GENE_Norm"]+1)))
data[,"SI"] = log2(((data[,"A_SEQ_Norm"]+1)/(data[,"A_GENE_Norm"]+1)) / ((data[,"B_SEQ_Norm"]+1)/(data[,"B_GENE_Norm"]+1)))

#Determine which events exhibit a reciprocal shift ((Gene DE +ve AND Seq DE -ve) OR (Gene DE -ve AND Seq DE +ve))
data[,"Reciprocal"] = 0
i = which((data[,"GENE_Log2_Diff"] > 0 & data[,"SEQ_Log2_Diff"] < 0) | (data[,"GENE_Log2_Diff"] < 0 & data[,"SEQ_Log2_Diff"] > 0))  
data[i,"Reciprocal"] = 1

#Calculate the degree of reciprocity
#Reciprocity Score = abs((abs(GENE_DE) + abs(SEQ_DE)) / (abs(GENE_DE) -abs(SEQ_DE)))
data[,"Reciprocity"] = (abs(data[,"GENE_Log2_Diff"]) + abs(data[,"SEQ_Log2_Diff"])) / (abs(data[,"GENE_Log2_Diff"]) - abs(data[,"SEQ_Log2_Diff"]))



#Calculate the percentage of differential expression contributed by the SEQ relative to the GENE
#percent_SEQ_DE = (abs(SEQ_DE)/(abs(GENE_DE)+abs(SEQ_DE)))*100
data[,"percent_SEQ_Log2_DE"] = (abs(data[,"SEQ_Log2_Diff"]) / (abs(data[,"SEQ_Log2_Diff"]) + abs(data[,"GENE_Log2_Diff"])))*100


#Calculate p-values for the Feature DE
print ("Calculating fisher.test p-values")
fisher_data=cbind((data[,"A_SEQ_RAW"]), (sum(data[,"A_SEQ_RAW"])-data[,"A_SEQ_RAW"]), (data[,"B_SEQ_RAW"]), (sum(data[,"B_SEQ_RAW"])-data[,"B_SEQ_RAW"]))

#Define function to calculate fisher exact statistic for each cluster/term combination
cont_table_fun=function(x){
  cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
  fisher_test_result=fisher.test(cont_table, alternative = "two.sided")
  result=matrix(c(fisher_test_result$p.value, fisher_test_result$estimate[[1]]), nr=1)
  return(result)
}

#Apply the fisher function to all rows of the data for the appropriate columns
fisher_results=t(apply(fisher_data, 1, cont_table_fun))

#Give some column names to the output
colnames(fisher_results)=c("pvalue", "odds_ratio")

#Correct p-values
print ("Correcting p-value for multiple testing")
fisher_pvalues=as.numeric(fisher_results[,"pvalue"])
fisher_pvalues_adj=mt.rawp2adjp(fisher_pvalues, proc=c("Bonferroni","BH"))
fisher_pvalues_adj_orig_order=fisher_pvalues_adj$adjp[order(fisher_pvalues_adj$index),]

#Add these corrected p-values to the fisher_results data.frame
#fisher_results=cbind(fisher_results, format(fisher_pvalues_adj_orig_order[,2:3], scientific=TRUE))
fisher_results=cbind(fisher_results, fisher_pvalues_adj_orig_order[,2:3])
data = cbind(data[,], fisher_results)


#Write to file the log2 DE, foldchange, and SI values for every input element (gene, exon, junction, etc.) - also write a version sorted it on SI value
print ("Writing fold-change output text files")
output=cbind(data)
write.table (output, sep="\t", file=text_si_values, quote=FALSE, row.names=FALSE)
i = order(abs(output[,"SI"]), decreasing=TRUE, na.last=TRUE)
write.table (output[i,], sep="\t", file=text_si_values_sorted, quote=FALSE, row.names=FALSE)

#Now write a filtered list that passes the cutoffs
#CUTOFF DATA ONLY - NOTE that to be considered of interest BOTH the SI AND the SEQ DE (DE of junction, exon, etc.) must exceed the cutoff
i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
output=cbind(data[i,])
i = order(abs(output[,"SI"]), decreasing=TRUE, na.last=TRUE)
write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff, quote=FALSE, row.names=FALSE)

#Now write a filtered list that passes the cutoff AND involves reciprocal DE events (Gene up and Seq down or vice versa)
#Sort this list according to the degree of reciprocity
i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,"Reciprocal"] > 0 & data[,cutoff_type] < cutoff)
output=cbind(data[i,])
i = order(abs(output[,"Reciprocity"]), decreasing=TRUE, na.last=TRUE)
write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_recip, quote=FALSE, row.names=FALSE)

#Now write a filtered list that passes the cutoff
#Sort this list according to percent of DE contributed by the SEQ as apposed to the GENE
i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
output=cbind(data[i,])
i = order(abs(output[,"percent_SEQ_Log2_DE"]), decreasing=TRUE, na.last=TRUE)
write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_percent_SEQ_DE, quote=FALSE, row.names=FALSE)

#Now generate two histograms, one showing the distribution of DE values for all elements, the other for only those meeting a cutoff
print ("Generating histograms to display distribution of SI values")
data1 = data[,"SI"]
units = "Log2 SI"
string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
print(string)

#ALL DATA
if (length(data1)){
  main_title = paste("Splicing Index values - ", shortname, " (", units, ")", sep="");
  y_label = paste("Frequency (n = ", length(data1), ")", sep="")
  x_label = paste(dataname, " (", units, ")", sep="")

  Cairo(width=750, height=750, file=svg_si_values_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = data1, col="plum", main=main_title, xlab=x_label, ylab=y_label,
       col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks = 100)
  N = length(data1)
  min = min(data1)
  lower_cutoff = (log2(fc_cutoff))*-1
  upper_cutoff = log2(fc_cutoff)
  max = max(data1)
  abline(v=lower_cutoff, col="red", lty=2, lwd=2)
  abline(v=upper_cutoff, col="red", lty=2, lwd=2)
  legend_text = c(paste("N = ", N, sep=""),
                  paste("Min. = ", round(min, digits=1), " ", units, sep=""),
                  paste("Lower cutoff = ", round(lower_cutoff, digits=1), " ", sep=""),
                  paste("Upper cutoff = ", round(upper_cutoff, digits=1), " ", sep=""),
                  paste("Max = ", round(max, digits=1), " ", units, sep=""))
  legend_location = "topright"
  legend(legend_location, legend=legend_text, lty=c(0,0,2,2,0), lwd=2, col=c("black","black","red","red","black"))
  zz=dev.off();
}

#CUTOFF DATA ONLY - NOTE that to be considered of interest BOTH the SI AND the SEQ DE (DE of junction, exon, etc.) must exceed the cutoff
i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
data1 = data[i,"SI"]

if (length(data1)){
  units = "Log2 SI"
  string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
  print(string)

  main_title = paste("Splicing Index values - ", shortname, " (", units, ")", sep="");
  y_label = paste("Frequency (n = ", length(data1), ")", sep="")
  x_label = paste(dataname, " (", units, ")", " (>=" , round(log2(fc_cutoff), digits=1), ")", sep="")

  Cairo(width=750, height=750, file=svg_si_values_cutoff_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)

  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = data1, col="purple", main=main_title, xlab=x_label, ylab=y_label,
       col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks = 100)
  N = length(data1)
  min = min(data1)
  lower_cutoff = (log2(fc_cutoff))*-1
  upper_cutoff = log2(fc_cutoff)
  max = max(data1)
  abline(v=lower_cutoff, col="red", lty=2, lwd=2)
  abline(v=upper_cutoff, col="red", lty=2, lwd=2)
  legend_text = c(paste("N = ", N, sep=""),
                  paste("Min. = ", round(min, digits=1), " ", units, sep=""),
                  paste("Lower cutoff = ", round(lower_cutoff, digits=1), " ", sep=""),
                  paste("Upper cutoff = ", round(upper_cutoff, digits=1), " ", sep=""),
                  paste("Max = ", round(max, digits=1), " ", units, sep=""))
  legend_location = "topright"
  legend(legend_location, legend=legend_text, lty=c(0,0,2,2,0), lwd=2, col=c("black","black","red","red","black"))
  zz=dev.off();
}

#Note sure how to get p-values for SI data.
#Fisher's Exact does not work/make sense for this kind of data...

#Calculate Gene normalized SEQ expression values and apply a Fisher's exact test to generate a p-value for differences between the two libraries...
#data[,"A_GNE"] = (data[,"A_SEQ_RAW"]+1) / (data[,"A_GENE_RAW"]+1)
#data[,"B_GNE"] = (data[,"B_SEQ_RAW"]+1) / (data[,"B_GENE_RAW"]+1)
#print ("Calculating fisher.test p-values")
#fisher_data=cbind((data[,"B_GNE"]), (sum(data[,"B_GNE"])-data[,"B_GNE"]), (data[,"A_GNE"]), (sum(data[,"A_GNE"])-data[,"A_GNE"]))

#Apply a Fisher's exact test to the values WITHIN each ROW of data... (i.e. A_GENE_RAW, A_SEQ_RAW, B_GENE_RAW, B_SEQ_RAW)
#This is probably not valid...
#print ("Calculating fisher.test p-values")
#fisher_data=cbind((data[,"A_SEQ_RAW"]), data[,"A_GENE_RAW"], (data[,"B_SEQ_RAW"]), data[,"B_GENE_RAW"])

#Define function to calculate fisher exact statistic for each cluster/term combination
#cont_table_fun=function(x){
#  cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
#  fisher_test_result=fisher.test(cont_table, alternative = "two.sided")
#  result=matrix(c(fisher_test_result$p.value, fisher_test_result$estimate[[1]]), nr=1)
#  return(result)
#}

#Apply the fisher function to all rows of the data for the appropriate columns
#fisher_results=t(apply(fisher_data, 1, cont_table_fun))

#Give some column names to the output
#colnames(fisher_results)=c("pvalue", "odds_ratio")

#Correct p-values
#print ("Correcting p-value for multiple testing")
#fisher_pvalues=as.numeric(fisher_results[,"pvalue"])
#fisher_pvalues_adj=mt.rawp2adjp(fisher_pvalues, proc=c("Bonferroni","BH"))
#fisher_pvalues_adj_orig_order=fisher_pvalues_adj$adjp[order(fisher_pvalues_adj$index),]

#Add these corrected p-value to the fisher_results data.frame
#fisher_results=cbind(fisher_results, format(fisher_pvalues_adj_orig_order[,2:3], scientific=TRUE))
#fisher_results=cbind(fisher_results, fisher_pvalues_adj_orig_order[,2:3])
#final_results = cbind(data, fisher_results)































