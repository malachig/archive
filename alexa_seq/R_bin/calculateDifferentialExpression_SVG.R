#!/usr/bin/env Rscript
#Written by Malachi Griffith and Obi Griffith
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
library(geneplotter)
library(RColorBrewer)
library(gcrma)

#Options:
#[1] dataname
#[2] datafile
#[3] outdir
#[4] data_type_name
#[5] libraryA_name
#[6] libraryB_name

#Example usage:
#/home/malachig/svn/solexa_analysis/R_bin/calculateDifferentialExpression_SVG.R Mip101_vs_Mip5FuR_Gene /projects/malachig/solexa/figures_and_stats/DE/5FU/ENST_v53/Gene_temp_data.txt /projects/malachig/solexa/figures_and_stats/DE/5FU/ENST_v53/ "Gene" "MIP/5FU" "MIP101"


#NOTES:
#The input file must contain: a unique ID (gene, exon, junction), EnsEMBL Gene ID, Gene Name, and two columns of expression data (library A and B) to be compared
#These expression values should be raw expression values (normalization to the smaller library will be done by this script)
#The output file will contain these two data columns, normalized data values, the log2 difference, fold change, pvalue, and corrected pvalues (qvalue)
args = (commandArgs(TRUE))
dataname = args[1];
datafile = args[2];
outdir = args[3];
shortname = args[4];
libraryA_name = args[5];
libraryB_name = args[6];

#Determine the number of 'significant' DE features (qvalue < 0.05 and absolute fold-change > 1.5) 
fc_cutoff= 2
qvalue_cutoff = 0.05
pvalue_cutoff = 0.05

bg = "white"

#DEBUG
#dataname = "MIP101_vs_MIP5FUR"
#datafile = "/projects/malachig/solexa/read_records/DE/ENST_v53/temp_data.txt";
#outdir = "/projects/malachig/solexa/read_records/DE/ENST_v53/";
#shortname = "Genes"

#Files to be created
text_de_values = paste(gsub(" +", "", dataname), "_DE_Values.txt", sep ="")
text_de_values_sorted = paste(gsub(" +", "", dataname), "_DE_Values_Sorted.txt", sep ="")
text_de_sig_sorted_mtc = paste(gsub(" +", "", dataname), "_DE_Values_Significant_MTC.txt", sep ="")
text_de_sig_sorted_nomtc = paste(gsub(" +", "", dataname), "_DE_Values_Significant_NoMTC.txt", sep ="")
svg_de_values_hist = paste(gsub(" +", "", dataname), "_DE_Values_Hist.svg", sep ="")
svg_de_sig_mtc_hist = paste(gsub(" +", "", dataname), "_DE_Values_Significant_Hist_MTC.svg", sep ="")
svg_de_sig_nomtc_hist = paste(gsub(" +", "", dataname), "_DE_Values_Significant_Hist_NoMTC.svg", sep ="")
jpeg_exp_scatter = paste(gsub(" +", "", dataname), "_Expression_Norm_Scatter.jpeg", sep ="")
jpeg_exp_smoothscatter = paste(gsub(" +", "", dataname), "_Expression_Norm_SmoothScatter.jpeg", sep ="")


#Read in data (expecting a tab-delimited file with header line and rownames)
data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:3))
setwd(outdir)
old_names = names(data)
names(data) = c(old_names[1], old_names[2], old_names[3], "A_Average_Coverage_RAW", "A_Expressed", "B_Average_Coverage_RAW", "B_Expressed")

#Create normalized versions of each data column - Adjust scale of values in the larger library to that of the smaller library 
print ("Normalizing data to the smaller library")
libraryA_size = sum(data[,"A_Average_Coverage_RAW"], na.rm=TRUE) 
libraryB_size = sum(data[,"B_Average_Coverage_RAW"], na.rm=TRUE)
larger_library_size = max(libraryA_size, libraryB_size, na.rm=TRUE)
smaller_library_size = min(libraryA_size, libraryB_size, na.rm=TRUE)

data[,"A_Norm"] = data[,"A_Average_Coverage_RAW"]*(smaller_library_size/libraryA_size)
data[,"B_Norm"] = data[,"B_Average_Coverage_RAW"]*(smaller_library_size/libraryB_size)
z = data[,c("A_Norm","B_Norm")]

#Calculate fold change and log2 difference for all genes/exons/junctions
#Add '1' to all values to stabilize variance and prevent divion by 0
#Define function to calculate fold change - Fold changes will always be >= 1 (no change) and +ve versus -ve indicates the direction of the change
fold.change=function(x){
  dataA = x["A_Norm"]
  dataB = x["B_Norm"]

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
data[,"Fold_Change"] = apply(z, 1, fold.change)
data[,"Log2_Diff"] = log2(data[,"A_Norm"]+1) - log2(data[,"B_Norm"]+1)


#Create an array with values for each gene to be analysed by fisher.test
#Use the data values (read count, average coverage, etc.) for each library found in columns 2 and 3
#For the two libraries we will compare the average coverage for geneA versus the sum of average coverage for all other genes/exons/junctions
initial_summary = paste("Sum of data for library A = ", round(libraryA_size, digits=1), "     Sum of data for library B = ", round(libraryB_size, digits=1), sep="")
print(initial_summary)


#Now define the subset of elements with evidence for expression in at least one of the two conditions
#Calculate p-values for only this subset of the data
ei = which(data[,"A_Expressed"] == 1 | data[,"B_Expressed"] == 1)

print ("Calculating fisher.test p-values")
fisher_data=cbind((data[ei,"A_Average_Coverage_RAW"]), (sum(data[ei,"A_Average_Coverage_RAW"])-data[ei,"A_Average_Coverage_RAW"]), (data[ei,"B_Average_Coverage_RAW"]), (sum(data[ei,"B_Average_Coverage_RAW"])-data[ei,"B_Average_Coverage_RAW"]))

#Define function to calculate fisher exact statistic for each cluster/term combination
cont_table_fun=function(x){
  cont_table=matrix(x,nr = 2,dimnames=list(Expression=c("geneA", "Others"),Library=c("Sample1", "Sample2")))
  fisher_test_result=suppressWarnings(fisher.test(cont_table, alternative = "two.sided"))
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
final_results = cbind(data[ei,], fisher_results)

#Print out only the significant results (sorted by fold-change) BEFORE multiple testing correction
sig_elements_nomtc = which(abs(final_results[,"Fold_Change"]) >= fc_cutoff & final_results[,"pvalue"] < pvalue_cutoff)
sig_data_nomtc = final_results[sig_elements_nomtc,]
i = order(abs(sig_data_nomtc[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
write.table (sig_data_nomtc[i,], sep="\t", file=text_de_sig_sorted_nomtc, quote=FALSE, row.names=FALSE)

#Print out only the significant results (sorted by fold-change) AFTER multiple testing correction
sig_elements_mtc = which(abs(final_results[,"Fold_Change"]) >= fc_cutoff & final_results[,"BH"] < qvalue_cutoff)
sig_data_mtc = final_results[sig_elements_mtc,]
i = order(abs(sig_data_mtc[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
write.table (sig_data_mtc[i,], sep="\t", file=text_de_sig_sorted_mtc, quote=FALSE, row.names=FALSE)

#Write to file the log2 DE and foldchange values for every input element (gene, exon, junction, etc.) - also write a sorted version of this file
print ("Writing fold-change output text files")
data[,"BH"] = "NA";
data[ei,"BH"] = fisher_results[,"BH"]
output=cbind(data)
write.table (output, sep="\t", file=text_de_values, quote=FALSE, row.names=FALSE)
i = order(abs(output[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
write.table (output[i,], sep="\t", file=text_de_values_sorted, quote=FALSE, row.names=FALSE)

#Create a plot of normalized expression values from library A vs. library B - simple scatter plot
x1 = log2(final_results[,"A_Norm"]+1)
y1 = log2(final_results[,"B_Norm"]+1)
x2 = log2(final_results[sig_elements_mtc,"A_Norm"]+1)
y2 = log2(final_results[sig_elements_mtc,"B_Norm"]+1)

if (length(x1)){
  xlab=paste("Library A (", libraryA_name, ") expression (log2[expression+1])", sep="")
  ylab=paste("Library B (", libraryB_name, ") expression (log2[expression+1])", sep="")
  CairoJPEG(filename=jpeg_exp_scatter, width=750, height=750, pointsize=12, quality=100, bg=bg)
  par(bg=bg, font.main = 2, font.lab = 2)
  plot(x=x1, y=y1, col="blue", pch=20, xlab=xlab, ylab=ylab, main="Correlation between expression values from libraries A and B",
       col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
  points(x=x2, y=y2, col="magenta", pch=20)
  lmfit = lm(y1~x1)
  lines(loess.smooth(x1, y1, span = 2/3, degree = 1, family = "gaussian", evaluation = 100), col="black", lty=2, lwd=2)
  abline(lmfit$coefficients, col="red", lwd=2, lty=2)
  corr = cor(x1,y1, method="spearman")
  legend_text = c(paste("Correlation = ", round(corr, digits=4), " (Spearman)", sep=""), "Differentially expressed", "NOT differentially expressed", "Loess fit")
  legend_location = "topleft"
  legend(legend_location, legend=legend_text, lty=c(2,1,1,2), lwd=c(2,0,0,2), pch=c(NA,20,20,NA), col=c("red","magenta","blue","black"),cex=1.3, pt.cex=1.3)
  zz=dev.off()
}

#Create a plot of normalized expression values from library A vs. library B - density scatter plot
if (length(x1) > 0){
  CairoJPEG(filename=jpeg_exp_smoothscatter, width=750, height=750, pointsize=12, quality=100, bg=bg)
  colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  par(bg=bg, font.main = 2, font.lab = 2)
  smoothScatter(x1, y1, xlab=xlab, ylab=ylab, main="Correlation between expression values from libraries A and B", colramp=colors, nbin=275,
                col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
  lmfit = lm(y1~x1)
  lines(loess.smooth(x1, y1, span = 2/3, degree = 1, family = "gaussian", evaluation = 100), col="black", lty=2, lwd=2)
  abline(lmfit$coefficients, col="red", lwd=2, lty=2)
  corr = cor(x1,y1, method="spearman")
  legend_text = c(paste("Correlation = ", round(corr, digits=4), " (Spearman)", sep=""), "Loess fit")
  legend_location = "topleft"
  legend(legend_location, legend=legend_text, lty=c(2,2), lwd=c(2,2), col=c("red","black"),cex=1.3, pt.cex=1.3)
  zz=dev.off()
}

#Now generate two histograms, one showing the distribution of DE values for all elements, the other for only the significant elements
#All data
print ("Generating histograms to display distribution of DE values")
data1 = final_results[,"Log2_Diff"]
if (length(data1) > 0){
  units = "Log2 DE"
  string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
  print(string)

  main_title = paste("DE values - ", shortname, " values ", "(", units, ")", sep="");
  y_label = paste("Frequency (n = ", length(data1), ")", sep="")
  x_label = paste(dataname, " (", units, ")", sep="")

  Cairo(width=750, height=750, file=svg_de_values_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = data1, col="blue", main=main_title, xlab=x_label, ylab=y_label,
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

#Sig Data MTC
data1 = sig_data_mtc[,"Log2_Diff"]
if (length(data1) > 0){
  units = "Log2 DE"
  string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
  print(string)

  main_title = paste("Significant DE MTC values - ", shortname, " (", units, ")", sep="");
  y_label = paste("Frequency (n = ", length(data1), ")", sep="")
  x_label = paste(dataname, " (", units, ")", sep="")

  Cairo(width=750, height=750, file=svg_de_sig_mtc_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = data1, col="blue", main=main_title, xlab=x_label, ylab=y_label,
      col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
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

#Sig Data No MTC
data1 = sig_data_nomtc[,"Log2_Diff"]
if (length(data1) > 0){
  units = "Log2 DE"
  string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
  print(string)

  main_title = paste("Significant DE values - ", shortname, " (", units, ")", sep="");
  y_label = paste("Frequency (n = ", length(data1), ")", sep="")
  x_label = paste(dataname, " (", units, ")", sep="")

  Cairo(width=750, height=750, file=svg_de_sig_nomtc_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = data1, col="blue", main=main_title, xlab=x_label, ylab=y_label,
      col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
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




