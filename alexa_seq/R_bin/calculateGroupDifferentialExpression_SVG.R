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
#[3] expressedfile
#[4] outdir
#[5] data_type_name
#[6] groupMemberLists

#Example usage:
#/scratch/obig/ALEXA/alexa_seq/R_bin/calculateGroupDifferentialExpression_SVG.R Basal_vs_Luminal_Gene /scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/Gene_Basal_vs_Luminal_temp_data.txt /scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/Gene_Basal_vs_Luminal_temp_expressed.txt /scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/ "Gene Basal:HCC1806,HCC1143,SUM149PT,HCC1937,HCC1187;Luminal:600MPE,ZR75B,HCC1428,HCC202,MCF7,CAMA1,ZR7530,BT474,HCC1419,MDAMB13v1"

#NOTES:
#The input file must contain: a unique ID (gene, exon, junction), EnsEMBL Gene ID, Gene Name, and columns of expression data to be compared
#These expression values (from matrix files) should be normalized raw expression values
#The output file will contain these two data columns, normalized data values, the log2 difference, fold change, pvalue, and corrected pvalues (qvalue)
args = (commandArgs(TRUE))
dataname = args[1];
datafile = args[2];
expressedfile = args[3];
outdir = args[4];
shortname = args[5];
groupMemberLists = args[6];

#Determine the number of 'significant' DE features (qvalue < 0.05 and absolute fold-change > 1.5) 
fc_cutoff= 2
qvalue_cutoff = 0.05
pvalue_cutoff = 0.05
pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 0.7 #Minimum coefficient of variation
cov_max = 10 #Maximum cov

bg = "white"

#DEBUG
#dataname = "Basal_vs_Luminal_Gene"
#datafile = "/scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/Gene_Basal_vs_Luminal_temp_data.txt"
#expressedfile = "/scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/Gene_Basal_vs_Luminal_temp_expressed.txt"
#outdir = "/scratch/obig/ALEXA/analysis/figures_and_stats/DE/BCCL/ENST_v53/"
#shortname = "Gene" 
#groupMemberLists = "Basal:HCC1806,HCC1143,SUM149PT,HCC1937,HCC1187;Luminal:X600MPE,ZR75B,HCC1428,HCC202,MCF7,CAMA1,ZR7530,BT474,HCC1419,MDAMB13v1"

#Files to be created
text_de_values = paste(gsub(" +", "", dataname), "_GDE_Values.txt", sep ="")
text_de_values_sorted = paste(gsub(" +", "", dataname), "_GDE_Values_Sorted.txt", sep ="")
text_de_sig_sorted_mtc = paste(gsub(" +", "", dataname), "_GDE_Values_Significant_MTC.txt", sep ="")
text_de_sig_sorted_nomtc = paste(gsub(" +", "", dataname), "_GDE_Values_Significant_NoMTC.txt", sep ="")
svg_de_values_hist = paste(gsub(" +", "", dataname), "_GDE_Values_Hist.svg", sep ="")
svg_de_sig_mtc_hist = paste(gsub(" +", "", dataname), "_GDE_Values_Significant_Hist_MTC.svg", sep ="")
svg_de_sig_nomtc_hist = paste(gsub(" +", "", dataname), "_GDE_Values_Significant_Hist_NoMTC.svg", sep ="")
jpeg_exp_scatter = paste(gsub(" +", "", dataname), "_Expression_Norm_Scatter.jpeg", sep ="")

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
exp_status=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
setwd(outdir)
header=colnames(data)

groups=strsplit(groupMemberLists, split=";")[[1]]
group_libs=vector("list", length=length(groups))
group_names=vector(length=length(groups))

libcount=0
for (i in 1:length(groups)){
  group_name=strsplit(groups[i],split=":")[[1]][1]
  libraries_str=strsplit(groups[i],split=":")[[1]][2]
  libraries=strsplit(libraries_str,split=",")[[1]]
  libcount=libcount+length(libraries)
  group_names[i]=group_name
  group_libs[[i]]=libraries
}

#Make vectors of libraries and their groups in the same order for look up purposes
libs=vector(length=libcount)
lib_group=vector(length=libcount)
k=0
for (i in 1:length(group_libs)){
  for (j in 1:length(group_libs[[i]])){
    k=k+1
    libs[k]=group_libs[[i]][j]
    lib_group[k]=group_names[i]
  }
}


#Define a percent expressed function and filter out features with less than minimum
w=exp_status[,libs]
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}
pe_data=apply(w, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data=raw_data[passed_pe,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
y=data[,libs]
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}
cov_data=apply(y, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data=data[passed_cov,]

#Grab just library data for genes remaining after filtering
z=data[,libs]

#Define function to calculate fold change - Fold changes will always be >= 1 (no change) and +ve versus -ve indicates the direction of the change
#Add '1' to all values to stabilize variance and prevent divion by 0
fold.change=function(x){
  dataA = x[group_libs[[1]]]
  dataB = x[group_libs[[2]]]
  mean_dataA = mean(dataA)
  mean_dataB = mean(dataB)

  #calculate fold change value
  if (mean_dataA >= mean_dataB){
      result = (mean_dataA+1)/(mean_dataB+1)
    }else{
      result = ((mean_dataB+1)/(mean_dataA+1))*-1
    }
  return(result)
}
#Define a function to calculate log2 difference
log2.diff=function(x){
  dataA = x[group_libs[[1]]]
  dataB = x[group_libs[[2]]]
  mean_dataA = mean(dataA)
  mean_dataB = mean(dataB)
  result = log2(mean_dataA+1) - log2(mean_dataB+1)
  return(result)
}

#If number of groups is 2 calculate FC and log2 Diff on mean values, otherwise set to NA:
if (length(groups)==2){
 print ("Calculating fold change values and log2 differences")
  data[,"Fold_Change"] = apply(z, 1, fold.change)
  data[,"Log2_Diff"] = apply(z, 1, log2.diff)
  data[,"groupA_mean"] = apply(z[,group_libs[[1]]], 1, mean)
  data[,"groupB_mean"] = apply(z[,group_libs[[2]]], 1, mean)
}else{
  data[,"Fold_Change"] = "NA"
  data[,"Log2_Diff"] = "NA"
  data[,"groupA_mean"] = "NA"
  data[,"groupB_mean"] = "NA"
}

#Calculate ANOVA p-values 
#Consider rank-transforming data to create non-parametric test
#suppressWarnings?
print ("Calculating ANOVA p-values")

#Define function to calculate ANOVA statistic
ANOVA_fun=function(x){
  xdf=data.frame(x,lib_group)
  names(xdf)=c("expression", "lib_group")
  ANOVA_result=aov(expression~lib_group, data=xdf)
  result=summary(ANOVA_result)[[1]]["lib_group","Pr(>F)"]
  return(result)
}

#Apply the ANOVA function to all rows of the data for the appropriate columns
ANOVA_results=apply(z, 1, ANOVA_fun)


#Correct p-values
print ("Correcting p-value for multiple testing")
ANOVA_pvalues_adj=mt.rawp2adjp(as.numeric(ANOVA_results), proc="BH")
ANOVA_pvalues_adj_orig_order=ANOVA_pvalues_adj$adjp[order(ANOVA_pvalues_adj$index),]

#Create final results table for printing to output
final_results = cbind(data[,1:4],data[,c("Fold_Change","Log2_Diff","groupA_mean","groupB_mean")], ANOVA_pvalues_adj_orig_order)


#Print out only the significant results (sorted by FC or p-value) BEFORE multiple testing correction
#If number of groups is 2 include FC as a filtering criteria, otherwise just use p-value:
print ("Writing output text files")
if (length(groups)==2){
  sig_elements_nomtc = which(abs(final_results[,"Fold_Change"]) >= fc_cutoff & final_results[,"rawp"] < pvalue_cutoff)
  sig_data_nomtc = final_results[sig_elements_nomtc,]
  i = order(abs(sig_data_nomtc[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
  write.table (sig_data_nomtc[i,], sep="\t", file=text_de_sig_sorted_nomtc, quote=FALSE, row.names=FALSE)
}else{
  sig_elements_nomtc = which(final_results[,"rawp"] < pvalue_cutoff)
  sig_data_nomtc = final_results[sig_elements_nomtc,]
  i = order(sig_data_nomtc[,"rawp"], decreasing=FALSE, na.last=TRUE)
  write.table (sig_data_nomtc[i,], sep="\t", file=text_de_sig_sorted_nomtc, quote=FALSE, row.names=FALSE)
}

#Print out only the significant results (sorted by FC or p-value) AFTER multiple testing correction
if (length(groups)==2){
  sig_elements_mtc = which(abs(final_results[,"Fold_Change"]) >= fc_cutoff & final_results[,"BH"] < qvalue_cutoff)
  sig_data_mtc = final_results[sig_elements_mtc,]
  i = order(abs(sig_data_mtc[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
  write.table (sig_data_mtc[i,], sep="\t", file=text_de_sig_sorted_mtc, quote=FALSE, row.names=FALSE)
}else{
  sig_elements_mtc = which(final_results[,"BH"] < qvalue_cutoff)
  sig_data_mtc = final_results[sig_elements_mtc,]
  i = order(abs(sig_data_mtc[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (sig_data_mtc[i,], sep="\t", file=text_de_sig_sorted_mtc, quote=FALSE, row.names=FALSE)
}

#Write to file the log2 DE and foldchange values for every input element (gene, exon, junction, etc.) - also write a sorted version of this file
write.table (final_results, sep="\t", file=text_de_values, quote=FALSE, row.names=FALSE)
if (length(groups)==2){
  i = order(abs(final_results[,"Fold_Change"]), decreasing=TRUE, na.last=TRUE)
  write.table (final_results[i,], sep="\t", file=text_de_values_sorted, quote=FALSE, row.names=FALSE)
}else{
  i = order(abs(final_results[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (final_results[i,], sep="\t", file=text_de_values_sorted, quote=FALSE, row.names=FALSE)
}

#Create a plot of normalized expression values from library A vs. library B - simple scatter plot
if (length(groups)==2){
  x1 = log2(final_results[,"groupA_mean"]+1)
  y1 = log2(final_results[,"groupB_mean"]+1)
  x2 = log2(final_results[sig_elements_mtc,"groupA_mean"]+1)
  y2 = log2(final_results[sig_elements_mtc,"groupB_mean"]+1)
  xlab=paste("group A (", group_names[1], ") expression (log2[mean expression+1])", sep="")
  ylab=paste("group B (", group_names[2], ") expression (log2[mean expression+1])", sep="")
  CairoJPEG(filename=jpeg_exp_scatter, width=750, height=750, pointsize=12, quality=100, bg=bg)
  par(bg=bg, font.main = 2, font.lab = 2)
  plot(x=x1, y=y1, col="blue", pch=20, xlab=xlab, ylab=ylab, main="Correlation between expression values from groups A and B",
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

#Generate histogram showing the distribution of DE values
#All data
if (length(groups)==2){
  data1 = final_results[,"Log2_Diff"]
  if (length(data1) > 0){
    units = "mean Log2 DE"
    string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
    print(string)

    main_title = paste("DE values - ", shortname, " (", units, ")", sep="");
    y_label = paste("Frequency (n = ", length(data1), ")", sep="")
    x_label = paste(dataname, " (", units, ")", sep="")

    Cairo(width=750, height=750, file=svg_de_values_hist, type="svg", pointsize=12, bg="white", units="px", dpi=72)
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
}


#Sig Data after MTC
if (length(groups)==2){
  data1 = sig_data_mtc[,"Log2_Diff"]
  if (length(data1) > 0){
    units = "mean Log2 DE"
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
}


#Sig Data No MTC
if (length(groups)==2){
  data1 = sig_data_nomtc[,"Log2_Diff"]
  if (length(data1) > 0){
    units = "mean Log2 DE"
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
}


