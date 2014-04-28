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
#[4] genedatafile
#[5] geneexpressedfile
#[6] outdir
#[7] data_type_name
#[8] groupMemberLists

#Example usage:
#/scratch/obig/ALEXA/alexa_seq/R_bin/calculateGroupDifferentialSplicing_SVG.R Basal_vs_Luminal_Transcript /scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL/Transcripts_v53/Transcript_Basal_vs_Luminal_temp_data_GSI.txt /scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL/Transcripts_v53/Transcript_Basal_vs_Luminal_temp_expressed_GSI.txt /scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL/Transcripts_v53/Transcript_Basal_vs_Luminal_temp_gene_data_GSI.txt /scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL/Transcripts_v53/Transcript_Basal_vs_Luminal_temp_gene_expressed_GSI.txt /scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL/Transcripts_v53/ "Transcript" "Basal:HCC1806,HCC1143,SUM149PT,HCC1937,HCC1187,HCC1954,HCC3153,HCC1569,HCC70,HCC1500;Luminal:X600MPE,ZR75B,HCC1428,HCC202,MCF7,CAMA1,ZR7530,BT474,HCC1419,MDAMB13v1"

#NOTES:
#The input file must contain: a unique ID (gene, exon, junction), EnsEMBL Gene ID, Gene Name, and columns of expression data to be compared
#These expression values (from matrix files) should be normalized raw expression values
#The output file will contain these two data columns, normalized data values, the log2 difference, fold change, pvalue, and corrected pvalues (qvalue)

args = (commandArgs(TRUE))
dataname = args[1];
datafile = args[2];
expressedfile = args[3];
genedatafile = args[4];
geneexpressedfile = args[5];
outdir = args[6];
shortname = args[7];
groupMemberLists = args[8];

#DEBUG
#dataname = "Null_Female_vs_Null_Male_SilentIntronRegion"
#datafile = "/gscmnt/gc2142/techd/analysis/alexa_seq/figures_and_stats/SI/NF1/Introns_v54/SilentIntronRegion_Null_Female_vs_Null_Male_temp_data_GSI.txt"
#expressedfile = "/gscmnt/gc2142/techd/analysis/alexa_seq/figures_and_stats/SI/NF1/Introns_v54/SilentIntronRegion_Null_Female_vs_Null_Male_temp_expressed_GSI.txt"
#genedatafile = "/gscmnt/gc2142/techd/analysis/alexa_seq/figures_and_stats/SI/NF1/Introns_v54/SilentIntronRegion_Null_Female_vs_Null_Male_temp_gene_data_GSI.txt"
#geneexpressedfile = "/gscmnt/gc2142/techd/analysis/alexa_seq/figures_and_stats/SI/NF1/Introns_v54/SilentIntronRegion_Null_Female_vs_Null_Male_temp_gene_expressed_GSI.txt"
#outdir = "/gscmnt/gc2142/techd/analysis/alexa_seq/figures_and_stats/SI/NF1/Introns_v54/"
#shortname = "SilentIntronRegion" 
#groupMemberLists = "Null_Female:MM5298,MM5299,MM5300;Null_Male:MM5295,MM5296,MM5297"

#Determine the number of 'significant' DE features (qvalue < 0.05 and absolute fold-change > 1.5) 
fc_cutoff= 2
cutoff = 0.05
#cutoff_type = "BH"
cutoff_type = "rawp"
pe_thresh = 0.2 #Minimum percent libraries "expressed"
cov_min = 0.7 #Minimum coefficient of variation
cov_max = 10 #Maximum cov
bg = "white"

#Files to be created
text_si_values = paste(gsub(" +", "", dataname), "_GSI_Values.txt", sep ="")
text_si_values_sorted = paste(gsub(" +", "", dataname), "_GSI_Values_Sorted.txt", sep ="")
text_si_values_sorted_cutoff = paste(gsub(" +", "", dataname), "_GSI_Values_Sorted_Cutoff.txt", sep ="")
text_si_values_sorted_cutoff_recip = paste(gsub(" +", "", dataname), "_GSI_Values_Sorted_Cutoff_Recip.txt", sep ="")
text_si_values_sorted_cutoff_percent_SEQ_DE = paste(gsub(" +", "", dataname), "_GSI_Values_Sorted_Cutoff_Percent_SEQ_DE.txt", sep ="")
text_si_values_sorted_cutoff_qvalue = paste(gsub(" +", "", dataname), "_GSI_Values_Sorted_Cutoff_qvalue.txt", sep ="")

svg_si_values_hist = paste(gsub(" +", "", dataname), "_GSI_Values_Hist.svg", sep ="")
svg_si_values_cutoff_hist = paste(gsub(" +", "", dataname), "_GSI_Values_Cutoff_Hist.svg", sep ="")

#Read in data (expecting a tab-delimited file with header line and rownames)
raw_data=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
exp_status=read.table(expressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
gene_raw_data=read.table(genedatafile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))
gene_exp_status=read.table(geneexpressedfile, header = TRUE, na.strings = "NA", sep="\t", as.is=c(1:4))

#If for some reason a feature could not be assigned to a gene ID, remove it from the candidates...
not_na=which(!is.na(raw_data[,"EnsEMBL_Gene_ID"]))
raw_data=raw_data[not_na,]

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

#Define a percent expressed (PE) function and filter out features with less than minimum
pe_fun=function(x){
 pe=sum(x)/length(x)
 return(pe)
}

#Apply PE fun to feature-level (or both feature- and gene-level?)
w=exp_status[,libs]
pe_data=apply(w, 1, pe_fun)
passed_pe = which(pe_data >= pe_thresh)
data=raw_data[passed_pe,]

#w2=gene_exp_status[,libs]
#pe_data2=apply(w2, 1, pe_fun)
#passed_pe2 = which(pe_data2 >= pe_thresh)
#gene_data=gene_data[passed_pe2,]

#Define function for coefficient of variation (sd/mean) and filter out features not within min/max
cov_fun=function(x){
  cov=sd(x)/mean(x)
  return(cov)
}

#Apply COV function to feature expression only (cases where gene is invariant can still be interesting)
y=data[,libs]
cov_data=apply(y, 1, cov_fun)
passed_cov = which(cov_data < cov_max & cov_data > cov_min)
data=data[passed_cov,]

#Get corresponding gene expression data for feature/genes that survived filtering
feat_genes=data[,"EnsEMBL_Gene_ID"]
gene_data=gene_raw_data
genes=gene_data[,"EnsEMBL_Gene_ID"]
rownames(gene_data)=genes
gene_data=gene_data[feat_genes,]

#Grab just library data for genes/features remaining after filtering
z=data[,libs] #feature data
z2=gene_data[,libs] #gene data


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

#Define a function to calculate SI value
SI_fun=function(x){
  A_SEQ_Norm = mean(as.numeric(z[x,group_libs[[1]]]))
  A_GENE_Norm = mean(as.numeric(z2[x,group_libs[[1]]]))
  B_SEQ_Norm = mean(as.numeric(z[x,group_libs[[2]]]))
  B_GENE_Norm = mean(as.numeric(z2[x,group_libs[[2]]]))
  SI = log2( ((A_SEQ_Norm+1)/(A_GENE_Norm+1)) / ((B_SEQ_Norm+1)/(B_GENE_Norm+1)) )
  return(SI)
}

#Define function to calculate ANOVA statistic (NOTE: Currently not using interaction p-value "gene_expression:lib_group". Something to consider)
ANCOVA_fun=function(x){
  xdf=data.frame(t(z[x,]), t(z2[x,]), lib_group)
  names(xdf)=c("feat_expression", "gene_expression", "lib_group")
  ANCOVA_result=aov(feat_expression~gene_expression+lib_group, data=xdf)
  result=summary(ANCOVA_result)[[1]]["lib_group","Pr(>F)"]
  return(result)
}

#Grab indexes of all remaining data (use to pull out corresponding feature/gene from two datafiles)
#index=as.matrix(rownames(z))
index=as.matrix(1:length(rownames(z)))

#If number of groups is 2 calculate FC and log2 Diff on mean values, otherwise set to NA:
#Calculate fold change and log2 difference for both GENES and SEQS
#Then calculate SI value
if (length(groups)==2){
  print ("Calculating seq fold change values, log2 differences, and group means")
  data[,"SEQ_Fold_Change"] = apply(z, 1, fold.change)
  data[,"SEQ_Log2_Diff"] = apply(z, 1, log2.diff)
  data[,"SEQ_groupA_mean"] = apply(z[,group_libs[[1]]], 1, mean)
  data[,"SEQ_groupB_mean"] = apply(z[,group_libs[[2]]], 1, mean)
  print ("Calculating gene fold change values, log2 differences, and group means")
  data[,"GENE_Fold_Change"] = apply(z2, 1, fold.change)
  data[,"GENE_Log2_Diff"] = apply(z2, 1, log2.diff)
  data[,"GENE_groupA_mean"] = apply(z2[,group_libs[[1]]], 1, mean)
  data[,"GENE_groupB_mean"] = apply(z2[,group_libs[[2]]], 1, mean)
  print ("Calculating SI values")
  data[,"SI"] = apply(index, 1, SI_fun)

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

}else{
  data[,"SEQ_Fold_Change"] = "NA"
  data[,"SEQ_Log2_Diff"] = "NA"
  data[,"SEQ_groupA_mean"] = "NA"
  data[,"SEQ_groupB_mean"] = "NA"
  data[,"GENE_Fold_Change"] = "NA"
  data[,"GENE_Log2_Diff"] = "NA"
  data[,"GENE_groupA_mean"] = "NA"
  data[,"GENE_groupB_mean"] = "NA"
  data[,"SI"] = "NA"
  data[,"Reciprocal"] = "NA"
  data[,"Reciprocity"] = "NA"
  data[,"percent_SEQ_Log2_DE"] = "NA"
}

#Calculate ANCOVA p-values with feature expression as dependent variable and gene expression as a covariate
#Consider rank-transforming data to create non-parametric test
print ("Calculating ANCOVA p-values")

#Apply the ANCOVA function to all rows of the data for the appropriate columns
ANCOVA_results=apply(index, 1, ANCOVA_fun)

#Correct p-values
print ("Correcting p-value for multiple testing")
ANCOVA_pvalues_adj=mt.rawp2adjp(as.numeric(ANCOVA_results), proc="BH")
ANCOVA_pvalues_adj_orig_order=ANCOVA_pvalues_adj$adjp[order(ANCOVA_pvalues_adj$index),]
data=cbind(data[,1:4],data[,c("SEQ_groupA_mean","SEQ_groupB_mean","SEQ_Fold_Change","SEQ_Log2_Diff","GENE_groupA_mean","GENE_groupB_mean","GENE_Fold_Change","GENE_Log2_Diff","SI","Reciprocal","Reciprocity","percent_SEQ_Log2_DE")], ANCOVA_pvalues_adj_orig_order)

#Write to file the log2 DE, foldchange, and SI values for every input element (gene, exon, junction, etc.) - also write a version sorted it on SI/pvalue
print ("Writing fold-change output text files")
output=cbind(data)
write.table (output, sep="\t", file=text_si_values, quote=FALSE, row.names=FALSE)

if (length(groups)==2){
  i = order(abs(output[,"SI"]), decreasing=TRUE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted, quote=FALSE, row.names=FALSE)
}else{
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted, quote=FALSE, row.names=FALSE)
}

#Now write a filtered list that passes the cutoffs
#CUTOFF DATA ONLY - NOTE that to be considered of interest BOTH the SI AND the SEQ DE (DE of junction, exon, etc.) must exceed the cutoff
if (length(groups)==2){
  i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"SI"]), decreasing=TRUE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff, quote=FALSE, row.names=FALSE)
}else{
  i = which(data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff, quote=FALSE, row.names=FALSE)
}

#Now write a filtered list that passes the cutoff AND involves reciprocal DE events (Gene up and Seq down or vice versa)
#Sort this list according to the degree of reciprocity
if (length(groups)==2){
  i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,"Reciprocal"] > 0 & data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"Reciprocity"]), decreasing=TRUE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_recip, quote=FALSE, row.names=FALSE)
}else{
  i = which(data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_recip, quote=FALSE, row.names=FALSE)
}

#Now write a filtered list that passes the cutoff
#Sort this list according to percent of DE contributed by the SEQ as apposed to the GENE
if (length(groups)==2){
  i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"percent_SEQ_Log2_DE"]), decreasing=TRUE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_percent_SEQ_DE, quote=FALSE, row.names=FALSE)
}else{
  i = which(data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_percent_SEQ_DE, quote=FALSE, row.names=FALSE)
}

#Now write a filtered list that passes the cutoff
#Sort this list according to qvalue
if (length(groups)==2){
  i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_qvalue, quote=FALSE, row.names=FALSE)
}else{
  i = which(data[,cutoff_type] < cutoff)
  output=cbind(data[i,])
  i = order(abs(output[,"rawp"]), decreasing=FALSE, na.last=TRUE)
  write.table (output[i,], sep="\t", file=text_si_values_sorted_cutoff_qvalue, quote=FALSE, row.names=FALSE)
}


#Now generate two histograms, one showing the distribution of DE values for all elements, the other for only those meeting a cutoff
print ("Generating histograms to display distribution of SI values")

#ALL DATA
if (length(groups)==2){
  data1 = data[,"SI"]
  units = "Log2 SI"
  string = paste ("PLOTTING DATA.  length:", length(data1), " min:", min(data1), " max:", max(data1))
  print(string)

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
if (length(groups)==2){
  i = which(abs(data[,"SI"]) > log2(fc_cutoff) & abs(data[,"SEQ_Log2_Diff"]) > log2(fc_cutoff) & data[,cutoff_type] < cutoff)
  data1 = data[i,"SI"]

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

