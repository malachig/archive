#!/usr/bin/env Rscript
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes a tab delimited files of gene and exon summary data and summarizes basic stats
#These data represent paired end reads that have been mapped to an EnsEMBL transcriptome by BLAST

#Options:
#[1] datafile1
#[2] datafile2
#[3] results_dir1
#[4] results_dir2
#[5] ensembl_version
#[6] library_name

#data_file1 = "/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_GeneExpression_v53.txt"
#data_file2 = "/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_ExonRegionExpression_v53.txt"
#results_dir1 = "/projects/malachig/solexa/figures_and_stats/Generic/ENST_v53/"
#results_dir2 = "/projects/malachig/solexa/figures_and_stats/HS04391/ENST_v53/"
#ensembl_version = "53"
#library_name = "MIP101"


#Example usage:
#/home/malachig/svn/solexa_analysis/R_bin/summarizeGeneExonCoverage.R /projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_GeneExpression_v53.txt  /projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_ExonRegionExpression_v53.txt  /projects/malachig/solexa/figures_and_stats/Generic/ENST_v53/  /projects/malachig/solexa/figures_and_stats/HS04391/ENST_v53/  53  MIP101

args = (commandArgs(TRUE))
data_file1 = args[1]
data_file2 = args[2]
results_dir1 = args[3]
results_dir2 = args[4]
ensembl_version = args[5]
library_name = args[6]

print ("ARGS supplied are as follows:")
print (data_file1)
print (data_file2)
print (results_dir1)
print (results_dir2)
print (ensembl_version)
print (library_name)

#Generic stats 
print ("Loading data files ...")
gene_data = read.table(data_file1, header=T, quote="",sep="\t", comment.char="", na.strings='NA', as.is=c(2:8,17,30))
exon_data = read.table(data_file2, header=T, quote="",sep="\t", comment.char="", na.strings='NA')
setwd(results_dir1)

bg = "white"

print ("Data import (sizes):")
print (length(gene_data[,1]))
print (length(exon_data[,1]))


#A.) EXON SIZES
units = "bp"
divider = quantile(exon_data[,"Base_Count"], probs=0.95)

data = exon_data[,"Base_Count"]
title = paste("Distribution of exon sizes (EnsEMBL v", ensembl_version, ")", sep="")
x_label = paste("Exon Size (", units, ")", sep="")
y_label = paste("Exon Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonSizes.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))

legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()

lower = which(exon_data[,"Base_Count"] <= divider)
data = exon_data[lower,"Base_Count"]
x_label = paste("Exon Size (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("Exon Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonSizes_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()


#B.) NUMBER OF EXONS PER GENE
units = "exons"
data = table(exon_data[,"Gene_ID"])
divider = quantile(data, probs=0.95)
title = paste("Distribution of exons per gene (EnsEMBL v", ensembl_version, ")", sep="")

x_label = paste("Exons Per Gene (", units, ")", sep="")
y_label = paste("Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonsPerGene.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))

legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()

temp = table(exon_data[,"Gene_ID"])
lower = which(temp <= divider)
data = temp[lower]
title = paste("Distribution of exons per gene (EnsEMBL v", ensembl_version, ")", sep="")
x_label = paste("Exons Per Gene (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonsPerGene_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))

legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()


#C.) GENE EXONIC BASE COUNTS (I.E. GENE SIZE)
units = "bp"
data = gene_data[,"Base_Count"]
divider = quantile(data, probs=0.95)
title = paste("Distribution of gene sizes (EnsEMBL v", ensembl_version, ")", sep="")
x_label = paste("Gene Size (", units, ")", sep="")
y_label = paste("Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneSizes.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))

legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()

lower = which(gene_data[,"Base_Count"] <= divider)
data = gene_data[lower,"Base_Count"]
x_label = paste("Gene Size (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneSizes_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()


#D.) GENE TYPE  
data = gene_data[,"Gene_Type"]
x_label = paste("Gene Type (n = ", length(data), ")", sep="") 
title = paste("Distribution of gene types (EnsEMBL v", ensembl_version, ")", sep="")

tiff(file="GeneTypes.tiff", width=1000, height=1000, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
pie(x=table(data), xlab=x_label, col=rainbow(length(table(data))), main=title, col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0) 
dev.off()


#E.) GENE EVIDENCE
data = gene_data[,"Gene_Evidence"]
x_label = paste("Gene Evidence (n = ", length(data), ")", sep="") 
title = paste("Distribution of gene evidence codes (EnsEMBL v", ensembl_version, ")", sep="")

tiff(file="GeneEvidence.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
pie(x=table(data), xlab=x_label, col=rainbow(length(table(data))), main=title, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0) 
dev.off()


#F.) GENES PER CHROMOSOME
data = gene_data[,"Chromosome"]
x_label = paste("Genes per Chromosome (n = ", length(data), ")", sep="") 
title = paste("Distribution of genes across chromosomes (EnsEMBL v", ensembl_version, ")", sep="")

tiff(file="GeneChromosomes.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
pie(x=table(data), xlab=x_label, col=rainbow(length(table(data))), main=title, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0) 
dev.off()


#LIBRARY SPECIFIC SUMMARIES
setwd(results_dir2)


#G.) AVERAGE EXON COVERAGE RAW - 'expressed' exons only (coverage > 0)
summary(exon_data[,"Average_Coverage_RAW"])  #Dynamic range of exon expression levels is ~10^5
units = "X"
non_zero = which(exon_data[,"Average_Coverage_RAW"] > 0)
data = exon_data[non_zero,"Average_Coverage_RAW"]
divider = quantile(data, probs=0.95)
title = paste(library_name, " - Distribution of Average Exon Base Coverage", sep="")
x_label = paste("Mean Base Coverage (", units, ")", sep="")
y_label = paste("'Expressed' Exon Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonCoverage.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", round(min, digits=4), " ", units, sep=""),
                paste("Median = ", round(median, digits=1), " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", round(max, digits=1), " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()

lower = which(exon_data[,"Average_Coverage_RAW"] > 0 & exon_data[,"Average_Coverage_RAW"] <= divider)
data = exon_data[lower,"Average_Coverage_RAW"]
x_label = paste("Mean Base Coverage (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("'Expressed' Exon Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="ExonCoverage_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", round(min, digits=4), " ", units, sep=""),
                paste("Median = ", round(median, digits=1), " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", round(max, digits=1), " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()


#H.) GENE READ COUNT
summary(gene_data[,"Read_Count"])  #Dynamic range of gene read counts is : 1 read - 2.1 million reads (for H19 in MIP101 cells...)
units = "Reads"
non_zero = which(gene_data[,"Read_Count"] > 0)
data = gene_data[non_zero,"Read_Count"]
divider = quantile(data, probs=0.95)
title = paste(library_name, " - Distribution of Gene Read Counts", sep="")
x_label = paste("Read Count (", units, ")", sep="")
y_label = paste("'Expressed' Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneReads.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, seq=""),
                paste("Max = ", max, " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red", "orange", "black"))
dev.off()

lower = which(gene_data[,"Read_Count"] > 0 & gene_data[,"Read_Count"] <= divider)
data = gene_data[lower,"Read_Count"]
x_label = paste("Read Count (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("'Expressed' Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneReads_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", round(min, digits=4), " ", units, sep=""),
                paste("Median = ", round(median, digits=1), " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, seq=""),
                paste("Max = ", round(max, digits=1), " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red", "orange", "black"))
dev.off()


#I.) PERCENT GENE COVERAGE 1X, 5X, 10X, 50X, 100X, 500X
units = "% Coverage"
non_zero = which(gene_data[,"Percent_Coverage_1x"] > 0)
data1 = gene_data[non_zero,"Percent_Coverage_1x"]
data2 = gene_data[non_zero,"Percent_Coverage_5x"]
data3 = gene_data[non_zero,"Percent_Coverage_10x"]
data4 = gene_data[non_zero,"Percent_Coverage_50x"]
data5 = gene_data[non_zero,"Percent_Coverage_100x"]
data6 = gene_data[non_zero,"Percent_Coverage_500x"]
names = c(">= 1X", ">= 5X", ">= 10X", ">= 50X", ">= 100X", ">= 500X")
title = paste(library_name, " - Distribution of percent gene coverage", sep="")
x_label = "X coverage cutoff"
y_label = paste("Percent Gene Coverage (n = ", prettyNum(length(data1), big.mark=",")," genes)", sep="")
tiff(file="PercentGeneCoverage_1x-500x.tiff", width=700, height=700, compression="none")
cols=rainbow(9)
par(bg=bg, font.main = 2, font.lab = 2); boxplot(x = list(data1,data2,data3,data4,data5,data6), names=names, main=title, xlab=x_label, ylab=y_label, col=cols[1:6], col.main = "black", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
dev.off();


#J.) AVERAGE GENE COVERAGE RAW - 'detected' genes only (coverage > 0)
summary(gene_data[,"Average_Coverage_RAW"])  #Dynamic range of gene expression levels is 0 to 10^5
units = "X"
non_zero = which(gene_data[,"Average_Coverage_RAW"] > 0)
data = gene_data[non_zero,"Average_Coverage_RAW"]
divider = quantile(data, probs=0.95)
title = paste(library_name, " - Distribution of Average Gene Base Coverage", sep="")
x_label = paste("Mean Base Coverage (", units, ")", sep="")
y_label = paste("'Expressed' Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneCoverage.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", round(min, digits=4), " ", units, sep=""),
                paste("Median = ", round(median, digits=1), " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", round(max, digits=1), " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()

lower = which(gene_data[,"Average_Coverage_RAW"] > 0 & gene_data[,"Average_Coverage_RAW"] <= divider)
data = gene_data[lower,"Average_Coverage_RAW"]
x_label = paste("Mean Base Coverage (<= ", round(divider, digits=1), " ", units, "; 95th percentile)", sep="")
y_label = paste("'Expressed' Gene Count  (n = ", length(data), ")", sep="")
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
tiff(file="GeneCoverage_95percentile.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); hist(x=data, main=title, xlab=x_label, ylab=y_label, col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=100)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", round(min, digits=4), " ", units, sep=""),
                paste("Median = ", round(median, digits=1), " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", round(max, digits=1), " ", units, sep=""))
legend("topright", legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off()


##########################################################################################################
#Now try some graphs where only the Expressed features are summarized...


#Limit analysis to only 'Expressed' genes and exons ('Expressed==1')
expressed_genes=which(gene_data[,"Expressed"]==1)

units = "% Coverage"
data1 = gene_data[expressed_genes,"Percent_Coverage_1x"]
data2 = gene_data[expressed_genes,"Percent_Coverage_5x"]
data3 = gene_data[expressed_genes,"Percent_Coverage_10x"]
data4 = gene_data[expressed_genes,"Percent_Coverage_50x"]
data5 = gene_data[expressed_genes,"Percent_Coverage_100x"]
data6 = gene_data[expressed_genes,"Percent_Coverage_500x"]
names = c(">= 1X", ">= 5X", ">= 10X", ">= 50X", ">= 100X", ">= 500X")
title = paste(library_name, " - Distribution of percent gene coverage", sep="")
x_label = "X coverage cutoff"
y_label = paste("Percent Gene Coverage (n = ", prettyNum(length(data1), big.mark=",")," genes)", sep="")
tiff(file="PercentGeneCoverage_1x-500x_Expressed.tiff", width=700, height=700, compression="none")
cols=rainbow(9)
par(bg=bg, font.main = 2, font.lab = 2); boxplot(x = list(data1,data2,data3,data4,data5,data6), names=names, main=title, xlab=x_label, ylab=y_label, col=cols[1:6], col.main = "black", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
dev.off();


#Print out some summary stats for percent gene coverage
#Median % coverage at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X
#Number of genes sequenced over 50% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X
#Number of genes sequenced over 75% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X
#Number of genes sequenced over 90% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X
stats_file = "GeneCoverageStats_ExpressedGenes.txt"
string = "\nMedian % coverage at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X"; write(string, file=stats_file);
string = paste ("Median at 1X = ", round(median(data1),2), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("Median at 5X = ", round(median(data2),2), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("Median at 10X = ", round(median(data3),2), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("Median at 50X = ", round(median(data4),2), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("Median at 100X = ", round(median(data5),2), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("Median at 500X = ", round(median(data6),2), sep=""); write(string, file=stats_file, append=TRUE);

string = paste ("\nTotal # genes detected by at least one read = ", length(non_zero), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("\nTotal # genes detected as expressed = ", length(expressed_genes), sep=""); write(string, file=stats_file, append=TRUE);

string = "\n# of genes sequenced of 50% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X"; write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 1X = ", length(which(data1 > 50)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 5X = ", length(which(data2 > 50)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 10X = ", length(which(data3 > 50)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 50X = ", length(which(data4 > 50)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 100X = ", length(which(data5 > 50)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 500X = ", length(which(data6 > 50)), sep=""); write(string, file=stats_file, append=TRUE);

string = "\n# of genes sequenced of 75% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X"; write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 1X = ", length(which(data1 > 75)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 5X = ", length(which(data2 > 75)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 10X = ", length(which(data3 > 75)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 50X = ", length(which(data4 > 75)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 100X = ", length(which(data5 > 75)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 500X = ", length(which(data6 > 75)), sep=""); write(string, file=stats_file, append=TRUE);

string = "\n# of genes sequenced of 90% of their length at X cutoffs: 1X, 5X, 10X, 50X, 100X, 500X"; write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 1X = ", length(which(data1 > 90)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 5X = ", length(which(data2 > 90)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 10X = ", length(which(data3 > 90)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 50X = ", length(which(data4 > 90)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 100X = ", length(which(data5 > 90)), sep=""); write(string, file=stats_file, append=TRUE);
string = paste ("# genes at 500X = ", length(which(data6 > 90)), sep=""); write(string, file=stats_file, append=TRUE);




