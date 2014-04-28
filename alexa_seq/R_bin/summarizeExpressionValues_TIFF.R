#!/usr/bin/env Rscript
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

library(geneplotter)
library(RColorBrewer)
library(gcrma)
library(Cairo)

#Options:
#[1] datafile1
#[2] datafile2
#[3] datafile3
#[4] results_dir
#[5] library_name

#Example usage:
#/home/malachig/svn/solexa_analysis/R_bin/summarizeExpressionValues.R /projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_joined.txt /projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_RAW_joined.txt /projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_ER.txt /projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_SIR.txt /projects/malachig/solexa/figures_and_stats/HS04391/Expression/ HS04391

#NOTES:
#The input file must contain: a unique ID (gene, exon, junction), EnsEMBL Gene ID, Gene Name, and two columns of expression data (library A and B) to be compared
#These expression values should be raw expression values (normalization to the smaller library will be done by this script)
#The output file will contain these two data columns, normalized data values, the log2 difference, fold change, pvalue, and corrected pvalues (qvalue)
args = (commandArgs(TRUE))
datafile1 = args[1]
datafile2 = args[2]
datafile3 = args[3]
datafile4 = args[4]
results_dir = args[5]
library_name = args[6]

#datafile1 = "/projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_joined.txt"
#datafile2 = "/projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_RAW_joined.txt"
#datafile3 = "/projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_ER.txt"
#datafile4 = "/projects/malachig/solexa/figures_and_stats/HS04391/temp/Average_Coverage_NORM1_SIR.txt"
#results_dir = "/projects/malachig/solexa/figures_and_stats/HS04391/Expression/"
#library_name = "HS04391"

data=read.table(datafile1, header = TRUE, na.strings = "NA", sep=",")
data2=read.table(datafile2, header = TRUE, na.strings = "NA", sep=",")
data_er=read.table(datafile3, header = TRUE, na.strings = "NA", sep="\t")
data_sr=read.table(datafile4, header = TRUE, na.strings = "NA", sep="\t")
setwd(results_dir)

#Define an output stats file
stats_file = paste(library_name, "_stats.txt", sep = "")

#1.) Write out a simple summary of normalized coverage values for each sequence type
string = paste ("Summary statistics for library: ", library_name, sep="")
write(string, file=stats_file)

for (i in 1:length(data[1,])){
  x = formatC(summary(data[,i]), format="fg")
  string = paste("\nCoverge values for: ", names(data[i]), " (Library ", library_name, ")", sep="")
  write(string, file=stats_file, append=TRUE)
  write.table(x[], file = stats_file, append=TRUE, col.names=FALSE, quote=FALSE)
}

#(Gene Transcript ExonRegion Junction KnownJunction NovelJunction Boundary KnownBoundary NovelBoundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion)
cols = rainbow(7)
col_groups1 = c(cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[6], cols[6], cols[7], cols[7], cols[7])
col_groups2 = c(cols[1], cols[2], cols[3], cols[4], cols[4], cols[4], cols[5], cols[5], cols[5], cols[6], cols[6], cols[6], cols[7], cols[7], cols[7])


#2.) Determine the percentiles for each seq type (normal and log2 scale)
#Plot these on a log2 scale and note the expression value that corresponds to the 95th percentile for each data type
percentiles = data.frame(seq(0, 1, 0.01))
for (i in 1:length(data[1,])){
  percentiles[,i+1] = quantile(data[,i], probs=seq(0, 1, 0.01), na.rm=TRUE)
}
names(percentiles) = c("percentile", names(data))

percentiles_log2 = data.frame(seq(0, 1, 0.01))
for (i in 1:length(data[1,])){
  percentiles_log2[,i+1] = quantile(log2(data[,i]), probs=seq(0, 1, 0.01), na.rm=TRUE)
}
names(percentiles_log2) = c("percentile", names(data))

#bizarre bug in R
#If value is '0.95' must be in quotes.  For all other values must not be in quotes...
target = "0.95"
pos = which(percentiles_log2[,"percentile"] == target)

tiff_percentiles = "PercentilesPlot_Log2.tiff"
CairoTIFF(filename=tiff_percentiles, width=800, height=800, pointsize=12, bg="white")
par(font.main = 2, font.lab = 2)
plot (x=percentiles_log2[,"percentile"], y=percentiles_log2[,"ExonRegion"], type="l", lty=1, lwd=2, col="dark green", 
      main="Percentiles of Expression values for Exon Regions and Intronic/Intergenic regions",
      xlab="Percentile", ylab="Expression value (Log2 average coverage)", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5)
lines (x=percentiles_log2[,"percentile"], y=percentiles_log2[,"SilentIntronRegion"], pch=20, col="orange", lwd=2)
lines (x=percentiles_log2[,"percentile"], y=percentiles_log2[,"SilentIntergenicRegion"], pch=20, col="red", lwd=2)

lines(x=rep(target,pos), y=percentiles_log2[1:pos,"SilentIntronRegion"], col="black", lty=2, lwd=2)
silent_intron_cutoff = percentiles_log2[pos,"SilentIntronRegion"]
silent_intergenic_cutoff = percentiles_log2[pos,"SilentIntergenicRegion"]
abline (h=silent_intron_cutoff, col="orange", lty=2, lwd=2)
abline (h=silent_intergenic_cutoff, col="red", lty=2, lwd=2)

legend_location = "topleft"
legend_text = c(paste("Exon Region (ER)", sep=""),
                paste("'Silent' Intronic (95th percetile = ", round(silent_intron_cutoff, digits=2), ")", sep=""),
                paste("'Silent' Intergenic (95th percentile = ", round(silent_intergenic_cutoff, digits=2), ")", sep=""),
                paste("95th percentile", sep = ""))
legend(legend_location, legend=legend_text, lty=c(1,1,1,2), lwd=c(2,2,2,2), col=c("dark green","orange","red", "black"), cex=1.5)
dev.off()
2^silent_intron_cutoff
2^silent_intergenic_cutoff

#Write the chosen 95th percentile cutoff to the stats output file 
string1 = paste ("\n\n95th percentile of silent intron regions for library (", library_name, ") ", "is: ", round(2^silent_intron_cutoff, digits=2), " (log2 = ", round(silent_intron_cutoff, digits=2), ")", sep="")
string2 = paste ("\n\n95th percentile of silent intergenic regions for library (", library_name, ") ", "is: ", round(2^silent_intergenic_cutoff, digits=2), " (log2 = ", round(silent_intergenic_cutoff, digits=2), ")", sep="")
write(string1, file=stats_file, append=TRUE)
write(string2, file=stats_file, append=TRUE)


#3-A.) Generate box-plots of normalized expression data for each data type - log2 scale

#Limit display to some types only
display_types = c("Gene","Transcript","ExonRegion","Junction","Boundary", "Intron", "ActiveIntronRegion","SilentIntronRegion","Intergenic","ActiveIntergenicRegion","SilentIntergenicRegion");
display_names = c("Gene","Transcript","Exon Region","Junction","Boundary", "Intron", "Active Intron","Silent Intron","Intergenic","Active Intergenic","Silent Intergenic");
temp_data=log2(data[display_types])
names(temp_data)=display_names

tiff_boxplot_all_types = "ExpressionValues_NORM1_Log2_AllTypes_bplot.tiff";
CairoTIFF(filename=tiff_boxplot_all_types, width=1000, height=1000, pointsize=12, bg="white")
par(font.main = 2, font.lab = 2, mar=(c(12, 4, 4, 2)+0.1))
boxplot(temp_data, main="Expression by Sequence Type", col=col_groups1, ylab="Normalized Log2 Expression (Average Coverage)", 
        col.lab = gray(.1), cex.main = 1.7, cex.lab = 1.6, cex.axis=1.5, las=2)
abline (h=silent_intergenic_cutoff, col="black", lty=2, lwd=2)
legend_location = "topright"
legend_text = "95th Percentile of Silent Intergenic Regions"
legend(legend_location, legend=legend_text, lty=c(2), lwd=c(2), col=c("black"), cex=1.5)
dev.off();

#3-B.) Now do the same thing but for RAW expression data - log2 scale
tiff_boxplot_all_types = "ExpressionValues_RAW_Log2_AllTypes_bplot.tiff";
CairoTIFF(filename=tiff_boxplot_all_types, width=1000, height=1000, pointsize=12, bg="white")
par(font.main = 2, font.lab = 2, mar=(c(12, 4, 4, 2)+0.1))
temp_data=log2(data2[display_types])
names(temp_data)=display_names
boxplot(temp_data, main="Expression by Sequence Type", col=col_groups1, ylab="RAW Log2 Expression (Average Coverage)", 
        col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5, las=2)
dev.off();


#4.) Create a histogram of the distribution of expression values (log2) for all data types - all on one plot
tiff_hist_all_types = "ExpressionValues_Log2_AllTypes_hist.tiff";
CairoTIFF(filename=tiff_hist_all_types, width=1920, height=1200, pointsize=12, bg="white")

par(mfrow=c(4,4))
par(font.main = 2, font.lab = 2)
for (i in 1:length(data[1,])){
  hist(log2(data[,i]), breaks=100, main=names(data[i]), col=col_groups2[i], xlab="Log2 Expression", ylab="Frequency", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
  abline (v=silent_intergenic_cutoff, col="black", lty=2, lwd=2)
}
dev.off();

#5.) Create a histogram of the distribution of expression values (log2) for all data types - one plot for each
for (i in 1:length(data[1,])){
  filename = paste("ExpressionValue_Log2_hist_", names(data[i]), ".tiff", sep="")
  CairoTIFF(filename=filename, width=700, height=700, pointsize=12, bg="white")
  par(font.main = 2, font.lab = 2)
  hist(log2(data[,i]), breaks=100, main=names(data[i]), col=col_groups2[i], xlab="Normalized Log2 Expression (Average Coverage)", ylab="Frequency", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
  abline (v=silent_intergenic_cutoff, col="black", lty=2, lwd=2)
  legend_location = "topleft"
  legend_text = "95th Percentile of Silent Intergenic Regions"
  legend(legend_location, legend=legend_text, lty=c(2), lwd=2, col=c("black"))
  dev.off();
}

#6.) Create scatter plots comparing the expression values of genes compared to ExonRegions or SilentIntronRegions within them

#Exon regions first - expect a very good correlation, since the gene expression estimate and exon expression estimates are derived from the same mappings to known transcripts
filename = paste("ExonRegionVsGeneExpression_Log2.tiff", sep="")
CairoTIFF(filename=filename, width=700, height=700, pointsize=12, bg="white")

x=log2(data_er[,"Gene_Expression"])
y=log2(data_er[,"ER_Expression"])

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
par(font.main = 2, font.lab = 2)
smoothScatter(x, y, xlab="Gene Expression (log2)", ylab="Exon Region Expression (log2)", main="Correlation between gene and exon expression estimates", colramp=colors, nbin=275,
              col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(loess.smooth(x, y, span = 2/3, degree = 1, family = "gaussian", evaluation = 100), col="black", lty=2, lwd=2)
legend_location = "topleft"
corr = cor(x,y, method="spearman") 
legend_text = c(paste("Correlation = ", round(corr, digits=2), " (Spearman)", sep=""))
legend(legend_location, legend=legend_text, lty=c(2), lwd=2, col=c("black"), cex=1.3)

dev.off()

#Silent Intron regions now - expect a correlation if intron signal primarily comes from pre-mRNA contamination and this is related to expression level of the gene
filename = paste("SilentIntronicVsGeneExpression_Log2.tiff", sep="")
CairoTIFF(filename=filename, width=700, height=700, pointsize=12, bg="white")

x=log2(data_sr[,"Gene_Expression"])
y=log2(data_sr[,"SR_Expression"])

colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
par(font.main = 2, font.lab = 2)
smoothScatter(x, y, xlab="Gene Expression (log2)", ylab="'Silent' Intron Region Expression (log2)", main="Correlation between gene and intron expression estimates", colramp=colors, nbin=275,
              col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

#Get the coefficients of a simple line fit to the silent intron expression values.
#These can be used to estimate the intronic sequence contamination of any gene based on the gene expression level
#This would allow us to decide whether introns or exon-boundaries are expressed above background for a particular gene and the cutoff will depend on the expression level of that gene
fit = lm(y~x);
fitted.y = fit$fitted.values
lines(x, fitted.y, col="black", lty=1, lwd=2)

#Now get the coefficients of a simple line to the 95th percentiles of sliding windows of the silent intron expression values
join = data.frame(x,y)
join_sort = join[order(join[,"x"]),]
divider = 1000
d = floor(length(x)/divider)
lower = seq(1, length(x), d) 
upper = c(seq((1+d)-1, length(x), d), length(x))
if (length(upper) == length(lower)+1){
  upper = seq((1+d)-1, length(x), d)
}

cutoff_p = data.frame(lower,upper)
upper_p = as.numeric(target)
lower_p = 1-as.numeric(target)
for (i in 1:length(lower)){
  y_upper = quantile(join_sort[lower[i]:upper[i],"y"], probs=upper_p)
  y_lower = quantile(join_sort[lower[i]:upper[i],"y"], probs=lower_p)
  xmed = median(join_sort[lower[i]:upper[i],"x"])
  cutoff_p[i,"xmed"] = xmed
  cutoff_p[i,"y_upper"] = y_upper
  cutoff_p[i,"y_lower"] = y_lower
}
fit_upper = lm(cutoff_p[,"y_upper"]~cutoff_p[,"xmed"])
fit_lower = lm(cutoff_p[,"y_lower"]~cutoff_p[,"xmed"])
lines(cutoff_p[,"xmed"], fit_upper$fitted.values, col="black", lty=2, lwd=2)
lines(cutoff_p[,"xmed"], fit_lower$fitted.values, col="black", lty=2, lwd=2)

legend_location = "topleft"
corr = cor(x,y, method="spearman") 
legend_text = c(paste("Linear Fit. Correlation = ", round(corr, digits=2), " (Spearman)", sep=""),
                "5th and 95th percentile bounds")
legend(legend_location, legend=legend_text, lty=c(1,2), lwd=2, col=c("black","black"), cex=1.4)
dev.off() 

#Print out the coefficients of the upper line fit - these will be used to estimating the gene-by-gene expression cutoff level
string1 = paste ("\n\nCoefficients of Line fit for sliding window of 95th percentiles of silent intron region expression", sep="")
string2 = paste ("\nSlope (log2) = ", fit_upper$coef[2], sep = "")
string3 = paste ("\nIntercept (log2) = ", fit_upper$coef[1], sep = "")
write(string1, file=stats_file, append=TRUE)
write(string2, file=stats_file, append=TRUE)
write(string3, file=stats_file, append=TRUE)

