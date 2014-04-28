#!/usr/bin/env Rscript
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Example usage:
#Options: 
#[1] input_data
#[2] data_length
#[3] data_name
#[4] units
#[5] output_dir
#[6] breaks

bg = "white"

#/home/malachig/svn/solexa_analysis/R_bin/hist.R /projects/malachig/solexa/figures_and_stats/HS04391/temp/temp_data.txt 3637354 "Fragment Size" bp /projects/malachig/solexa/figures_and_stats/HS04391/ENST_v49/ 100


#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1];
record_count = as.numeric(args[2]);
data_name = args[3];
units = args[4];
output_dir = args[5];
breaks = as.numeric(args[6]);

print ("User specified the following options:")
print (infile)
print (record_count)
print (data_name)
print (units)
print (output_dir)
print (breaks)

#Import the data
#data = scan(file=infile, what=integer(0), nmax=record_count)
data = scan(file=infile, what=numeric(0), nmax=record_count)

#Relative skewness > 0 indicates positive skewness (a longer right tail) and relative skewness < 0 indicates negative skewness (a longer left tail).
library(e1071)
skew = skewness(data)

filename=paste(output_dir, gsub(" +", "", data_name), "_Hist.tiff", sep="")
filename_lower = "lower_undef.txt"
filename_upper = "upper_undef.txt"
legend_location = "topright"

#If the data skew is >= 0 (a long right tail) then we want to display: (<= 95th percentile; bulk of data is in 'lower') and (> 95th percentile; minority of data is in 'upper')
data_split = 0;
percentile = "percentile_undef";

if (skew >= 0){
  data_split = quantile(data, probs=0.95)
  filename_lower=paste(output_dir, gsub(" +", "", data_name), "_below_95th_percentile_Hist.tiff", sep="")
  filename_upper=paste(output_dir, gsub(" +", "", data_name), "_above_95th_percentile_Hist.tiff", sep="")
  percentile = "95th"
  legend_location = "topright"
}

#If the data skew is < 0 (a long left tail) then we want to display: (<= 5th percentile; minority of data is in 'lower') and (>5th percentile; bulk of data is in 'upper')
if (skew < 0){
  data_split = quantile(data, probs=0.05)
  filename_lower=paste(output_dir, gsub(" +", "", data_name), "_below_05th_percentile_Hist.tiff", sep="")
  filename_upper=paste(output_dir, gsub(" +", "", data_name), "_above_05th_percentile_Hist.tiff", sep="")
  percentile = "5th"
  legend_location = "topleft"
}

lower = which(data <= data_split)
upper = which(data > data_split)


#A.) All data without dividing into two pieces
string = paste ("DATA.  length:", length(data), " min:", min(data), " max:", max(data))
print(string)

main_title = paste("Distribution of ", data_name, " values ", "(", units, ")", sep="");
y_label = paste("Frequency (n = ", prettyNum(length(data), big.mark=","), ")", sep="")
x_label = paste(data_name, " (", units, ")", sep="")

tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data, col="blue", main=main_title, xlab=x_label, ylab=y_label, breaks=breaks, cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
min = min(data)
median = median(data)
mean = mean(data)
max = max(data)
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"), cex=1.4)
dev.off();


#B.) Lower range of data - Up to 95th percentile
skew = skewness(data[lower])
if (skew >= 0){legend_location = "topright"}
if (skew <= 0){legend_location = "topleft"} 

string = paste ("DATA LOWER.  length:", length(data[lower]), " min:", min(data[lower]), " max:", max(data[lower]))
print(string)

main_title = paste("Distribution of ", data_name, " values ", "(", units, ")", sep="");
y_label = paste("Frequency (n = ", prettyNum(length(data[lower]), big.mark=","), ")", sep="")
x_label = paste(data_name, " (<= ", round(data_split, digits=1), " ", units, "; ", percentile," percentile)", sep="")

tiff(file=filename_lower, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data[lower], col="blue", main=main_title, xlab=x_label, ylab=y_label, breaks=breaks, cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
min = min(data[lower])
median = median(data[lower])
mean = mean(data[lower])
max = max(data[lower])
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"), cex=1.4)
dev.off();

#C.) Upper range of data - Above 95th percentile
skew = skewness(data[upper])
if (skew >= 0){legend_location = "topright"}
if (skew <= 0){legend_location = "topleft"} 

string = paste ("DATA UPPER.  length:", length(data[upper]), " min:", min(data[upper]), " max:", max(data[upper]))
print(string)

main_title = paste("Distribution of ", data_name, " values ", "(", units, ")", sep="");
y_label = paste("Frequency (n = ", prettyNum(length(data[upper]), big.mark=","), ")", sep="")
x_label = paste(data_name, " (> ", round(data_split, digits=1), " ", units, "; ", percentile," percentile)", sep="")


tiff(file=filename_upper, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data[upper], col="blue", main=main_title, xlab=x_label, ylab=y_label, breaks=breaks, cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
min = min(data[upper])
median = median(data[upper])
mean = mean(data[upper])
max = max(data[upper])
abline(v=median, col="red", lty=2, lwd=2)
abline(v=mean, col="orange", lty=2, lwd=2)
legend_text = c(paste("Min. = ", min, " ", units, sep=""),
                paste("Median = ", median, " ", units, sep=""),
                paste("Mean = ", round(mean, digits=1), " ", units, sep=""),
                paste("Max = ", max, " ", units, sep=""))
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"), cex=1.4)
dev.off();








