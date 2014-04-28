#!/usr/bin/env Rscript
#Written by Malachi Griffith 
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Generate basic plots describing exon-junction data
#A.) Number of known exon-junctions involving no exon skip (S0), 1 exon skip (S1), 2 exons skipped (S2), etc.
#B.) Number of expressed junctions that are S0, S1, S2, etc.
#C.) Plot both A and B on the same plot as percents (i.e. percent of all known events that are S0, S1, ... AND the percent of expressed events that are S0, S1, ...

#D.) Plot percent of all S0, S1, S2, etc. events that are known to occur
#    - Also plot percent of all S0, S1, S2, etc. events that are expressed

#E.) Expression-level (as histogram) for all exon-junctions by exon skip #(S0, S1, S2, S3, etc.)

#Options:
#[1] datafile1
#[2] results_dir

#Example usage:

#/home/malachig/svn/solexa_analysis/R_bin/summarizeExonJunctionData.R  /projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_JunctionExpression_v53.txt  /projects/malachig/solexa/figures_and_stats/HS04391/Junctions_v53/
args = (commandArgs(TRUE))
datafile = args[1]
results_dir = args[2]

print ("ARGS supplied by user")
print (datafile)
print (results_dir)

#datafile = "/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_JunctionExpression_v53_SELECTED_COLUMNS.txt"
#results_dir = "/projects/malachig/solexa/figures_and_stats/HS04391/Junctions_v53/"

junctions = read.table(datafile, header = TRUE, na.strings = "na", sep="\t")
setwd(results_dir)

bg="white"

#A.) Number of known exon-junctions involving no exon skip (S0), 1 exon skip (S1), 2 exons skipped (S2), etc.

#A-1.) All known junctions (S0-Sn)
known_junctions = junctions[which(junctions[,"Supporting_EnsEMBL_Count"] >= 1),]
novel_junctions = junctions[which(junctions[,"Supporting_EnsEMBL_Count"] == 0),]

data = known_junctions[,"Exons_Skipped"]
length(data)

units = "Exons skipped"
y_label = paste("Frequency (n = ", length(data), ")", sep="")
x_label = units

tiff(file="Exons_Skipped_KnownJunctions_S0-Sn.tiff", width=700, height=700, compression="none")
main_title = "Distribution of Exons Skipped (known exon junctions; S0-Sn)"
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data, col="blue", main=main_title, xlab=x_label, ylab=y_label, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=50)
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
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off();

#A-2.) All known skipping junctions (S1-Sn)
tiff(file="Exons_Skipped_KnownJunctions_S1-Sn.tiff", width=700, height=700, compression="none")
known_junctions_S1_Sn = known_junctions[which(known_junctions[,"Exons_Skipped"] > 0),]
data = known_junctions_S1_Sn[,"Exons_Skipped"] 
length(data)

main_title = "Distribution of Exons Skipped (known exon junctions; S1-Sn)"
y_label = paste("Frequency (n = ", length(data), ")", sep="")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data, col="blue", main=main_title, xlab=x_label, ylab=y_label, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=50)
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
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off();


#B.) Number of expressed junctions that are S0, S1, S2, etc. (expressed junctions can be known or otherwise)
coverage_cutoff = 75;
expressed_i = which(junctions[,"Expressed"] == 1 & junctions[,"Percent_Coverage_1x"] > coverage_cutoff)
expressed_junctions = junctions[expressed_i,]

expressed_novel_i = which(junctions[,"Expressed"] == 1 & junctions[,"Percent_Coverage_1x"] > coverage_cutoff & junctions[,"Supporting_EnsEMBL_Count"] == 0)
expressed_novel_junctions = junctions[expressed_novel_i,]

length(expressed_i) #Expressed
length(which(expressed_junctions[,"Supporting_EnsEMBL_Count"] >=1)) #Expressed and known
length(expressed_i) - length(which(expressed_junctions[,"Supporting_EnsEMBL_Count"] >=1)) #Expressed and 'novel'
#Expressed and 'novel' but have mRNA or EST support
length(which(expressed_junctions[,"Supporting_EnsEMBL_Count"] ==0 & (expressed_junctions[,"Supporting_mRNA_Count"] >= 1 | expressed_junctions[,"Supporting_EST_Count"] >= 1)))


#B-1.) All expressed junctions (S0-Sn)
tiff(file="Exons_Skipped_ExpressedJunctions_S0-Sn.tiff", width=700, height=700, compression="none")
data=expressed_junctions[,"Exons_Skipped"]
length(data)
main_title = "Distribution of Exons Skipped (expressed exon junctions; S0-Sn)"
y_label = paste("Frequency (n = ", length(data), ")", sep="")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data, col="blue", main=main_title, xlab=x_label, ylab=y_label, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=50)
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
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off();

#B-2.) All expressed skipping junctions (S1-Sn)
expressed_junctions_S1_Sn = expressed_junctions[which(expressed_junctions[,"Exons_Skipped"] > 0),]
data = expressed_junctions_S1_Sn[,"Exons_Skipped"] 
length(data)

tiff(file="Exons_Skipped_ExpressedJunctions_S1-Sn.tiff", width=700, height=700, compression="none")
data=expressed_junctions_S1_Sn[,"Exons_Skipped"]
length(data)
main_title = "Distribution of Exons Skipped (expressed exon junctions; S1-Sn)"
y_label = paste("Frequency (n = ", length(data), ")", sep="")
par(bg=bg, font.main = 2, font.lab = 2)
hist(x = data, col="blue", main=main_title, xlab=x_label, ylab=y_label, col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks=50)
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
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(0,2,2,0), lwd=2, col=c("black","red","orange","black"))
dev.off();


#How many of the expressed exon skipping events are known?
length(expressed_junctions_S1_Sn[,"Supporting_EnsEMBL_Count"]) #Expressed skipping junctions
length(which(expressed_junctions_S1_Sn[,"Supporting_EnsEMBL_Count"] >=1)) #Expressed skipping junctions that are known
length(which(expressed_junctions_S1_Sn[,"Supporting_EnsEMBL_Count"] ==0)) #Expressed skipping junctions that are novel
length(which(expressed_junctions_S1_Sn[,"Supporting_EnsEMBL_Count"] ==0 & (expressed_junctions_S1_Sn[,"Supporting_mRNA_Count"] >= 1 | expressed_junctions_S1_Sn[,"Supporting_EST_Count"] >= 1)))


#C.) Plot both A and B on the same plot as percents (i.e. percent of all known events that are S0, S1, ... AND the percent of expressed events that are S0, S1, ...
known_skip_counts = table(known_junctions[,"Exons_Skipped"])
known_skip_sum = sum(known_skip_counts[2:length(known_skip_counts)])
known_skip_p = (known_skip_counts/known_skip_sum)*100
expressed_skip_counts = table(expressed_junctions[,"Exons_Skipped"])
expressed_skip_sum = sum(expressed_skip_counts[2:length(expressed_skip_counts)])
expressed_skip_p = (expressed_skip_counts/expressed_skip_sum)*100
expressed_novel_skip_counts = table(expressed_novel_junctions[,"Exons_Skipped"])
expressed_novel_skip_sum = sum(expressed_novel_skip_counts[2:length(expressed_novel_skip_counts)])
expressed_novel_skip_p = (expressed_novel_skip_counts/expressed_novel_skip_sum)*100

#Report the N values for each category to be summarized
print ("\nSummary of N values")
string = paste("Known junctions = ", known_skip_sum, sep="")
print (string)
string = paste("Expressed junctions = ", expressed_skip_sum, sep="")
print (string)
string = paste("Novel Expressed junctions = ", expressed_novel_skip_sum, sep="")
print (string)

#Produce a random sampling of exon-exon junctions - This should be the same size as all 'known' junctions
i = (1:length(junctions[,1]))
i_rand = sample(i, length(known_junctions[,1]), replace=FALSE)
random_junctions=junctions[i_rand,]
random_skip_counts = table(random_junctions[,"Exons_Skipped"])
random_skip_sum = sum(random_skip_counts[2:length(random_skip_counts)])
random_skip_p = (random_skip_counts/random_skip_sum)*100

main_title = "Exon-skipping of all known, all expressed, and novel expressed junctions"
tiff(file="Exons_Skipped_KnownVsExpressedJunctions_S1-S10.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=names(known_skip_p[2:11]), y=known_skip_p[2:11], ylim=c(0,80), type="l", lwd=3, col="dark green", lty=2, xlab="Number of exons skipped", ylab="Percentage of exon-skipping junctions", main=main_title, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(x=names(expressed_skip_p[2:11]), y=expressed_skip_p[2:11], lwd=3, col="blue", lty=1)
lines(x=names(expressed_novel_skip_p[2:10]), y=expressed_novel_skip_p[2:10], lwd=3, col="magenta", lty=3)
lines(x=names(random_skip_p[2:11]), y=random_skip_p[2:11], lwd=3, col="red", lty=4)
abline(h=0, col=gray(0.5), lty=2, lwd=1)

legend_text = c("All Known Junctions","All Expressed Junctions", "Novel Expressed Junctions", "Random")
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(2,1,3,4), lwd=3, col=c("dark green","blue", "magenta", "red"), cex=1.4)
dev.off()


#D.)  Plot percent of all S0, S1, S2, etc. events that are known to occur
#    - Also plot percent of all S0, S1, S2, etc. events that are expressed
skip_counts = table(junctions[,"Exons_Skipped"])               #Number of junctions that are S0-Sn, known or not (i.e. all theoretical)
skip_counts_known = table(known_junctions[,"Exons_Skipped"])   #Number of known junctions that are S0-Sn

#Known percents - i.e. percent of all junctions of each skip number that are known
s0_kp = (skip_counts_known["0"]/skip_counts["0"])*100
s1_kp = (skip_counts_known["1"]/skip_counts["1"])*100
s2_kp = (skip_counts_known["2"]/skip_counts["2"])*100
s3_kp = (skip_counts_known["3"]/skip_counts["3"])*100
s4_kp = (skip_counts_known["4"]/skip_counts["4"])*100
s5_kp = (skip_counts_known["5"]/skip_counts["5"])*100
s6_sN_kp = (sum(skip_counts_known[7:length(skip_counts_known)]) / sum(skip_counts[7:length(skip_counts)]))*100 

kp = c(s0_kp,s1_kp,s2_kp,s3_kp,s4_kp,s5_kp,s6_sN_kp)
names(kp) = c("S0", "S1", "S2", "S3", "S4", "S5", "S6+")

#Observed/Expressed percents
skip_counts_expressed = table(expressed_junctions[,"Exons_Skipped"])
s0_ep = (skip_counts_expressed["0"]/skip_counts["0"])*100
s1_ep = (skip_counts_expressed["1"]/skip_counts["1"])*100
s2_ep = (skip_counts_expressed["2"]/skip_counts["2"])*100
s3_ep = (skip_counts_expressed["3"]/skip_counts["3"])*100
s4_ep = (skip_counts_expressed["4"]/skip_counts["4"])*100
s5_ep = (skip_counts_expressed["5"]/skip_counts["5"])*100
s6_sN_ep = (sum(skip_counts_expressed[7:length(skip_counts_expressed)]) / sum(skip_counts[7:length(skip_counts)]))*100 
ep = c(s0_ep,s1_ep,s2_ep,s3_ep,s4_ep,s5_ep,s6_sN_ep)
names(ep) = c("S0", "S1", "S2", "S3", "S4", "S5", "S6+")

tiff(file="Exons_Skipped_PercentOfAllJunctions_KnownVsExpressed_S0-Sn.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(y=kp, x=(0:6), type="l", lwd=2, col="blue", lty=1, xlab="Number of exons skipped", ylab="Percentage of ALL junctions", main=main_title, 
     col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
lines(y=ep, x=(0:6), lwd=2, col="red", lty=2)
abline(h=0, col=gray(0.5), lty=2, lwd=1)
legend_text = c("Known Junctions","Expressed Junctions")
legend_location = "topright"
legend(legend_location, legend=legend_text, lty=c(1,2), lwd=2, col=c("blue","red"))
dev.off()

#E.) Expression-level (as box plot) for all exon-junctions by exon skip #(S0, S1, S2, S3, etc.)

#E-1.) Expressed junctions by exon skip #(S0, S1, S2, S3, etc.)  - First limit to only those junctions actually expressed
coverage_cutoff = 70;
expressed_i = which(junctions[,"Expressed"] == 1 & junctions[,"Percent_Coverage_1x"] > coverage_cutoff)
x = junctions[expressed_i,]

s0 = log2(x[which(x[,"Exons_Skipped"] == 0),"Average_Coverage_NORM1"])
s1 = log2(x[which(x[,"Exons_Skipped"] == 1),"Average_Coverage_NORM1"])
s2 = log2(x[which(x[,"Exons_Skipped"] == 2),"Average_Coverage_NORM1"])
s3 = log2(x[which(x[,"Exons_Skipped"] == 3),"Average_Coverage_NORM1"])
s4 = log2(x[which(x[,"Exons_Skipped"] == 4),"Average_Coverage_NORM1"])
s5 = log2(x[which(x[,"Exons_Skipped"] == 5),"Average_Coverage_NORM1"])
sN = log2(x[which(x[,"Exons_Skipped"] >= 6),"Average_Coverage_NORM1"])
data_list = list(s0,s1,s2,s3,s4,s5,sN)
names(data_list) = c("S0","S1","S2","S3","S4","S5","S6+")
tiff(file="ExonJunctionExpressionLevels_S0-S6.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
boxplot(data_list, col=rainbow(7), xlab="Number of Exons Skipped", ylab="Log2 Average Coverage (Expression)", 
        main="Exon junction expression level (by number of exons skipped)",
        col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()

#E-1.) Expressed junctions by exon skip #(S0, S1, S2, S3, etc.)  - Now all those junction whether expressed or not (add 1 to prevent -INF log2 values)
expressed_i = which(junctions[,"Average_Coverage_NORM1"] > -1 & junctions[,"Percent_Coverage_1x"] > -1)
x = junctions[expressed_i,]

s0 = log2(x[which(x[,"Exons_Skipped"] == 0),"Average_Coverage_NORM1"]+1)
s1 = log2(x[which(x[,"Exons_Skipped"] == 1),"Average_Coverage_NORM1"]+1)
s2 = log2(x[which(x[,"Exons_Skipped"] == 2),"Average_Coverage_NORM1"]+1)
s3 = log2(x[which(x[,"Exons_Skipped"] == 3),"Average_Coverage_NORM1"]+1)
s4 = log2(x[which(x[,"Exons_Skipped"] == 4),"Average_Coverage_NORM1"]+1)
s5 = log2(x[which(x[,"Exons_Skipped"] == 5),"Average_Coverage_NORM1"]+1)
s6 = log2(x[which(x[,"Exons_Skipped"] == 6),"Average_Coverage_NORM1"]+1)
s7 = log2(x[which(x[,"Exons_Skipped"] == 7),"Average_Coverage_NORM1"]+1)
s8 = log2(x[which(x[,"Exons_Skipped"] == 8),"Average_Coverage_NORM1"]+1)
s9 = log2(x[which(x[,"Exons_Skipped"] == 9),"Average_Coverage_NORM1"]+1)
s10 = log2(x[which(x[,"Exons_Skipped"] == 10),"Average_Coverage_NORM1"]+1)
sN = log2(x[which(x[,"Exons_Skipped"] >= 11),"Average_Coverage_NORM1"]+1)
data_list = list(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,sN)
names(data_list) = c("S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10", "S11+")
tiff(file="ExonJunctionExpressionLevels_ALL_JUNCTIONS_S0-S11.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
boxplot(data_list, col=rainbow(12), border=rainbow(12), xlab="Number of Exons Skipped", ylab="Log2 Average Coverage (Expression)", 
        main="Exon junction expression level (by number of exons skipped)",
        col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()

