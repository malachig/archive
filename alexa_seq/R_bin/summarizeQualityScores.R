#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#R script that works in conjunction with summarizeQualityScores.pl

#HS04391 LANE-BY-LANE PLOTS
setwd("/projects/malachig/solexa/maq_analysis/HS04391/quality_stats/results/")

#1.) Import Illumina base quality data according to read position - broken down for all lanes of data
#Get the average read score for each position for every lane
pos_data_r12 = read.table("R12_lanes.pos_by_data.txt", header=TRUE, na.strings = "NA")

bg = "white"

names=(1:((length(pos_data_r12[1,])-1)))
cols=c("green","light blue")
title = paste("MIP101 - Read quality (R1 + R2) - position variation "," (", (length(pos_data_r12[1,])-1)/2, " lanes)",  sep = "")
tiff(file="MIP101_QualityByLane_R12.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(x=pos_data_r12[,2:length(pos_data_r12[1,])], names=names, main=title, xlab="Lane (R1 then R2)", ylab="Average Illumina Qualities (for each read position)", col=cols, col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()


#HS04401 LANE-BY-LANE PLOTS
setwd("/projects/malachig/solexa/maq_analysis/HS04401/quality_stats/results/")

#1.) Import Illumina base quality data according to read position - broken down for all lanes of data
#Get the average read score for each position for every lane
pos_data_r12 = read.table("R12_lanes.pos_by_data.txt", header=TRUE, na.strings = "NA")

names=(1:((length(pos_data_r12[1,])-1)))
cols=c("green","light blue")
title = paste("MIP/5FU Read quality (R1 + R2) - position variation "," (", (length(pos_data_r12[1,])-1)/2, " lanes)",  sep = "")
tiff(file="MIP5FU_QualityByLane_R12.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(x=pos_data_r12[,2:length(pos_data_r12[1,])], names=names, main=title, xlab="Lane (R1 then R2)", ylab="Average Illumina Qualities (for each read position)", col=cols, col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()



#ALL LANES DROM BOTH LIBRARIES
setwd("/projects/malachig/solexa/maq_analysis/combined_stats/results/")

dir()

#1.) Import Illumina base quality data according to read position - broken down for all lanes of data
#Get the average read score for each position for every lane
data_pos_r12 = read.table("R12_lanes.data_by_pos.txt", header=TRUE, na.strings = "NA")
data_pos_r1 = read.table("R1_lanes.data_by_pos.txt", header=TRUE, na.strings = "NA")
data_pos_r2 = read.table("R2_lanes.data_by_pos.txt", header=TRUE, na.strings = "NA")

#2.) Examine the relationship between read position (i.e. beginning to end of read) and the average quality of bases at that position
#create a plot showing the variability in average score across all lanes of data at each read position (where R1 and R2 are considered together)
tiff(file="QualityVsPosition_R12combined.tiff", width=700, height=700, compression="none")
title = paste("Illumina read quality (R1 + R2 combined) - lane variation "," (", length(data_pos_r12[,1])/2, " lanes)",  sep = "")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(data_pos_r12[,2:length(data_pos_r12[1,])], names=(1:(length(data_pos_r12[1,])-1)), main=title, xlab="Read Position", ylab="Average Illumina Qualities (for each lane)", col="green", col.main = "black", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
dev.off()

#3.) Now do the same thing but consider R1 and R2 seperately... 
tiff(file="QualityVsPosition_R12seperate.tiff", width=700, height=700, compression="none")
names=(1:((length(data_pos_r1[1,])-1) + (length(data_pos_r2[1,])-1)))
cols = c(rep("green", length(data_pos_r1[1,])-1), rep("light blue", length(data_pos_r2[1,])-1))
title = paste("Illumina read quality (R1 then R2) - lane variation "," (", length(data_pos_r12[,1])/2, " lanes)",  sep = "")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(x=c(data_pos_r1[,2:length(data_pos_r1[1,])], data_pos_r2[,2:length(data_pos_r2[1,])]), names=names, main=title, xlab="Read Position", ylab="Average Illumina Qualities (for each lane)", col=cols, col.main = "black", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
dev.off()


#4.) Now import data representing the distribution of Illumina base quality scores (% of bases with quality score of 1, 2, ..., 40)
data_qual_r12 = read.table("R12_lanes.data_by_qual.txt", header=TRUE, na.strings = "NA")
data_qual_r1 = read.table("R1_lanes.data_by_qual.txt", header=TRUE, na.strings = "NA")
data_qual_r2 = read.table("R2_lanes.data_by_qual.txt", header=TRUE, na.strings = "NA")

#5.) Examine the distribution of base qualities (for R12, then just R1 then just R2)
tiff(file="BaseQualityDistribution_R12combined.tiff", width=700, height=700, compression="none")
title = paste("Illumina bases quality distribution (R1 + R2 combined) - lane variation "," (", length(data_qual_r12[,1])/2, " lanes)",  sep = "")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(data_qual_r12[,2:length(data_qual_r12[1,])], names=(1:(length(data_qual_r12[1,])-1)), main=title, xlab="Quality Score", ylab="Percentage of Bases (for each lane)", col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()

tiff(file="BaseQualityDistribution_R1.tiff", width=700, height=700, compression="none")
title = paste("Illumina bases quality distribution (R1) - lane variation "," (", length(data_qual_r1[,1]), " lanes)",  sep = "")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(data_qual_r1[,2:length(data_qual_r1[1,])], names=(1:(length(data_qual_r1[1,])-1)), main=title, xlab="Quality Score", ylab="Percentage of Bases (for each lane)", col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()

tiff(file="BaseQualityDistribution_R2.tiff", width=700, height=700, compression="none")
title = paste("Illumina bases quality distribution (R2) - lane variation "," (", length(data_qual_r2[,1]), " lanes)",  sep = "")
par(bg=bg, font.main = 2, font.lab = 2); boxplot(data_qual_r2[,2:length(data_qual_r2[1,])], names=(1:(length(data_qual_r2[1,])-1)), main=title, xlab="Quality Score", ylab="Percentage of Bases (for each lane)", col="blue", col.main = "black", col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
dev.off()






