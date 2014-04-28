#!/usr/bin/env Rscript
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Example usage:
#/home/malachig/svn/alexa_seq/R_bin/FragmentSizeDistribution.R /projects/malachig/alexa_seq/figures_and_stats/HS04391/temp/temp_fragsize_data.txt /projects/malachig/alexa_seq/figures_and_stats/HS04391/ENST_v53/  
library(Cairo)

#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1]
output_dir = args[2]

#Debug
#infile = "/projects/malachig/alexa_seq/figures_and_stats/HS04391/temp/temp_data.txt"
#output_dir = "/projects/malachig/alexa_seq/figures_and_stats/HS04391/ENST_v53/"

main_title = "Fragment size distribution for paired-end reads"
bg = "white"

#Import the data but skip the header line
data = read.table(file=infile, sep="\t", header=TRUE, check.names=FALSE)

filename=paste(output_dir, "FragmentSize.jpeg", sep="")

#Convert NAs (if any) to 0
data[is.na(data[,"count"]),"count"]=0
perc=(data[,"count"]/sum(data[,"count"]))*100

#Make sure there are actually some values to plot 
distinct_values=length(unique(data[,"count"]))

if (distinct_values>1){
  #Filter out fragment sizes observed only once
  data=data[perc>0.01,]
  #Plot the distribution
  Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
  plot(x=data[,"fragment_size"], y=data[,"count"], type="l", xlab="Fragment Size (bp)", ylab="Count", main=main_title, col="blue",col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
  zz=dev.off()
}

