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
#/home/malachig/svn/solexa_analysis/R_bin/pie.R /projects/malachig/solexa/figures_and_stats/HS04391/temp/temp_data.txt 1000000 "Hit Type" "Reads" /projects/malachig/solexa/figures_and_stats/HS04391/ENST_v49/  

#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1]
record_count = as.numeric(args[2])
data_name = args[3]
units = args[4]
output_dir = args[5]

bg = "white"

#Import the data
data = scan(file=infile, what=character(0), nmax=record_count)

#Display as a pie chart
filename=paste(output_dir, gsub(" +", "", data_name), "_Pie.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg)
main_title = paste("Counts of ", data_name, " values", sep="");
x_label = paste(units, " (n = ", record_count, ")", sep="");
pie(table(data), col=rainbow(length(table(data))), main=main_title, xlab=x_label)
dev.off()

#Display same data as a dot chart
filename=paste(output_dir, gsub(" +", "", data_name), "_Dot.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg)
main_title = paste("Counts of ", data_name, " values", sep="");
x_label = paste(data_name, " (n = ", record_count, ")", sep="");
dotchart(table(data), col=rainbow(length(table(data))), main=main_title, xlab=x_label, pch=17)
dev.off()

#Display the same data as a bar chart
filename=paste(output_dir, gsub(" +", "", data_name), "_Bar.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg)
main_title = paste("Counts of ", data_name, " values", sep="");
x_label = paste(data_name, sep="");
y_label = paste("Count (", units, "; n = ", record_count, ")", sep="") 
barplot(table(data), col=rainbow(length(table(data))), main=main_title, xlab=x_label, ylab=y_label)
dev.off()




