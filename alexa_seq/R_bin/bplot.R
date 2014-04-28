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
#/home/malachig/svn/solexa_analysis/R_bin/bplot.R /projects/malachig/solexa/figures_and_stats/HS04391/temp/temp_data.txt 10 1000000 "Alignment Length" "bp" /projects/malachig/solexa/figures_and_stats/HS04391/ENST_v49/  

#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1]
column_count = as.numeric(args[2])
record_count = as.numeric(args[3])
data_name = args[4]
units = args[5]
output_dir = args[6]

bg = "white"

#Import the data

#NOTE: The input data file must be space delimited
what_list = vector(mode="list", column_count)
for (i in 1:column_count){
  what_list[[i]] = numeric(0)
}

print("Importing data with scan()")
data = scan(file=infile, what=what_list, nlines=record_count, fill=TRUE, sep=",")

#Experiment with filehash library
#library(filehash)
#setwd(output_dir)
#print("Importing data with scan() and dumping to database file with dumpDF()")
#data = scan(file=infile, what=what_list, nlines=record_count, fill=TRUE, sep=",")
#names(data) = c(1:column_count)
#db = dumpList(data,  dbName="db01")
#print ("Creating environment creating an environment to access the data on file")
#env = new.env()
#dbLoad(db, env)
#ls(env)


print("Creating boxplot with this data structure")
filename=paste(output_dir, gsub(" +", "", data_name), "_BoxPlot.tiff", sep="")
main_title = paste("Distribution of ", data_name, " values", sep="");
y_label = paste(data_name, " (", units, ")", sep="")
x_label = "Lane Number"
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg)
boxplot(x=data[1:length(data)], main=main_title, xlab=x_label, ylab=y_label, border=c(grey(0),grey(0.5)), col=c("blue","green"))
dev.off()



