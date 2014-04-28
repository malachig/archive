#!/usr/bin/env Rscript
#Written by Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

library(Cairo)

#Options:
#[1] phred_dir
#[2] flowcell_lanes
#[3] result_dir 
#[4] library_name 

#Example usage:
#/scratch/obig/ALEXA/alexa_seq/R_bin/summarizePhred.R /csb/home/obig/ALEXA/analysis/read_records/600MPE/Summary/ /csb/home/obig/ALEXA/figures_and_stats/600MPE/Expression_v53/ "42JVHAAXX_6,42JVHAAXX_7" 600MPE

#Get command-line arguments
args = (commandArgs(TRUE))
phred_dir = args[1]
result_dir = args[2]
flowcell_lanes = args[3]
library_name = args[4]

#Get list of lanes to process
lane_list=strsplit(flowcell_lanes,split=",")[[1]]

#All lanes will be presented in a single figure
result_file = paste(result_dir,library_name,"_phred.jpeg", sep="")

#Go through each lane and combine data from different lanes
data_cum_list_R1=NULL
data_cum_list_R2=NULL

for (i in 1:length(lane_list)){
  print(lane_list[i])

  phred_file = paste(phred_dir,lane_list[i],"_phred.txt", sep="")
  data=read.table(phred_file, header=TRUE, na.strings = "NA")
  
  #Convert NAs to 0
  data[is.na(data[,"R1_count"]),"R1_count"]=0
  data[is.na(data[,"R2_count"]),"R2_count"]=0

  #Make sure there are actually some values to plot 
  distinct_values_R1=length(unique(data[,"R1_count"]))
  distinct_values_R2=length(unique(data[,"R2_count"]))

  if (distinct_values_R1>1 && distinct_values_R2>1){
    #Determine cumulative percent reads for each phred score
    data_cum_R1=0
    data_cum_R2=0
    R1_cum_dist=NULL
    R2_cum_dist=NULL
    phred_scores=data[,"phred_score"]
    for (i in 1:length(data[,"phred_score"])){
      data_cum_R1=data_cum_R1+data[i,"R1_count"]
      data_cum_R2=data_cum_R2+data[i,"R2_count"]
      R1_cum_dist[i]=data_cum_R1/sum(data[,"R1_count"])
      R2_cum_dist[i]=data_cum_R2/sum(data[,"R2_count"])
    }
    data_cum_R1=cbind(phred_scores,R1_cum_dist)
    data_cum_R2=cbind(phred_scores,R2_cum_dist)
    data_cum_list_R1=c(data_cum_list_R1,list(data_cum_R1))
    data_cum_list_R2=c(data_cum_list_R2,list(data_cum_R2))
  }
}

#Plot cumulative distribution for all lanes and both reads
setwd(result_dir)

#Set colors
bg="white"
colors=rainbow(8) #assumes never more than 8 lanes per library

#Create legend
legend_text = lane_list
legend_location = "topleft"
legend_colors = colors[1:length(lane_list)]

CairoJPEG(filename=result_file, width=750, height=750, pointsize=12, quality=100, bg=bg)

par(mfrow=c(1,2), bg=bg, font.main = 2, font.lab = 2, oma=c(4, 2, 4, 1))

#R1 - Check to make sure some data was actually found
if (length(data_cum_list_R1)>0){
  plot(x=c(min(data_cum_list_R1[[1]][,"phred_scores"]),max(data_cum_list_R1[[1]][,"phred_scores"])), y=c(0,1), type="n", main="R1", xlab="Mean Phred Score", ylab="Cumulative Fraction of Reads",cex.main = 1.4, cex.lab = 1.2)
  for (i in 1:length(data_cum_list_R1)){
    lines(x=data_cum_list_R1[[i]][,"phred_scores"], y=data_cum_list_R1[[i]][,"R1_cum_dist"], col=colors[i], lty=1, lwd=3)
  }
  legend(legend_location, legend=legend_text, lty=1, lwd=3, col=legend_colors)
}

#R2 - Check to make sure some data was actually found
if (length(data_cum_list_R2)>0){
  plot(x=c(min(data_cum_list_R2[[1]][,"phred_scores"]),max(data_cum_list_R2[[1]][,"phred_scores"])), y=c(0,1), type="n", main="R2", xlab="Mean Phred Score", ylab="Cumulative Fraction of Reads",cex.main = 1.4, cex.lab = 1.2)
  for (i in 1:length(data_cum_list_R2)){
    lines(x=data_cum_list_R2[[i]][,"phred_scores"], y=data_cum_list_R2[[i]][,"R2_cum_dist"], col=colors[i], lty=1, lwd=3)
  }
  legend(legend_location, legend=legend_text, lty=1, lwd=3, col=legend_colors)
}
#Create overall title for multi-plot figure
title_text=paste("Phred score distribution for",library_name,sep=" ")
title(title_text, outer=TRUE, cex.main=1.6)

dev.off()

