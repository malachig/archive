#!/usr/bin/env Rscript
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

library(Cairo)

#Options:
#[1] cutoffs_file
#[2] results_dir 

#Example usage:
#/home/malachig/svn/solexa_analysis/R_bin/summarizeCutoffs.R /projects/malachig/solexa/figures_and_stats/MM0472/Expression_v54/MM0472_NORM1_average_coverage_cutoffs.txt  /projects/malachig/solexa/figures_and_stats/MM0472/Expression_v54/

args = (commandArgs(TRUE))
datafile = args[1];
results_dir = args[2];

bg="white"

data=read.table(datafile, header=TRUE, na.strings = "N/A")

#Make sure there is actually some values to plot (in some libraries, all cutoffs are equal to the intergenic cutoff and therefore there is no distribution to plot)
distinct_values=length(unique(data[,"Cutoff"]))

#Note that most genes will have the miniumum cutoff (intergenic based) but some expressed genes will have a higher intragenic based cutoff
#Summarize the number of genes in each group but only plot those with a cutoff higher than the base
if (distinct_values > 5){
  min_cutoff = min(data[,"Cutoff"])
  max_cutoff = max(data[,"Cutoff"])
  upper = which(data[,"Cutoff"] > min_cutoff)
  n_all = length(data[,"Cutoff"])
  n_upper = length(upper)
  n_lower = n_all-n_upper

  legend_text = c("Intergenic cutoff",
                paste("Total number of genes = ", prettyNum(n_all, big.mark=","), sep=""),
                paste("Genes with intergenic cutoff = ", prettyNum(n_lower, big.mark=","), sep=""),
                paste("Genes with intragenic cutoff (summarized below) = ", prettyNum(n_upper, big.mark=","), sep=""))

  legend_location = "topleft"

  outfile1 = "DistributionOfCutoffs.jpeg"
  setwd(results_dir)
  CairoJPEG(filename=outfile1, width=750, height=750, pointsize=12, quality=100, bg=bg)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = log2(data[upper,"Cutoff"]), col="blue", main="Distribution of cutoff values", xlab="Intragenic cutoff value (log2)", ylab="Frequency", xlim=c(0,log2(max_cutoff)),
      col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks = 100)
  abline(v=log2(min_cutoff), col="red", lty=2, lwd=3)
  legend(legend_location, legend=legend_text, lty=c(2,NA,NA,NA), lwd=c(3,NA,NA,NA), col=c("red",NA,NA,NA))
  zz=dev.off()
}else{
  min_cutoff = data[1,"Cutoff"]
  max_cutoff = data[1,"Cutoff"]
  n_all = length(data[,"Cutoff"])

  legend_text = c("Intergenic cutoff",
                paste("Total number of genes = ", prettyNum(n_all, big.mark=","), sep=""),
                paste("Genes with intergenic cutoff = ", prettyNum(n_all, big.mark=","), sep=""))

  legend_location = "topleft"

  outfile1 = "DistributionOfCutoffs.jpeg"
  setwd(results_dir)
  CairoJPEG(filename=outfile1, width=750, height=750, pointsize=12, quality=100, bg=bg)
  par(bg=bg, font.main = 2, font.lab = 2)
  hist(x = log2(data[,"Cutoff"]), col="blue", main="Distribution of cutoff values", xlab="Cutoff values (log2)", ylab="Frequency", xlim=c(0,log2(max_cutoff)),
      col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, breaks = 100)
  abline(v=log2(min_cutoff), col="red", lty=2, lwd=3)
  legend(legend_location, legend=legend_text, lty=c(2,NA,NA), lwd=c(3,NA,NA), col=c("red",NA,NA))
  zz=dev.off()
}



