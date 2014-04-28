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
#[1] dir
#[2] filebase

#Example usage:
#/home/malachig/svn/alexa_seq/R_bin/summarizeComplexity.R /projects/malachig/alexa_seq/figures_and_stats/HS04391/LibraryQuality/
#working_dir = "/projects/malachig/alexa_seq/figures_and_stats/HS04391/LibraryQuality/"
#filebase = "LibraryComplexity_SEQ"

args = (commandArgs(TRUE))
working_dir = args[1]
filebase = args[2]
infile = paste(working_dir, "temp.txt", sep="")
data=read.table(infile, sep="\t", na.strings=c("NA"), header=TRUE, as.is=c(1))

#Initiate three panel figure
title_text = paste("Library complexity metrics for ", data[1,"LibraryID"], " compared to ", length(data[,1]), " other libraries", sep="")
outfile1 = paste(filebase, ".jpeg", sep="")
setwd(working_dir)
Cairo(width=750, height=750, file=outfile1, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
par(mfrow=c(1,3), font.main=2, font.lab=2, font.axis=2, cex.axis=1.5, cex.lab=1.5, oma=c(4, 2, 4, 1))
legend_text=c(data[1,"LibraryID"],
              paste("Total libraries = ", length(data[,1]), sep=""))

#PercentUniqueReadsPerMillion
x=data[,"PercentUniqueReadsPerMillion"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Unique Reads", col="light green", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentUniqueReadsPerMillion"], col="red", pch=c(20), cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)

#PercentRedundantReadsPerMillion
x = data[,"PercentRedundantReadsPerMillion"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Redundant Reads", col="orange", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentRedundantReadsPerMillion"], col="red", pch=20, cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)

#PercentDistinctRedundantReads
x = data[,"PercentDistinctRedundantReads"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Distinct Redundant", col="yellow", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentDistinctRedundantReads"], col="red", pch=20, cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)
title(title_text, outer=TRUE,cex.main=2, font.main=2)
zz = dev.off()


#Same figure as an SVG
outfile2 = paste(filebase, ".svg", sep="")
Cairo(width=750, height=750, file=outfile2, type="svg", pointsize=12, bg="white", units="px", dpi=72)
par(mfrow=c(1,3), font.main=2, font.lab=2, font.axis=2, cex.axis=1.5, cex.lab=1.5, oma=c(4, 2, 4, 1))
legend_text=c(data[1,"LibraryID"],
              paste("Total libraries = ", length(data[,1]), sep=""))

#PercentUniqueReadsPerMillion
x=data[,"PercentUniqueReadsPerMillion"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Unique Reads", col="light green", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentUniqueReadsPerMillion"], col="red", pch=c(20), cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)

#PercentRedundantReadsPerMillion
x = data[,"PercentRedundantReadsPerMillion"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Redundant Reads", col="orange", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentRedundantReadsPerMillion"], col="red", pch=20, cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)

#PercentDistinctRedundantReads
x = data[,"PercentDistinctRedundantReads"]
min_x=min(x)
max_x=max(x)
diff_x=max_x-min_x
boxplot(x, ylab="Percent", ylim=c(min_x, max_x+(diff_x*0.1)), xlab="Distinct Redundant", col="yellow", pch=20, cex=1.5)
points(x=1, y=data[1,"PercentDistinctRedundantReads"], col="red", pch=20, cex=2)
legend("topleft", legend_text, pch=c(20,NA), col="red", cex=1.5)
title(title_text, outer=TRUE,cex.main=2, font.main=2)
zz = dev.off()



  
