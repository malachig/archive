#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Test sets of coords
exons = as.data.frame(list(name=c("e1","e2","e3"), start=c(10,43,1000), end=c(20,100,1200)))
reads = as.data.frame(list(name=c("r1","r2","r3"), start=c(15,30,240), end=c(45,60,270)))

exons = as.data.frame(list(start=c(10,43,1000), end=c(20,100,1200)))
reads = as.data.frame(list(start=c(15,30,240), end=c(45,60,270)))
rownames(exons) = c("e1","e2","e3")
rownames(reads) = c("r1","r2","r3")

exons
reads

#Find the exons that each read overlaps
#Answer.  r1 overlaps e1,e2; r2 overlap e2; r3 overlaps nothing

#Method 1.  Using parallel maxima and minima
#Does not tell you when something in list1 matches multiple things in list2
#Not sure about importance of BOTH lists being ordered.  This is not always possible in the case of overlapping sets of coordinates WITHIN each list
overlaping = function(list1,list2){
  Intersec = data.frame(start=pmax(list1$start, list2$start), end=pmin(list1$end, list2$end)) 
  i = which(Intersec$start <= Intersec$end)
  list1[i,]
}
overlaping(exons,reads)


#Method 2.  Using intervals package
install.packages("intervals")
library(intervals)

exons = Intervals(matrix(c(10,20,43,100,1000,1200), byrow = TRUE, ncol = 2), closed = c( TRUE, TRUE ), type = "Z")
rownames(exons) = c("e1", "e2", "e3")
reads = Intervals(matrix(c(15,45,30,60,240,270), byrow = TRUE, ncol = 2), closed = c( TRUE, TRUE ), type = "Z")
rownames(reads) = c("r1", "r2", "r3")
interval_overlap(exons,reads)

as.matrix(interval_overlap(exons,reads))

#The 'which_nearest' function seems to return a nicer output
which_nearest(exons, reads)
cbind(reads, which_nearest( reads, exons ))
x = cbind(reads, which_nearest( reads, exons ))

write.table(x, "test.txt", sep="\t")

#Now try with input files containing coords
setwd("/projects/malachig/solexa/temp/")
r_data = read.table("chr1_reads.txt")
e_data = read.table("chr1_exons.txt")

exons = Intervals(matrix(as.matrix(e_data), byrow = TRUE, ncol = 2), closed = c( TRUE, TRUE ), type = "Z")
reads = Intervals(matrix(as.matrix(r_data), byrow = TRUE, ncol = 2), closed = c( TRUE, TRUE ), type = "Z")
x = cbind(reads, which_nearest( reads, exons ))
write.table(x, "test.txt", sep="\t")




#Method 3.  Using genomeIntervals package
source("http://bioconductor.org/biocLite.R")
biocLite("genomeIntervals")


