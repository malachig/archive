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
#/home/malachig/svn/alexa_seq/R_bin/positionBias.R /projects/malachig/alexa_seq/figures_and_stats/HS04391/temp/temp_data.txt "Alignment Length" "bp" /projects/malachig/alexa_seq/figures_and_stats/HS04391/ENST_v53/  
library(Cairo)

#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1]
output_dir = args[2]

#Debug
#infile = "/projects/malachig/alexa_seq/figures_and_stats/HS04391/temp/temp_data.txt"
#output_dir = "/projects/malachig/alexa_seq/figures_and_stats/HS04391/ENST_v53/"

data_name = "PositionBias"
main_title = "Relative position of reads mapping within transcripts of increasing size"
bg = "white"

#Import the data but skip the header line
data = read.table(file=infile, sep=" ", header=TRUE, check.names=FALSE)

#Assign names from header to the data object
names(data) = prettyNum(names(data), big.mark=",")
bin_names = names(data)[2:10]

column_count=length(data[1,])
colors = rainbow(column_count-1)

#Get max observed frequencies
max_observed = 0
for (i in 2:column_count){
  x=max(data[,i])
  if (x > max_observed){
    max_observed=x
  }
}


#1.) Draw as a series of line plots
filename=paste(output_dir, gsub(" +", "", data_name), "_Lines.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=data[,1], y=data[,2], type="l", xlim=c(0,105), ylim=c(0,max_observed), col=colors[1],
    xlab="Relative read position", ylab="Frequency", main=main_title, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 3:column_count){
  x=(data[,1])  
  y=data[,i]
  lines(x=x, y=y, col=colors[i-1], lwd=2)
}
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", names(data[2:column_count]), lty=1, lwd=2, col=colors)
zz=dev.off()


#2.) Now use loess to smooth these lines
span=0.2 
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesLoess.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=data[,1], y=data[,2], type="n", xlim=c(0,105), ylim=c(0,max_observed), col=colors[1],
    xlab="Relative read position", ylab="Frequency (loess smoothed)", main=main_title, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 3:column_count){
  x=data[,1]
  y=data[,i]
  lines(loess.smooth(x,y,span=span),col=colors[i-1],lwd=2)
}
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", names(data[2:column_count]), lty=1, lwd=2, col=colors)
zz=dev.off()


#Get the sum of each row (i.e. all size bins combined)
data[,"SUM"]=apply(data[,2:10], 1, sum)


#3.) Now try scaling all the frequencies to see what the pattern would look like if each bin had equal numbers
max_observed = 0
for (i in 2:(column_count+1)){
  grand_total=sum(data[,i])
  z=(data[,i]/grand_total)*100
  x=max(z)
  if (x > max_observed){
    max_observed=x
  }
}
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaled.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
all_y=(data[,"SUM"]/grand_total)*100
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 2:column_count){
  total=sum(data[,i])
  x=data[,1]
  y=(data[,i]/total)*100
  lines(x=x, y=y,col=colors[i-1],lwd=2)
}
lines(x=all_x, y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
zz=dev.off()


#4.) Now try scaling and smoothing with loess
max_observed=0.25
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaledLoess.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
all_y=(data[,"SUM"]/grand_total)*100
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency (loess smoothed)", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 2:column_count){
  total=sum(data[,i])
  x=data[,1]
  y=(data[,i]/total)*100
  lines(loess.smooth(x,y,span=span),col=colors[i-1],lwd=2)
}
lines(loess.smooth(all_x,all_y,span=span), col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
zz=dev.off()


#5.) Now try a running median
kval=21
max_observed = 0
for (i in 2:(column_count+1)){
  grand_total=sum(data[,i])
  z=(data[,i]/grand_total)*100
  z2=runmed(z, k=kval)
  x=max(z2)
  if (x > max_observed){
    max_observed=x
  }
}
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaledRunMed.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
z=(data[,"SUM"]/grand_total)*100
all_y=runmed(z, k=kval)

par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency (running median)", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 2:column_count){
  total=sum(data[,i])
  x=data[,1]
  z=(data[,i]/total)*100
  y=runmed(z, k=kval)
  lines(x=x, y=y,col=colors[i-1],lwd=2)
}
lines(x=all_x, y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
zz=dev.off()

filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaledRunMed.svg", sep="")
Cairo(width=750, height=750, file=filename, type="svg", pointsize=12, bg="white", units="px", dpi=72)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
z=(data[,"SUM"]/grand_total)*100
all_y=runmed(z, k=kval)

par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency (running median)", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 2:column_count){
  total=sum(data[,i])
  x=data[,1]
  z=(data[,i]/total)*100
  y=runmed(z, k=kval)
  lines(x=x, y=y,col=colors[i-1],lwd=2)
}
lines(x=all_x, y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
zz=dev.off()


#6.) Now try a running median and display only the alternating bins
alternates=c(seq(2,column_count,2),11)

max_observed = 0
for (i in alternates){
  grand_total=sum(data[,i])
  z=(data[,i]/grand_total)*100
  z2=runmed(z, k=kval)
  x=max(z2)
  if (x > max_observed){
    max_observed=x
  }
}
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaledRunMedAlternates.jpeg", sep="")
Cairo(width=750, height=750, file=filename, type="jpeg", pointsize=12, bg="white", units="px", dpi=72, quality=100)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
z=(data[,"SUM"]/grand_total)*100
all_y=runmed(z, k=kval)

par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in alternates){
  total=sum(data[,i])
  x=data[,1]
  z=(data[,i]/total)*100
  y=runmed(z, k=kval)
  lines(x=x, y=y,col=colors[i-1],lwd=2)
}
lines(x=all_x, y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names[alternates[1:5]-1]), lty=c(2,rep(1,length(alternates))), lwd=2, col=c("black",colors[alternates[1:5]-1]))
zz=dev.off()

filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaledRunMedAlternates.svg", sep="")
Cairo(width=750, height=750, file=filename, type="svg", pointsize=12, bg="white", units="px", dpi=72)
grand_total=sum(data[,"SUM"])
all_x=data[,1]
z=(data[,"SUM"]/grand_total)*100
all_y=runmed(z, k=kval)

par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x, y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed), col="black",
    xlab="Relative read position", ylab="Relative frequency", main=main_title, lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in alternates){
  total=sum(data[,i])
  x=data[,1]
  z=(data[,i]/total)*100
  y=runmed(z, k=kval)
  lines(x=x, y=y,col=colors[i-1],lwd=2)
}
lines(x=all_x, y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined", bin_names[alternates[1:5]-1]), lty=c(2,rep(1,length(alternates))), lwd=2, col=c("black",colors[alternates[1:5]-1]))
zz=dev.off()

