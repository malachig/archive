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
#/home/malachig/svn/solexa_analysis/R_bin/positionBias.R /projects/malachig/solexa/figures_and_stats/HS04391/temp/temp_data.txt 10 1000000 "Alignment Length" "bp" /projects/malachig/solexa/figures_and_stats/HS04391/ENST_v49/  

#Get user specified arguments
args = (commandArgs(TRUE))
infile = args[1]
column_count = as.numeric(args[2])
record_count = as.numeric(args[3])
data_name = args[4]
units = args[5]
output_dir = args[6]

#Debug
#infile = "/projects/malachig/solexa/figures_and_stats/HS04391/temp/temp_data.txt"
#column_count = 9
#record_count = 902035
#data_name = "MIP101 Library - Position Bias"
#units = "Percent Position"
#output_dir = "/projects/malachig/solexa/figures_and_stats/HS04391/ENST_v53/"

bg = "white"

#Import the data

#NOTE: The input data file must be space delimited
what_list = vector(mode="list", column_count)
for (i in 1:column_count){
  what_list[[i]] = numeric(0)
}
print("Importing data with scan()")

#Import only the header line
header = scan(file=infile, nlines=1, sep=",")

#Import the data but skip the header line
data = scan(file=infile, what=what_list, nlines=record_count, fill=TRUE, sep=",", skip=1)

#Assign names from header to the data object
names(data) = header
names(data) = prettyNum(names(data), big.mark=",")

#Draw a box plot for each bin of data
colors = rainbow(column_count)
print("Creating boxplot with this data structure")
filename=paste(output_dir, gsub(" +", "", data_name), "_BoxPlot.tiff", sep="")
main_title = paste("Distribution of ", data_name, " values", sep="");
y_label = paste(data_name, " (", units, ")", sep="")
x_label = "Transcript Size (bp; Upper limit of bin)"
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
boxplot(x=data[1:length(data)], main=main_title, xlab=x_label, ylab=y_label, col=colors,
        col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
dev.off()

#Draw a box plot for each bin of data - using the 'varwidth' option
filename=paste(output_dir, gsub(" +", "", data_name), "_BoxPlot_VarWidth.tiff", sep="")
main_title = paste("Distribution of ", data_name, " values", sep="");
y_label = paste(data_name, " (", units, ")", sep="")
x_label = "Transcript Size (bp; Upper limit of bin)"
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
boxplot(x=data[1:length(data)], varwidth=TRUE, main=main_title, xlab=x_label, ylab=y_label, col=colors,
        col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
dev.off()

#Get max observed frequencies
max_observed = 0
max_observed2 = 0
com = 0
for (i in 1:column_count){
  x=max(table(data[i]))
  x2=max(table(round(data[[i]], digits=0)))
  if (x > max_observed){
    max_observed=x
  }
  if (x2 > max_observed2){
    max_observed2=x2
  }
  com=c(com,data[[i]])
}

#Draw as a series of line plots
filename=paste(output_dir, gsub(" +", "", data_name), "_Lines.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=names(table(data[1])),y=table(data[1]), type="l", xlim=c(0,105), ylim=c(0,max_observed), col=colors[1],
    xlab="Relative read position", ylab="Frequency", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 2:column_count){
  z=table(data[i])
  x=names(z)
  y=z
  lines(x=x, y=y, col=colors[i], lwd=2)
}
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", names(data), lty=1, lwd=2, col=colors)
dev.off()

#Draw as a series of line plots but use smoothing
span=0.2 
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesSmooth.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=names(table(data[1])),y=table(data[1]), type="n", xlim=c(0,105), ylim=c(0,max_observed), col=colors[1],
    xlab="Relative read position", ylab="Frequency", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 1:column_count){
  z=table(data[i])
  x=names(z)
  y=z
  lines(loess.smooth(x,y,span=span),col=colors[i],lwd=2)
}
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", names(data), lty=1, lwd=2, col=colors)
dev.off()

#Try using hist function to get values to plot
breaks=100
max_count=0
for (i in 1:column_count){
  h=hist(data[[i]], breaks=breaks, plot=FALSE)
  max=max(h$counts)
  if(max>max_count){max_count=max}
}  
filename=paste(output_dir, gsub(" +", "", data_name), "_HistLines.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=names(table(data[1])),y=table(data[1]), type="n", xlim=c(0,115), ylim=c(0,max_count), col=colors[1],
    xlab="Relative read position", ylab="Frequency", col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 1:column_count){
  h=hist(data[[i]], breaks=breaks, plot=FALSE)
  x=h$mids
  y=h$counts
  lines(x=x, y=y, col=colors[i], lwd=2)
}
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c(names(data)), lty=1, lwd=2, col=c(colors))
dev.off()

#Simpler smoothing strategy
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesRounded.tiff", sep="")
grand_total=sum(table(round(com, digits=0)))
all_x=names(table(round(com, digits=0)))
all_y=table(round(com, digits=0))
 
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x,y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed2), col="black",
    xlab="Relative read position", ylab="Relative frequency", lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 1:column_count){
  z=table(round(data[[i]], digits=0))
  x=names(z)
  y=z
  lines(x=x, y=y, col=colors[i], lwd=2)
}
lines(x=all_x,y=all_y, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined",names(data)), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
dev.off()



#Try scaling all the frequencies to see what the pattern would look like if each bin had equal numbers
max_observed3 = 0.25
filename=paste(output_dir, gsub(" +", "", data_name), "_LinesScaled.tiff", sep="")
grand_total=sum(table(com))
all_x=names(table(com))
all_y=(table(com)/grand_total)*100
span=0.2 
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=all_x,y=all_y, type="n", xlim=c(0,105), ylim=c(0,max_observed3), col="black",
    xlab="Relative read position", ylab="Relative frequency", lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
for (i in 1:column_count){
  z=table(data[i])
  total=sum(z)
  x=names(z)
  y=(z/total)*100
  lines(loess.smooth(x,y,span=span),col=colors[i],lwd=2)
}
lines(loess.smooth(all_x,all_y,span=span), col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined",names(data)), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
dev.off()


#Try scaling all the frequencies to see what the pattern would look like if each bin had equal numbers
max_observed4 = 0
for (i in 1:column_count){
  hi=hist(data[[i]], breaks=breaks, plot=FALSE)
  scaled_counts=(hi$counts/sum(hi$counts))*100
  max_y = max(scaled_counts)
  if(max_y > max_observed4){
    max_observed4 = max_y
  }
}
max_observed4=max_observed4+(max_observed4*0.1)
filename=paste(output_dir, gsub(" +", "", data_name), "_HistLinesScaled.tiff", sep="")
breaks=100
h=hist(com, breaks=breaks, plot=FALSE)
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=h$mids,y=h$counts, type="n", xlim=c(0,105), ylim=c(0,max_observed4), col="black",
    xlab="Relative read position", ylab="Relative frequency", lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
alternates=seq(1,column_count,2)
for (i in 1:column_count){
  hi=hist(data[[i]], breaks=breaks, plot=FALSE)
  scaled_counts=(hi$counts/sum(hi$counts))*100
  lines(x=hi$mids, y=scaled_counts, col=colors[i],lwd=2, lty=1)
}
scaled_counts=(h$counts/sum(h$counts))*100
lines(x=h$mids, y=scaled_counts, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined",names(data)), lty=c(2,rep(1,column_count)), lwd=2, col=c("black",colors))
dev.off()

#Again but this time only show alternate bins
max_observed5 = 0
for (i in alternates){
  hi=hist(data[[i]], breaks=breaks, plot=FALSE)
  scaled_counts=(hi$counts/sum(hi$counts))*100
  max_y = max(scaled_counts)
  if(max_y > max_observed5){
    max_observed5 = max_y
  }
}
max_observed5=max_observed5+(max_observed5*0.1)
filename=paste(output_dir, gsub(" +", "", data_name), "_HistAlternateLinesScaled.tiff", sep="")
tiff(file=filename, width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2)
plot(x=h$mids,y=h$counts, type="n", xlim=c(0,105), ylim=c(0,max_observed5), col="black",
    xlab="Relative read position", ylab="Relative frequency", lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.25)
alternates=seq(1,column_count,2)
for (i in alternates){
  hi=hist(data[[i]], breaks=breaks, plot=FALSE)
  scaled_counts=(hi$counts/sum(hi$counts))*100
  lines(x=hi$mids, y=scaled_counts, col=colors[i],lwd=2, lty=1)
}
scaled_counts=(h$counts/sum(h$counts))*100
lines(x=h$mids, y=scaled_counts, col="black",lwd=3, lty=2)
abline(h=0,col="black",lty=3,lwd=1)
legend("topright", c("Combined",names(data[alternates])), lty=c(2,rep(1,length(alternates))), lwd=2, col=c("black",colors[alternates]))
dev.off()




#ks.test(x=sample(na.omit(data[["500"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["1000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["2000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["3000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["4000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["5000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["10000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["15000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["20000"]]), 10000), y=sample(na.omit(data[["2000"]]), 10000), alternative = "two.sided")$p.value

#ks.test(x=sample(na.omit(data[["500"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["1000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["2000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["3000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["4000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["5000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["10000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["15000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value
#ks.test(x=sample(na.omit(data[["20000"]]), 10000), y=sample(na.omit(data[["3000"]]), 10000), alternative = "two.sided")$p.value

#median(na.omit(data[["500"]]))
#median(na.omit(data[["1000"]]))
#median(na.omit(data[["2000"]]))
#median(na.omit(data[["3000"]]))
#median(na.omit(data[["4000"]]))
#median(na.omit(data[["5000"]]))
#median(na.omit(data[["10000"]]))
#median(na.omit(data[["15000"]]))
#median(na.omit(data[["20000"]]))


