#Written by Malachi Griffith
#Using the following 64-bit installation of R
#/home/malachig/R64/R-2.9.0/bin/R

#Run this script from an xhost
#To allow graphs to open, you need to log in to an xhost with the -X option
#You can then use the X11 display to open graph windows as well as create tiffs and jpegs
#Alternatively you can write these straight to file and using the Cairo package you can do this without X11 forwarding

#Load neccessary libraries
library(geneplotter)
library(RColorBrewer)
library(Cairo)

#Example datafile
#/projects/luterlot/analysis/HS1374-HS1373.pileups.merge

#Open the file as a 'dataframe' - with large files, this will take a while... - for very large datasets you will need to use different methods to load the data
datafile = "/projects/luterlot/analysis/HS1374-HS1373.pileups.merge"
data=read.table(datafile, header = FALSE, na.strings = "NA", sep=" ")

#Set working directory for output
results_dir = "/home/malachig/"
setwd(results_dir)

#Assign names to each column of data
names(data) = c("position", "Lib_A_Depth", "Lib_B_Depth")

#Plot x vs. y
x = data[,"Lib_A_Depth"]
y = data[,"Lib_B_Depth"]

#Fit a linear model to the data (i.e. y = mx+b).  The coefficients of this fit will be used to show the fitted line on the graph
fit = lm(y~x);

#What's the r-squared value for this fit (by spearman method)
r_squared = (cor(x,y, method="spearman"))^2

#Text for the legend including the r-squared value
legend_text = c(paste("R-squared = ", round(r_squared, digits=2), " (Spearman)", sep=""))

#Where to put the legend
legend_location = "topleft"

#Define the colors to be used for the density scatter plot
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#Now make the actual plots
#Define an output file name, open a Cairo tiff connection, output graphics including legend, close the connection

#Data has an unsual structure... Plot all data and then try dividing data into two seperate partitions - say above vs. below 1000
partition = 1000
lower = which(data[,"Lib_A_Depth"] <= partition)
upper = which(data[,"Lib_A_Depth"] > partition)

#A.) All data
filename = "ScatterPlot_LibA_vs_LibB.tiff"

Cairo(width=800, height=800, file=filename, type="tiff", pointsize=12, bg="white", compression=1)
par(font.main = 2, font.lab = 2)
smoothScatter(x, y, xlab="Library A Depth", ylab="Library B Depth", main="Correlation between depth at each position of the genome", colramp=colors, nbin=275,
              col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
abline(coef=fit$coef, col="black", lty=2, lwd=2)
legend(legend_location, legend=legend_text, lty=c(2), lwd=2, col=c("black"))

#Indicate the box of data that will be considered in the next plot
lines(seq(0,partition,1), rep(partition,partition+1), lwd=1, col="orange")
lines(rep(partition,partition+1), seq(0,partition,1), lwd=1, col="orange")
lines(rep(0,partition+1), seq(0,partition,1), lwd=1, col="orange")
lines(seq(0,partition,1), rep(0,partition+1), lwd=1, col="orange")
dev.off()


#B.) Lower - for this plot limit the view to the 0-1000 data space to zoom in
x = data[lower,"Lib_A_Depth"]
y = data[lower,"Lib_B_Depth"]
filename = "ScatterPlot_LibA_vs_LibB_lower.tiff"
fit = lm(y~x);
r_squared = (cor(x,y, method="spearman"))^2
legend_text = c(paste("R-squared = ", round(r_squared, digits=2), " (Spearman)", sep=""))

Cairo(width=800, height=800, file=filename, type="tiff", pointsize=12, bg="white", compression=1)
par(font.main = 2, font.lab = 2)
smoothScatter(x, y, xlab="Library A Depth", ylab="Library B Depth", main="Correlation between depth at each position of the genome", colramp=colors, nbin=275,
              col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0, xlim=c(0,partition), ylim=c(0,partition))
abline(coef=fit$coef, col="black", lty=2, lwd=2)
legend(legend_location, legend=legend_text, lty=c(2), lwd=2, col=c("black"))
dev.off()


#C.) Upper
x = data[upper,"Lib_A_Depth"]
y = data[upper,"Lib_B_Depth"]
filename = "ScatterPlot_LibA_vs_LibB_upper.tiff"
fit = lm(y~x);
r_squared = (cor(x,y, method="spearman"))^2
legend_text = c(paste("R-squared = ", round(r_squared, digits=2), " (Spearman)", sep=""))

Cairo(width=800, height=800, file=filename, type="tiff", pointsize=12, bg="white", compression=1)
par(font.main = 2, font.lab = 2)
smoothScatter(x, y, xlab="Library A Depth", ylab="Library B Depth", main="Correlation between depth at each position of the genome", colramp=colors, nbin=275,
              col.lab = gray(.1), cex.main = 1.4, cex.lab = 1.2, cex.axis=1.0)
abline(coef=fit$coef, col="black", lty=2, lwd=2)
legend(legend_location, legend=legend_text, lty=c(2), lwd=2, col=c("black"))
dev.off()






















