#!/usr/bin/env Rscript
#Written by Malachi Griffith

#Load R libraries needed
library(Cairo)

#Example usage:
#[1] Gene name
#[2] Input file containing gene based coordinates 
#[3] CDS Start Position
#[4] CDS End Position
#[5] Output file

#Example data
#gene_name="CA12"
#infile3="/Users/malachig/Desktop/R_bin/CA12.Features.txt"
#cds_start=160
#cds_end=55595
#outfile="CA12.svg"


#Load infile and other user specified arguments
args = (commandArgs(TRUE))
gene_name = args[1];
infile3 = args[2];
cds_start = as.numeric(args[3])
cds_end = as.numeric(args[4])
outfile= args[5]

#Example command:
#/home/malachig/svn/alexa_seq/R_bin/drawGeneModel.R   CA12  /home/malachig/svn/alexa_seq/R_bin/CA12.Features.txt  160  55595  CA12.svg


#Data needed to create simple gene model visualization
#Start   End     Name    Type    IsAE    SkippedExons
model=read.table(infile3, header=TRUE, sep="\t", na.strings=c("NA"), as.is=c(3,4))
model[,"Size"]=(model[,"End"]-model[,"Start"])+1

#Get indexes of Exons, Introns, Junctions
ei = which(model[,"Type"]=="ExonRegion")
ii = which(model[,"Type"]=="Intron")
ji = which(model[,"Type"]=="Junction")
eii = which(model[,"Type"]=="ExonRegion" | model[,"Type"]=="Intron")
ae_i = which(model[,"IsAE"]==1)
ae_ji = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1)
ae_ji_can = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1 & model[,"SkippedExons"]==0)
ae_ji_skip = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1 & model[,"SkippedExons"]>0)

#Each gene will have two models, one drawn to scale at the top and another drawn on log2/log10 scale in the center
#'Draw1' will be the scale drawing at the top
#'Draw2' will be the main drawing in the center

#Create scaled coordinates, modify sizes as follows:
#exons on log2 scale
#introns on log10 scale
model[,"Log2Size"]=log2(model[,"Size"])
model[,"Log10Size"]=log10(model[,"Size"])
model[,"DrawSize2"]=NA
model[ei,"DrawSize2"]=model[ei,"Log2Size"]
if (length(ii) > 0){
 model[ii,"DrawSize2"]=model[ii,"Log10Size"]
}
if (length(ji) > 0){
 model[ji,"DrawSize2"]=model[ji,"Log2Size"]
}

#What is the total length needed?
grand_length=sum(model[eii,"DrawSize2"])

#Get draw size1 relative to draw size2
grand_length_scale=sum(model[eii,"Size"])
model[,"DrawSize1"]=(model[,"Size"]/grand_length_scale)*grand_length

#Define some extra width for the plot
extra_width_ratio=0.05
extra_width=ceiling(grand_length*extra_width_ratio)

#Define some dummy values to help create a blank plot of the correct size
xrange=0:ceiling(grand_length+(extra_width*2))
yrange=rep(5,length(xrange))

#Width is the cumulative size of all exons/introns plus some extra width
#Height is a constant defined below
height=10
y_inset=3
center_main=0+(height/y_inset)
center_scale=height-(height/y_inset)

#How many features (exons+introns) are to be plotted
feature_count=length(eii)
exon_count=length(ei)

#Define the y positions to be used for various parts of the plot
#All y positions are in an arbitrary space of range 0-10
#'main' refers to the y positions of the main gene model
#'scale' refers to the y positions of the to-scale model
#'j' refers to the y positions of junction connecting lines
#'label' refers to the y positions of text labels
height_ratio_main=0.05
height_ratio_scale=0.05

height_main=(height*height_ratio_main)
height_scale=(height*height_ratio_scale)
bottom_main=(center_main)-height_main
bottom_scale=(center_scale)-height_scale
top_main=(center_main)+height_main
top_scale=(center_scale)+height_scale
bottom_main_j=(center_main)-(height_main*2)
bottom_label_lower=(center_main)-(height_main*5)
top_main_j=(center_main)+(height_main*2)
bottom_label_main=(center_main)-(height_main*2)
top_label_main=(center_main)+(height_main*2)
top_label_scale=(center_scale)+(height_scale*5)

#Define the y coordinates of the gene model to be drawn
bottoms_main=rep(bottom_main, exon_count)
bottoms_scale=rep(bottom_scale, exon_count)
tops_main=rep(top_main, exon_count)
tops_scale=rep(top_scale, exon_count)

#Define the x coordinates of the gene model to be drawn
#First get the 'to-scale' drawing x coordinates
model[,"DrawStart1"]=NA
model[eii,"DrawStart1"]=(cumsum(model[eii,"DrawSize1"]))-(model[eii,"DrawSize1"])
model[eii,"DrawStart1"]=model[eii,"DrawStart1"]+extra_width
model[,"DrawEnd1"]=NA
model[eii,"DrawEnd1"]=cumsum(model[eii,"DrawSize1"])
model[eii,"DrawEnd1"]=model[eii,"DrawEnd1"]+extra_width
model[,"DrawMid1"]=NA
model[eii,"DrawMid1"]=model[eii,"DrawStart1"]+((model[eii,"DrawEnd1"]-model[eii,"DrawStart1"])/2)

#Next get the mixed log2/log10 scale drawing x coordinates
model[,"DrawStart2"]=NA
model[eii,"DrawStart2"]=(cumsum(model[eii,"DrawSize2"]))-(model[eii,"DrawSize2"])
model[eii,"DrawStart2"]=model[eii,"DrawStart2"]+extra_width
model[,"DrawEnd2"]=NA
model[eii,"DrawEnd2"]=cumsum(model[eii,"DrawSize2"])
model[eii,"DrawEnd2"]=model[eii,"DrawEnd2"]+extra_width
model[,"DrawMid2"]=NA
model[eii,"DrawMid2"]=model[eii,"DrawStart2"]+((model[eii,"DrawEnd2"]-model[eii,"DrawStart2"])/2)

#Find the location of the cds start and end
first_coding_exon=which((cds_start >= model[,"Start"]) & (cds_start <= model[,"End"]) & model[,"Type"]=="ExonRegion")
last_coding_exon=which((cds_end >= model[,"Start"]) & (cds_end <= model[,"End"]) & model[,"Type"]=="ExonRegion")
cds_draw_start=0
cds_draw_end=0
if(length(first_coding_exon)){
 relative_pos1=(model[first_coding_exon,"Size"] - (model[first_coding_exon,"End"]-cds_start))/model[first_coding_exon,"Size"]
 cds_draw_start=model[first_coding_exon,"DrawStart2"]+model[first_coding_exon,"DrawSize2"]*relative_pos1
}
if(length(last_coding_exon)){
 relative_pos2=(model[last_coding_exon,"Size"] - (model[last_coding_exon,"End"]-cds_end))/model[last_coding_exon,"Size"]
 cds_draw_end=model[last_coding_exon,"DrawStart2"]+model[last_coding_exon,"DrawSize2"]*relative_pos2
}

#Get the drawing x positions for only those junctions that were AE
#Remember that the 'start' of a junction is the end of the first exon, and 'end' is the start of the second exon
for(i in ae_ji){
  #Start of junction line
  j_start=(model[i,"Start"])	
  start_i=which(model[,"End"]==j_start)
  draw_start=(model[start_i,"DrawEnd2"])
  model[i,"DrawStart2"]=draw_start

  #End of junction line		
  j_end=(model[i,"End"])
  end_i=which(model[,"Start"]==j_end)
  draw_end=(model[end_i,"DrawStart2"])
  model[i,"DrawEnd2"]=draw_end
	
  #Midpoint of junction line
  half_diff=(draw_end-draw_start)/2
  draw_mid=draw_start+half_diff
  model[i,"DrawMid2"]=draw_mid
}

#Set the color of each feature.  Introns=black, Exons=blue, all AE features are red
exon_color="blue"
junction_color="blue"
intron_color="black"
ae_color="red"

model[,"Color"]=intron_color
model[ei,"Color"]=exon_color
if (length(ji) > 0){
 model[ji,"Color"]=junction_color
}
if (length(ae_i) > 0){
 model[ae_i,"Color"]=ae_color
}

#Define a gene label including basic information about the gene (name, exon count, total bases, exon bases)
gene_bases=sum(model[eii,"Size"])
exon_bases=sum(model[ei,"Size"])
gene_label=paste(gene_name, ":  ", exon_count, " exon(s) | ", prettyNum(gene_bases, big.mark=","), " bases | ", prettyNum(exon_bases, big.mark=","), " exon bases", sep="")
bottom_label="Exons are depicted on a log2 scale, introns on a log10 scale"

#Define the line weights that will be used for lines in the plot
lwd0=0.5
lwd1=1
lwd2=2
lwd3=5
lwd4=10

#Define the text cex adjustment
text_cex=1

Cairo(width=800, height=(400), file=outfile, type="svg", pointsize=12, bg="white", units="px", dpi=72)

#Create a blank plot and add some basic labels.
par(mar=c(0.5,0.5,4,0.5))
main_title=paste("Gene model for: ", gene_name, sep="")
plot(x=xrange, y=yrange, xlim=c(0,ceiling(grand_length+(extra_width*2))), ylim=c(0,height), type="n", xlab="", ylab="", xaxt="n", yaxt="n", font.lab=2, main=main_title)
text(x=model[ei[1],"DrawStart1"], y=top_label_scale, labels=gene_label, cex=text_cex, pos=4)
text(x=model[ei[1],"DrawStart1"], y=bottom_label_lower, labels=bottom_label, cex=text_cex, pos=4)

#Draw the end-to-end intron line for the 'to-scale' and main gene models
segments(x0=model[ei[1],"DrawMid1"], y0=center_scale, x1=model[ei[length(ei)],"DrawMid1"],
		  y1=center_scale, col=intron_color, lty=1, lwd=lwd1)
segments(x0=model[ei[1],"DrawMid2"], y0=center_main, x1=model[ei[length(ei)],"DrawMid2"],
		  y1=center_main, col=intron_color, lty=1, lwd=lwd4)

#Then draw thin connecting lines between the exon positions of both gene models
segments(x0=model[ei,"DrawMid1"], y0=bottoms_scale, x1=model[ei,"DrawMid2"],
		  y1=tops_main, col="black", lty=3, lwd=lwd1)

#Now draw the introns individually so that those that are AE will be marked
if (length(ii)){
 segments(x0=model[ii,"DrawStart1"], y0=rep(center_scale, length(ii)), x1=model[ii,"DrawEnd1"],
		   y1=rep(center_scale, length(ii)), col=model[ii,"Color"], lty=1, lwd=lwd1)
 segments(x0=model[ii,"DrawStart2"], y0=rep(center_main, length(ii)), x1=model[ii,"DrawEnd2"],
		   y1=rep(center_main, length(ii)), col=model[ii,"Color"], lty=1, lwd=lwd4)
}

#Now draw the exons for the scale gene model
rect(xleft=model[ei,"DrawStart1"], ybottom=bottoms_scale, xright=model[ei,"DrawEnd1"], ytop=tops_scale, 
	 density=NULL, angle=45, col=model[ei,"Color"], border=model[ei, "Color"], lty=1, lwd=lwd1)


#Now draw the main exons on top of the intron line
rect(xleft=model[ei,"DrawStart2"], ybottom=bottoms_main, xright=model[ei,"DrawEnd2"], ytop=tops_main, 
	 density=NULL, angle=45, col=model[ei,"Color"], border="black", lty=1, lwd=lwd1)

#Draw canonical AE junctions on the bottom
if (length(ae_ji_can) > 0){
 segments(x0=model[ae_ji_can,"DrawStart2"], y0=bottom_main, x1=model[ae_ji_can,"DrawMid2"], y1=bottom_main_j, 
		   col=model[ae_ji_can,"Color"], lty=1, lwd=lwd2)
 segments(x0=model[ae_ji_can,"DrawMid2"], y0=bottom_main_j, x1=model[ae_ji_can,"DrawEnd2"], y1=bottom_main, 
		   col=model[ae_ji_can,"Color"], lty=1, lwd=lwd2)
}
#Draw skip AE junctions on the top
if (length(ae_ji_skip)){
 segments(x0=model[ae_ji_skip,"DrawStart2"], y0=top_main, x1=model[ae_ji_skip,"DrawMid2"], y1=top_main_j, 
	 	   col=model[ae_ji_skip,"Color"], lty=2, lwd=lwd2)
 segments(x0=model[ae_ji_skip,"DrawMid2"], y0=top_main_j, x1=model[ae_ji_skip,"DrawEnd2"], y1=top_main, 
		  col=model[ae_ji_skip,"Color"], lty=2, lwd=lwd2)
}

#Mark the CDS start and end positions 
x_ratio_major=0.1
x_ratio_minor=0.01
y_ratio_minor=0.01
y0=top_label_main-(height*y_ratio_minor)
y1=top_label_main
y2=top_label_main+(height*y_ratio_minor)
arrow_weight=2
if (cds_draw_start){
 x0=cds_draw_start
 x1=x0+(grand_length*(x_ratio_minor*2))
 x2=x0+(grand_length*(x_ratio_minor*3))
 segments(x0=x0, y0=y0, x1=x0, y1=y2, col="black", lty=1, lwd=lwd2) #Vertical line at left
 segments(x0=x0, y0=y1, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Horizontal line
 segments(x0=x1, y0=y0, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Bottom diagonal line
 segments(x0=x1, y0=y2, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Top diagonal line 
}
if (cds_draw_end){
 x0=cds_draw_end-(grand_length*(x_ratio_minor*3)) 
 x1=x0+(grand_length*(x_ratio_minor*2))
 x2=cds_draw_end
 x3=cds_draw_end+((grand_length*(x_ratio_minor/5)*arrow_weight))
 segments(x0=x0, y0=y1, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Horizontal line
 segments(x0=x1, y0=y0, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Bottom diagonal line
 segments(x0=x1, y0=y2, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Top diagonal line 
 segments(x0=x3, y0=y0, x1=x3, y1=y2, col="black", lty=1, lwd=lwd2) #Vertical line at right
}

#Add the labels for each exon region
text(x=model[ei,"DrawMid2"], y=rep(bottom_label_main, exon_count), labels=model[ei,"Name"], cex=text_cex)

#Add a simple legend
legend_pos = "topright"
legend_text = c("Alternative exon usage", "Alternative junction usage", "Alternative intron retention")
legend_cols = c(ae_color, ae_color, ae_color)
legend_ltys = c(NA, 1, 1)
legend_lwds = c(NA, 1, 5)
legend_pchs = c(15, NA, NA)
legend(legend_pos, legend_text, col=legend_cols, lty=legend_ltys, lwd=legend_lwds, pch=legend_pchs, cex=1)

#Close Cairo SVG device
dev.off()





