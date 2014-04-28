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
#Options: 
#[1] file1 (all gene expression values)
#[2] file2 (expression values for all features of a single gene)
#[3] file3 (gene model information)
#[4] gene_name
#[5] ensg_name
#[6] gene_id
#[7] output_dir
#[8] intragenic_cutoffs_string
#[9] intergenic_cutoffs_string
#[10] cds_start
#[11] cds_end

#/home/malachig/svn/alexa_seq/R_bin/geneFigure_SVG.R /projects/malachig/alexa_seq/temp/website/Breast/genes/GeneData_TEMP.txt /projects/malachig/alexa_seq/temp/website/Breast/genes/single_gene_data_TEMP.txt /projects/malachig/alexa_seq/temp/website/Breast/genes/gene_model_TEMP.txt "CLCN6" "ENSG00000011021" 281 /projects/malachig/alexa_seq/temp/website/Breast/genes/chr1_1/images/ "20.13524 15.24836" "9.46 8.87" 114 34074
#infile1 = "/projects/malachig/alexa_seq/temp/website/Breast/genes/GeneData_TEMP.txt"
#infile2 = "/projects/malachig/alexa_seq/temp/website/Breast/genes/single_gene_data_TEMP.txt"
#infile3 = "/projects/malachig/alexa_seq/temp/website/Breast/genes/gene_model_TEMP.txt"
#gene_name = "FRMD4B" 
#ensg_name = "ENSG00000114541" 
#gene_id = as.numeric(4381)
#results_dir = "/projects/malachig/alexa_seq/temp/website/Breast/genes/chr3_3/images/"
#intragenic_cutoffs_string = "16.81684 15.50229" 
#intergenic_cutoffs_string = "9.46 8.87" 
#cds_start = as.numeric(72763)
#cds_end = as.numeric(214419)

#Get user specified arguments
args = (commandArgs(TRUE))
infile1 = args[1];
infile2 = args[2];
infile3 = args[3];
gene_name = args[4];
ensg_name = args[5]
gene_id = as.numeric(args[6]);
results_dir = args[7];
intragenic_cutoffs_string = args[8];
intergenic_cutoffs_string = args[9];
cds_start = as.numeric(args[10])
cds_end = as.numeric(args[11])

#Define some hardcoded colors
bg = "white"
color1="#B3CDE3" #Exons, etc. in gene model
color2="#F03B20" #AE features in gene model 
color3="#C8E4B4" #Highlighting of DE events in line plot
color4="#FFFF99" #Highlighting of AE events in line plot
color5="#FF00FF" #Asterixes in bar plots 

all_gene_data = read.table(infile1, header=TRUE, sep = "\t", na.strings = "NA")
single_gene_data = read.table(infile2, header=TRUE, sep = "\t", na.strings = "NA", as.is=c(2,3))

#Data needed to create simple gene model visualization
#Start   End     Name    Type    IsDE   IsAE    SkippedExons
model=read.table(infile3, header=TRUE, sep="\t", na.strings=c("NA"), as.is=c(3,4))

intragenic_cutoffs = as.numeric(strsplit(intragenic_cutoffs_string, " ")[[1]])
intergenic_cutoffs = as.numeric(strsplit(intergenic_cutoffs_string, " ")[[1]])

#Create a plot to accomodate the number of libraries
lib_count = length(all_gene_data[1,])-1
lib_names = names(all_gene_data[2:length(all_gene_data[1,])])
feature_count = length(single_gene_data[,1])

#Create a vector of colors where colors are assinged for each feature based on their feature type
#Types and color
#Gene(1) Transcript(2) ExonRegion(3) KnownJunction(4) NovelJunction(4) KnownBoundary(5) NovelBoundary(5) Intron(6) ActiveIntronRegion(6) SilentIntronRegion(6) Intergenic(7) ActiveIntergenicRegion(7) SilentIntergenicRegion(7)
rain = rainbow(7)
cols = rep("black", feature_count)

cols[which(single_gene_data[,"SeqType"]=="Gene")] = rain[1]
cols[which(single_gene_data[,"SeqType"]=="Transcript")] = rain[2]
cols[which(single_gene_data[,"SeqType"]=="ExonRegion")] = "forestgreen"
cols[which(single_gene_data[,"SeqType"]=="KnownJunction")] = rain[6]
cols[which(single_gene_data[,"SeqType"]=="NovelJunction")] = rain[5]
cols[which(single_gene_data[,"SeqType"]=="KnownBoundary")] = rain[3]
cols[which(single_gene_data[,"SeqType"]=="NovelBoundary")] = rain[4]
cols[which(single_gene_data[,"SeqType"]=="Intron")] = "brown"
cols[which(single_gene_data[,"SeqType"]=="ActiveIntronRegion")] = "brown"
cols[which(single_gene_data[,"SeqType"]=="SilentIntronRegion")] = "brown"
cols[which(single_gene_data[,"SeqType"]=="Intergenic")] = "grey90"
cols[which(single_gene_data[,"SeqType"]=="ActiveIntergenicRegion")] = "grey90"
cols[which(single_gene_data[,"SeqType"]=="SilentIntergenicRegion")] = "grey90"

legend1_text = c(gene_name, "Intragenic Cutoff", "Intergenic Cutoff")
legend1_cols = c("red", "black", "black")
legend1_pos = "topright"

legend2_text = c("Gene", "Transcript", "Exon Region", "Known Junction", "Novel Junction", "Known Boundary", "Novel Boundary", "Intronic", "Intergenic", "Significant AE")
legend2_cols = c(rain[1], rain[2], "forestgreen", rain[6], rain[5], rain[3], rain[4], "brown", "grey90", "#FF00FF")
legend2_pch = c(15,15,15,15,15,15,15,15,15,8)
legend2_pos = "bottomright"

#Get the maximum value of Y for all features for all libraries
max_y = 0;
for(i in 1:lib_count){
  data = log2(single_gene_data[,i+5]+1)
  max_data = max(data, na.rm = TRUE)
  if (max_data > max_y){
    max_y = max_data
  }
}

#Get the number of features so that the font size for feature names can be altered accordingly
#cex.names of 1.5 works well for 25 features
#If the feature count is greater than limit, decrease this by 'vary_rate' for every 1 additional features above the limit
#Do not allow cex.names to be less than min_size
feature_count = (length(single_gene_data[,1]))
min_size = 0.5
names_size = 1.5
limit = 25
vary_rate = 0.01
if (feature_count > limit){
  excess = feature_count - limit
  value = excess*vary_rate
  names_size = round(names_size-value, digits=2)
}
if(names_size < min_size){
  names_size = min_size
}

#Define the figure height and width
image_width=1000
image_height=(400*((lib_count*2)+1))

#Using SVG functionality of the SVG package - OPEN SVG DEVICE
library(Cairo)
outfile = paste(results_dir, ensg_name, ".svg", sep="")
Cairo(width=image_width, height=image_height, file=outfile, type="svg", pointsize=12, bg=bg, units="px", dpi=72)

#Define a subset of all features consisting only of the exon regions and known or novel exon junctions
selection = which((single_gene_data[,"SeqType"]=="ExonRegion") | (single_gene_data[,"SeqType"]=="KnownJunction") | (single_gene_data[,"SeqType"]=="NovelJunction"))

#Define a subset of all features consisting only of the exon regions and known or novel exon junctions that were ALSO alternatively expressed
temp_x = single_gene_data[selection,]
selection_de = which(temp_x[,"IsDE"]==1)
selection_ae = which(temp_x[,"IsAE"]==1)

#Set the points character according to these three types:
pch_list = rep(0, length(selection))
temp_data = single_gene_data[selection,]
pch_list[which(temp_data[,"SeqType"]=="ExonRegion")] = 15;
pch_list[which(temp_data[,"SeqType"]=="KnownJunction")] = 16;
pch_list[which(temp_data[,"SeqType"]=="NovelJunction")] = 17;

#Set up plotting area
rows_needed = (lib_count*2)
if (length(selection) > 0){
  rows_needed=rows_needed+1
}
rows_needed=rows_needed+1
par(bg=bg, font.main=2, font.axis=2, font.lab=2, mfrow=c(rows_needed,1), mar=(c(10, 4, 4, 2)+0.1))



###############################################################################################################################################
#A.) Gene model drawing                                                                                                                       #
###############################################################################################################################################

#Eliminate any features that would be plotted outside the boundaries of the first and last exon (should not really happen normally)
ei = which(model[,"Type"]=="ExonRegion")
ref_start = model[ei[1],"Start"]
ref_end = model[ei[length(ei)],"End"]
oki = which(model[,"Start"] >= ref_start & model[,"End"] <= ref_end)
model = model[oki,]

#Determine size of each feature
model[,"Size"]=(model[,"End"]-model[,"Start"])+1

#Get indexes of Exons (E), Introns (I), Junctions (J), Boundaries (B), and ActiveIntronRegions (A)
ei = which(model[,"Type"]=="ExonRegion")
ii = which(model[,"Type"]=="Intron")
ji = which(model[,"Type"]=="Junction")
bi = which(model[,"Type"]=="Boundary")
ai = which(model[,"Type"]=="ActiveIntronRegion")
eii = which(model[,"Type"]=="ExonRegion" | model[,"Type"]=="Intron")
ae_i = which(model[,"IsAE"]==1)
ae_ji = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1)
ae_ji_can = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1 & model[,"SkippedExons"]==0)
ae_ji_skip = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1 & model[,"SkippedExons"]>0)
ae_bi = which(model[,"Type"]=="Boundary" & model[,"IsAE"]==1)
ae_ai = which(model[,"Type"]=="ActiveIntronRegion" & model[,"IsAE"]==1)

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
if (length(bi) > 0){
 model[bi,"DrawSize2"]=model[bi,"Log2Size"]
}

#What is the total length needed?
grand_length=sum(model[eii,"DrawSize2"])

#Add some padding to the intron sizes, so that very small introns will still be visible
intron_padding_ratio = 0.005
intron_padding = grand_length*intron_padding_ratio
if (length(ii) > 0){
 model[ii,"DrawSize2"] = model[ii,"DrawSize2"]+intron_padding
}

#Now recalculate the grand size after adding the padding
grand_length=sum(model[eii,"DrawSize2"])


#Get draw size1 relative to draw size2
grand_length_scale=sum(model[eii,"Size"])
model[,"DrawSize1"]=(model[,"Size"]/grand_length_scale)*grand_length

#Define some extra width for the plot
extra_width_ratio=0.01
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
bottom_label_main=(center_main)-(height_main*3)
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
  start_i=which((model[,"End"]==j_start) & (model[,"Type"]=="ExonRegion"))
  if (length(start_i) == 0){
    start_i=which((model[,"End"]==(j_start-1)) & (model[,"Type"]=="ExonRegion"))
  }
  if (length(start_i) == 0){
    start_i=which((model[,"End"]==(j_start+1)) & (model[,"Type"]=="ExonRegion"))
  }
  draw_start=(model[start_i,"DrawEnd2"])
  model[i,"DrawStart2"]=draw_start

  #End of junction line
  j_end=(model[i,"End"])
  end_i=which((model[,"Start"]==j_end) & (model[,"Type"]=="ExonRegion"))
  if(length(end_i) == 0){
    end_i=which((model[,"Start"]==(j_end-1)) & (model[,"Type"]=="ExonRegion"))
  }
  if(length(end_i) == 0){
    end_i=which((model[,"Start"]==(j_end+1)) & (model[,"Type"]=="ExonRegion"))
  }

  draw_end=(model[end_i,"DrawStart2"])
  model[i,"DrawEnd2"]=draw_end
	
  #Midpoint of junction line
  half_diff=(draw_end-draw_start)/2
  draw_mid=draw_start+half_diff
  model[i,"DrawMid2"]=draw_mid
}

#Get the drawing x positions for only those boundaries that were AE
#The x position will be the centre of the boundary sequence (this should correspond to the edge of an exon)
for(i in ae_bi){

  #Midpoint of boundary seq
  b_mid_l = model[i,"Start"] + (((model[i,"End"]-model[i,"Start"])+1)/2)-1
  b_mid_r = model[i,"Start"] + (((model[i,"End"]-model[i,"Start"])+1)/2)

  #Find mid point of boundary to draw (find 5' or 3' end of an exon that matches the centre of the boundary seq)
  #Exon acceptors
  coords=c(b_mid_l, b_mid_r, b_mid_l-1, b_mid_l+1, b_mid_r-1, b_mid_r+1)
  side=NULL
  mid_i=NULL
  for (j in coords){
    if(length(mid_i) == 0){
      mid_i=which((model[,"Start"]==j) & (model[,"Type"]=="ExonRegion"))
      side=1
    }
    if(length(mid_i) == 0){
      mid_i=which((model[,"End"]==j) & (model[,"Type"]=="ExonRegion"))
      side=2
    }
  }
  #If more than one exon match was found, chose one
  if (length(mid_i) > 1){
    mid_i=mid_i[1]
  }
  draw_mid1=NULL
  draw_mid2=NULL
  if (side == 1){
    draw_mid1=(model[mid_i,"DrawStart1"])
    draw_mid2=(model[mid_i,"DrawStart2"])
  }else{
    draw_mid1=(model[mid_i,"DrawEnd1"])
    draw_mid2=(model[mid_i,"DrawEnd2"])
  }
  model[i,"DrawMid1"]=draw_mid1
  model[i,"DrawMid2"]=draw_mid2
}

#Get the drawing x positions for only those 'active intron regions' were AE
for(i in ae_ai){

  #Get the intron that contains the active intron region
  ai_start = model[i,"Start"]
  ai_end = model[i,"End"]
  ai_center = ceiling(ai_start+(((ai_end-ai_start)+1)/2))

  ai_intron = which((ai_center >= model[,"Start"]) & (ai_center <= model[,"End"]) & model[,"Type"]=="Intron")

  draw_mid1=NULL
  draw_mid2=NULL
  if (length(ai_intron) > 0){
    if(length(ai_intron) > 1){
      ai_intron=ai_intron[1]
    }
    relative_pos=(model[ai_intron,"Size"] - (model[ai_intron,"End"]-ai_center))/model[ai_intron,"Size"]
    draw_mid1=model[ai_intron,"DrawStart1"]+model[ai_intron,"DrawSize1"]*relative_pos
    draw_mid2=model[ai_intron,"DrawStart2"]+model[ai_intron,"DrawSize2"]*relative_pos
  }
  model[i,"DrawMid1"]=draw_mid1
  model[i,"DrawMid2"]=draw_mid2
}

#Set the color of each feature.
exon_color=color1
junction_color=color1
boundary_color=color1
intron_color="black"
ae_color=color2

model[,"Color"]=intron_color
model[ei,"Color"]=exon_color
if (length(ji) > 0){
 model[ji,"Color"]=junction_color
}
if (length(bi) > 0){
 model[bi,"Color"]=boundary_color
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
model_text_cex=1

#Create a blank plot and add some basic labels.
#Reduce the margin sizes for this plot
par(mar=c(0.5,0.5,4,0.5))
main_title=paste("Gene model for '", gene_name, "'", sep="")
plot(x=xrange, y=yrange, xlim=c(0,ceiling(grand_length+(extra_width*2))), ylim=c(0,height), type="n", xlab="", ylab="", xaxt="n", yaxt="n", main=main_title, cex.main=2)
text(x=model[ei[1],"DrawStart1"], y=top_label_scale, labels=gene_label, cex=model_text_cex, pos=4, font=2)
text(x=model[ei[1],"DrawStart1"], y=bottom_label_lower, labels=bottom_label, cex=model_text_cex, pos=4, font=2)

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
	 density=NULL, angle=45, col=model[ei,"Color"], border="black", lty=1, lwd=lwd0)


#Now draw the main exons on top of the intron line
rect(xleft=model[ei,"DrawStart2"], ybottom=bottoms_main, xright=model[ei,"DrawEnd2"], ytop=tops_main, 
	 density=NULL, angle=45, col=model[ei,"Color"], border="black", lty=1, lwd=lwd1)

#Now mark AE boundaries at edges of exons
if (length(ae_bi) > 0){
  points(x=model[ae_bi,"DrawMid1"], y=rep(center_scale, length(ae_bi)), col=ae_color, pch=8)
  points(x=model[ae_bi,"DrawMid2"], y=rep(center_main, length(ae_bi)), col=ae_color, pch=8)
}

#Now draw AE active intron regions within introns
if (length(ae_ai) > 0){
  points(x=model[ae_ai,"DrawMid1"], y=rep(center_scale, length(ae_ai)), col=ae_color, pch=8)
  points(x=model[ae_ai,"DrawMid2"], y=rep(center_main, length(ae_ai)), col=ae_color, pch=8)
}

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
y3=top_label_main+((height*y_ratio_minor)*4)
arrow_weight=2
if (cds_draw_start){
 x0=cds_draw_start
 x1=x0+(grand_length*(x_ratio_minor*2))
 x2=x0+(grand_length*(x_ratio_minor*3))
 segments(x0=x0, y0=y0, x1=x0, y1=y2, col="black", lty=1, lwd=lwd2) #Vertical line at left
 segments(x0=x0, y0=y1, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Horizontal line
 segments(x0=x1, y0=y0, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Bottom diagonal line
 segments(x0=x1, y0=y2, x1=x2, y1=y1, col="black", lty=1, lwd=lwd2) #Top diagonal line 
 text(x=x0, y=y3, labels="CDS Start", cex=model_text_cex, font=2)
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
 text(x=x0, y=y3, labels="CDS End", cex=model_text_cex, font=2)
}

#Determine the subset of exon labels that should be used (as the number of exons grows, only display a subset of labels)
#Adjust the value of 'modulo' according to the number of exons (for each multiple of exon_count_breaks, start skipping more label)
exon_count_breaks=25
modulo=ceiling(exon_count/exon_count_breaks)
labels_i = c(1:length(ei))
display_labels_i=which(labels_i%%modulo == 0 | labels_i == 1)


#Add the labels for each exon region
text(x=model[ei[display_labels_i],"DrawMid2"], y=rep(bottom_label_main, length(display_labels_i)), labels=model[ei[display_labels_i],"Name"], cex=model_text_cex, font=2)

#Add a simple legend
legend_pos = "topright"
legend_text = c("Alternative exon usage", "Alternative junction usage", "Alternative intron retention", "Alternative boundary/cryptic exon")
legend_cols = c(ae_color, ae_color, ae_color, ae_color)
legend_ltys = c(NA, 1, 1, NA)
legend_lwds = c(NA, 1, 5, NA)
legend_pchs = c(15, NA, NA, 8)
legend(legend_pos, legend_text, col=legend_cols, lty=legend_ltys, lwd=legend_lwds, pch=legend_pchs, cex=model_text_cex)

#Reset margins to defaults
par(mar=c(10,4,4,2)+0.1)

###############################################################################################################################################
#B.) Line graph - One total                                                                                                                   #
###############################################################################################################################################
#    Line graph showing expression values of exons, and known+novel exon junctions for all libraries on one plot
#    To do this, create a blank barplot first with the first library.  Set the the border to 'NA' and the color to "white" to get a blank barplot with the desired labels
#    When generating the barplot, use a variable to capture the barplot midpoints (these will be your x values in the following lines() commands)
#    Use lines() to plot the expression levels for each library on this plot
lib_cols = c("blue", "red")
if (lib_count == 3){
  lib_cols = c("red", "blue", "magenta")
}
if (lib_count == 4){
  lib_cols = c("blue", "blue", "red", "red")
}
if (lib_count == 6){
  lib_cols = c("blue", "blue", "red", "red", "magenta", "magenta")
}
if (lib_count == 5 | lib_count > 6){
  lib_cols = rainbow(lib_count)
}

lib_ltys = rep(c(1,2),1000)
legend0_text = c(lib_names[1])
legend0_cols = lib_cols[1]
legend0_ltys = lib_ltys[1]
legend0_pos = "topright"

for(i in 2:lib_count){
  legend0_text = c(legend0_text, lib_names[i])
  legend0_cols = c(legend0_cols, lib_cols[i])
  legend0_ltys = c(legend0_ltys, lib_ltys[i])
}
if (length(selection) > 0){
  data = log2(single_gene_data[selection,6]+1)

  #Add some dummy values to the data to make room for the legend - try increasing number of entries by 15%
  extra_percent1 = 0.15
  extra_values=round(length(data)*extra_percent1, digits=0)
  if(extra_values<1){
    extra_values=1
  }
  dummy_data=c(rep(max_y, length(data)), rep(0, extra_values))

  #Set the colors for dummy bars (then define those that correspond to AE events to provide highlighting behind the line plots)
  bar_colors=(rep("white", length(dummy_data)))
  bar_colors[selection_de]=color3;
  bar_colors[selection_ae]=color4;

  #Set up the labels
  main_title = "Exon and junction expression levels (all libraries)"
  y_label = "log2 (expression level +1)"
  feature_names=single_gene_data[selection,"SeqName"]
  new_feature_names=c(feature_names, rep("", extra_values))

  #Create a blank barplot
  z = barplot(dummy_data, col=bar_colors, border=NA, names.arg=new_feature_names, main=main_title, ylab=y_label, ylim=c(0,max_y),
             las=2, col.lab = gray(.1), cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, cex.names=names_size)

  #Add a 0 line and tick marks to the base line
  abline(h=0, col="black", lty=1, lwd=1)
  points(z[1:length(data)], rep(0, length(selection)), pch=3, col="black", lwd=1) 

  #Plot each line of data for the selected features
  for(i in 1:lib_count){
    data = log2(single_gene_data[selection,i+5]+1)
    lines(z[1:length(data)], data, col=lib_cols[i], lwd=2, lty=lib_ltys[i])
    points(z[1:length(data)], data, col="black", lwd=3, pch=pch_list)
  }
  legend(legend0_pos, legend0_text, col=legend0_cols, lty=legend0_ltys, lwd=2, cex=1.5)
}

#Get the overall max and min x values for all genes all libraries
max_x = 0
min_x = 10000000000000000000;
for(i in 1:lib_count){
  x = log2(all_gene_data[,i+1])
  y = max(x, na.rm=TRUE)
  z = min(x, na.rm=TRUE)
  if (y > max_x){
    max_x=y
  }
  if (z < min_x){
    min_x=z
  }
}


###############################################################################################################################################
#C.) Histograms - One per library                                                                                                             #
###############################################################################################################################################
#Histograms
target_gene = which(all_gene_data[,1]==gene_id)
for(i in 1:lib_count){
  
  #Histogram of gene expression data showing position of current gene
  main_title = lib_names[i]
  x_label = "log2 gene expression level - all genes"
  y_label = "frequency"
  hist(x = log2(all_gene_data[,i+1]), col="blue", xlim=c((min_x-1),(max_x+1)), main=main_title, xlab=x_label, ylab=y_label,
       col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.5, cex.axis=1.5, breaks = 100)
  target_gene_level = log2(all_gene_data[target_gene,i+1])
  abline(v=target_gene_level, col="red", lwd=2, lty=2)
  abline(v=log2(intragenic_cutoffs[i]), col="black", lwd=2, lty=2)
  abline(v=log2(intergenic_cutoffs[i]), col="black", lwd=2, lty=3)
  legend(legend1_pos, legend1_text, col=legend1_cols, lty=c(2,2,3), lwd=2, cex=1.5)
}


###############################################################################################################################################
#D.) Bar Plots - One per library                                                                                                              #
###############################################################################################################################################
#     Barplot showing expression level of all individual features
#     Show legend and feature names for first library but not subsequent ones

main_title = paste(gene_name, " feature expression levels for library: ", lib_names[1], " (", feature_count," expressed features)", sep="")
y_label = "log2 (expression level +1)"

#Get data.  Convert values to log2, unless they are 0 in which case set them to 0
data = log2(single_gene_data[,6]+1)

#Get the subset of features that were significantly AE
ae_i = which(single_gene_data[,"IsAE"]==1)


#Add some dummy values to the data to make room for the legend - try increasing number of entries by 15%
extra_percent2=0.17
extra_values=round(length(data)*extra_percent2, digits=0)
if(extra_values<1){
  extra_values=1
}
new_data=c(data, rep(0, extra_values))
feature_names=single_gene_data[,"SeqName"]
new_feature_names=c(feature_names, rep("", extra_values))

#Plot the first library
z=barplot(new_data, col=cols, names.arg=new_feature_names, main=main_title, ylab=y_label, ylim=c(0,max_y),
        las=2, col.lab = gray(.1), cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, cex.names=names_size)
abline(h=log2(intragenic_cutoffs[1]+1), col="black", lwd=2, lty=2)
abline(h=log2(intergenic_cutoffs[1]+1), col="black", lwd=2, lty=3)
points(x=z[ae_i], y=new_data[ae_i], col="#FF00FF", pch=8, cex=1.5) 
legend(legend2_pos, legend2_text, col=legend2_cols, pch=legend2_pch, cex=1.5, pt.cex=1.4)

#Now for all other libraries
for(i in 2:lib_count){
  main_title = paste(gene_name, " feature expression levels for library: ", lib_names[i], sep="")
  y_label = "log2 (expression level +1)"

  #Get data.  Convert values to log2, unless they are 0 in which case set them to 0
  data = log2(single_gene_data[,i+5]+1)
  new_data=c(data, rep(0, extra_values))

  z=barplot(new_data, col=cols, main=main_title, ylab=y_label, names.arg=new_feature_names, ylim=c(0,max_y),
          las=2, col.lab = gray(.1), cex.main = 1.5, cex.lab=1.5, cex.axis=1.5, cex.names=names_size)
  abline(h=log2(intragenic_cutoffs[i]+1), col="black", lwd=2, lty=2)
  abline(h=log2(intergenic_cutoffs[i]+1), col="black", lwd=2, lty=3)
  points(x=z[ae_i], y=new_data[ae_i], col=color5, pch=8, cex=1.5) 
  legend(legend2_pos, legend2_text, col=legend2_cols, pch=legend2_pch, cex=1.5, pt.cex=1.4)
}

#CLOSE SVG DEVICE
x = dev.off()

