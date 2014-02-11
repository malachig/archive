#############################################################################################################################################################
#Draw model function                                                                                                                                        #
#############################################################################################################################################################
#Draw a gene model for the region immediately encompassing the splice mutation
#Use a genome wide object containing exon content blocks (i.e. merged exons)
#Draw exons and introns for this region, label exons with donor/acceptor sites
#Mark the mutation position
#Mark the canonical and alternative junctions and label them
drawModel = function(ec_data, jc_data, mut_position, chr, strand, grand_start, grand_end, j_list, j_names, plot_name, kj_color, aj_color, intron_color, exon_color, mutant_color, mut_lib_name_list, active){

#Get the exon content for the specified region of interest
eci=which((ec_data[,"chr"] == chr) & (ec_data[,"strand"] == strand) & ((ec_data[,"end"] >= grand_start-100) & (ec_data[,"start"] <= grand_end+100)))
model=ec_data[eci,c("start","end")]
names(model) = c("Start","End")
ref_start = model[1,"Start"]
ref_end = model[length(model[,1]),"End"]
oki = which(model[,"Start"] >= ref_start & model[,"End"] <= ref_end)
model = model[oki,]
model[, "Type"]="Exon"
model[, "Name"]=NA
model[, "JID"] = NA
model[, "IsAE"]=0
rownames(model) = 1:length(model[,1])
#reorder
model=model[order(model[,"Start"]),]
rownames(model) = 1:length(model[,1])

#Get the intron sizes (if applicable)
n = length(model[,1])
if (n > 1){
  for(j in 1:(n-1)){
    m = length(model[,1])
    intron_start=model[j,"End"]
    intron_end=model[j+1,"Start"]
    intron_size=intron_end-intron_start
    if(intron_size > 1){
      model[m+1, "Start"] = intron_start+1
      model[m+1, "End"] = intron_end-1
      model[m+1, "Type"] = "Intron"
      model[m+1, "Name"] = NA
      model[m+1, "JID"] = NA
      model[m+1, "IsAE"] = 0
    }
  }
  #reorder
  model=model[order(model[,"Start"]),]
  rownames(model) = 1:length(model[,1])
}

#Add in the mutation
n = length(model[,1])
model[n+1, "Start"] = mut_pos
model[n+1, "End"] = mut_pos
model[n+1, "Type"] = "Mutant"
model[n+1, "Name"] = NA
model[n+1, "JID"] = NA
model[n+1, "IsAE"] = 0
#reorder
model=model[order(model[,"Start"]),]
rownames(model) = 1:length(model[,1])

#Get the junction coords and add to the model
tmp=jc_data[j_list,]
tmp_names=names(j_list)
for(j in 1:length(j_list)){
  m = length(model[,1])
  model[m+1, "Start"] = tmp[j,"start"]
  model[m+1, "End"] = tmp[j,"end"]
  model[m+1, "Type"] = "Junction"
  model[m+1, "Name"] = tmp_names[j]
  model[m+1, "JID"] = tmp[j,"jid"]
  if (tmp[j,"type"] == "Alternate"){
    model[m+1, "IsAE"] = 1
  }else{
    model[m+1, "IsAE"] = 0
  }
}
model=model[order(model[,"Start"]),]
rownames(model) = 1:length(model[,1])

#Determine size of each feature
model[,"Size"]=(model[,"End"]-model[,"Start"])+1

#Get indexes of Exons (E), Introns (I), Junctions (J)
mi = which(model[,"Type"]=="Mutant")
ei = which(model[,"Type"]=="Exon")
ii = which(model[,"Type"]=="Intron")
ji = which(model[,"Type"]=="Junction")
kji = which(model[,"Type"]=="Junction" & model[,"IsAE"]==0)
aji = which(model[,"Type"]=="Junction" & model[,"IsAE"]==1)
eii = which(model[,"Type"]=="Exon" | model[,"Type"]=="Intron")

#Draw a single gene models, on a log2/log10 scale in the center of the plot
#Create scaled coordinates, modify sizes as follows:
#exons on log2 scale
#introns on log10 scale
model[,"Log2Size"]=log2(model[,"Size"])+1
model[,"Log10Size"]=log10(model[,"Size"])+1
model[,"DrawSize"]=NA
model[ei,"DrawSize"]=model[ei,"Log2Size"]
model[mi,"DrawSize"]=model[mi,"Log2Size"]
if (length(ii) > 0){
 model[ii,"DrawSize"]=model[ii,"Log10Size"]
}
if (length(ji) > 0){
 model[ji,"DrawSize"]=model[ji,"Log2Size"]
}

#What is the total length needed?
grand_length=sum(model[eii,"DrawSize"])

#Add some padding to the intron sizes, so that very small introns will still be visible
intron_padding_ratio = 0.005
intron_padding = grand_length*intron_padding_ratio
if (length(ii) > 0){
 model[ii,"DrawSize"] = model[ii,"DrawSize"]+intron_padding
}

#Now recalculate the grand size after adding the padding
grand_length=sum(model[eii,"DrawSize"])

#Define some extra width for the plot
extra_width_ratio=0.01
extra_width=ceiling(grand_length*extra_width_ratio)

#Define some dummy values to help create a blank plot of the correct size
xrange=0:ceiling(grand_length+(extra_width*2))
yrange=rep(5,length(xrange))

#Width is the cumulative size of all exons/introns plus some extra width
#Height is a constant defined below
height=10
y_inset=2 
center_main=0+(height/y_inset)

#How many features (exons+introns) are to be plotted
feature_count=length(eii)
exon_count=length(ei)

#Define the y positions to be used for various parts of the plot
#All y positions are in an arbitrary space of range 0-10
#'main' refers to the y positions of the main gene model
#'j' refers to the y positions of junction connecting lines
#'label' refers to the y positions of text labels
height_ratio_main=0.05

height_main=(height*height_ratio_main)
bottom_main=(center_main)-height_main
top_main=(center_main)+height_main
bottom_main_j=(center_main)-(height_main*2)
top_main_j=(center_main)+(height_main*2)
bottom_label_main=(center_main)-(height_main*3)
bottom_label_mid=(center_main)-(height_main*4)
bottom_label_lower=(center_main)-(height_main*6)
top_label_main=(center_main)+(height_main*3)
top_label_mid=(center_main)+(height_main*4)
top_label_upper=(center_main)+(height_main*6)


#Define the y coordinates of the gene model to be drawn
bottoms_main=rep(bottom_main, exon_count)
tops_main=rep(top_main, exon_count)

#Define the x coordinates of the gene model to be drawn
#i.e. get the mixed log2/log10 scale drawing x coordinates
model[,"DrawStart"]=NA
model[eii,"DrawStart"]=(cumsum(model[eii,"DrawSize"]))-(model[eii,"DrawSize"])
model[eii,"DrawStart"]=model[eii,"DrawStart"]+extra_width
model[,"DrawEnd"]=NA
model[eii,"DrawEnd"]=cumsum(model[eii,"DrawSize"])
model[eii,"DrawEnd"]=model[eii,"DrawEnd"]+extra_width
model[,"DrawMid"]=NA
model[eii,"DrawMid"]=model[eii,"DrawStart"]+((model[eii,"DrawEnd"]-model[eii,"DrawStart"])/2)

#Get the drawing x positions for junctions
#Remember that the 'start' of a junction is the end of the first exon, and 'end' is the start of the second exon
#Also remember that we are dealing with exon content blocks.  So the actual junction may be within the exon content...
#Figure out which exon each side of the junction falls within.  Then figure out the relative position within the exon
#Then use the relative position to determine the relative drawing start, end and midpoint
for(i in ji){
  #Start of junction line
  j_start=(model[i,"Start"])
  start_i=which((model[,"Start"]<=j_start) & (model[,"End"]>=j_start) & (model[,"Type"]=="Exon" | model[,"Type"]=="Intron"))
  start_r=(j_start - model[start_i[1],"Start"])/(model[start_i[1],"Size"])
  draw_start=(model[start_i[1],"DrawSize"]*start_r) + (model[start_i[1],"DrawStart"])
  model[i,"DrawStart"]=draw_start

  #End of junction line
  j_end=(model[i,"End"])
  end_i=which((model[,"Start"]<=j_end) & (model[,"End"]>=j_end) & (model[,"Type"]=="Exon" | model[,"Type"]=="Intron"))
  end_r=(j_end - model[end_i[1],"Start"])/(model[end_i[1],"Size"])
  draw_end=(model[end_i[1],"DrawSize"]*end_r) + (model[end_i[1],"DrawStart"])
  model[i,"DrawEnd"]=draw_end
	
  #Midpoint of junction line
  half_diff=(draw_end-draw_start)/2
  draw_mid=draw_start+half_diff
  model[i,"DrawMid"]=draw_mid
}

#Get the drawing x positions for the mutation(s) - the mutation could be within an intron or an exon
#mutation is a single point - so start and end are the same
for(i in mi){
  m_start=(model[i,"Start"])
  start_i=which((model[,"Start"]<=m_start) & (model[,"End"]>=m_start) & (model[,"Type"]=="Exon" | model[,"Type"]=="Intron"))
  start_r=(m_start - model[start_i[1],"Start"])/(model[start_i[1],"Size"])
  draw_start=(model[start_i[1],"DrawSize"]*start_r) + (model[start_i[1],"DrawStart"])
  model[i,"DrawStart"]=draw_start
  model[i,"DrawEnd"]=draw_start
  model[i,"DrawMid"]=draw_start
}

#Set the color of each feature.
model[,"Color"]=intron_color
model[ei,"Color"]=exon_color
model[mi,"Color"]=mutant_color
if (length(ji) > 0){model[ji,"Color"]=kj_color}
if (length(aji) >0){model[aji,"Color"]=aj_color}

#Define a gene label including basic information about the gene (name, exon count, total bases, exon bases)
gene_label=plot_name
bottom_label="Exons are depicted on a log2 scale, introns on a log10 scale"

#Define the line weights that will be used for lines in the plot
lwd0=0.5
lwd1=1
lwd2=2
lwd3=3
lwd4=10

#Define the text cex adjustment
model_text_cex=1

#Create a blank plot and add some basic labels.
#Reduce the margin sizes for this plot
par(mar=c(0.5,0.5,4,0.5))
main_title=paste("Splicing model for '", gene_label, "'", sep="")

plot(x=xrange, y=yrange, xlim=c(0,ceiling(grand_length+(extra_width*2))), ylim=c(0,height), type="n", xlab="", ylab="", xaxt="n", yaxt="n", main=main_title, cex.main=1.25)
text(x=model[ei[1],"DrawStart"], y=bottom_label_lower, labels=bottom_label, cex=1.15, pos=4, font=2)

#Draw the end-to-end intron line for the gene model
segments(x0=model[ei[1],"DrawMid"], y0=center_main, x1=model[ei[length(ei)],"DrawMid"],
		  y1=center_main, col=intron_color, lty=1, lwd=lwd4)


#Now draw the main exons on top of the intron line
rect(xleft=model[ei,"DrawStart"], ybottom=bottoms_main, xright=model[ei,"DrawEnd"], ytop=tops_main, 
	 density=NULL, angle=45, col=model[ei,"Color"], border="black", lty=1, lwd=lwd1)


#Now mark the mutation
if (length(mi) > 0){
  points(x=model[mi,"DrawMid"], y=rep(center_main, length(mi)), col=model[mi,"Color"], pch=8, cex=4)
}

#Draw canonical AE junctions on the bottom and label them (short name AND full junction ID)
if (length(kji) > 0){
  segments(x0=model[kji,"DrawStart"], y0=bottom_main, x1=model[kji,"DrawMid"], y1=bottom_main_j, 
		   col=model[kji,"Color"], lty=1, lwd=lwd3)
  segments(x0=model[kji,"DrawMid"], y0=bottom_main_j, x1=model[kji,"DrawEnd"], y1=bottom_main, 
		   col=model[kji,"Color"], lty=1, lwd=lwd3)
  text(x=model[kji,"DrawMid"], y=rep(bottom_label_main, length(kji)), labels=model[kji,"Name"], cex=model_text_cex, font=2)
  #text(x=model[kji,"DrawMid"], y=rep(bottom_label_mid, length(kji)), labels=model[kji,"JID"], cex=model_text_cex, font=2)
}
#Draw alternate junctions on the top and label them (short name AND full junction ID)
if (length(aji)){
  segments(x0=model[aji,"DrawStart"], y0=top_main, x1=model[aji,"DrawMid"], y1=top_main_j, 
	 	   col=model[aji,"Color"], lty=2, lwd=lwd3)
  segments(x0=model[aji,"DrawMid"], y0=top_main_j, x1=model[aji,"DrawEnd"], y1=top_main, 
		  col=model[aji,"Color"], lty=2, lwd=lwd3)
  text(x=model[aji,"DrawMid"], y=rep(top_label_main, length(aji)), labels=model[aji,"Name"], cex=model_text_cex, font=2)
  #text(x=model[aji,"DrawMid"], y=rep(top_label_mid, length(aji)), labels=model[aji,"JID"], cex=model_text_cex, font=2)
}

#Add SVG mouse-over and clickable box for each exon
if (active == 1){
  for(i in ei){
    url_text=paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=", chr, ":", model[i,"Start"], "-", model[i,"End"], sep="")
    title_text = "Exon"
    desc1_text = paste(chr, ":", model[i,"Start"], "-", model[i,"End"], sep="")
    desc2_text = paste("Exon Size: ", model[i,"Size"], " bp", sep="")
    setSVGShapeToolTip(title=title_text, desc1=desc1_text, desc2=desc2_text)
    setSVGShapeURL(url_text, target="_blank")  
    rect(xleft=model[i,"DrawStart"], ybottom=bottom_main, xright=model[i,"DrawEnd"], ytop=top_main, col=NULL, density=NULL, border=NA)
  }
}

#Add SVG mouse-over and clickable names for each exon-exon junction
if (active == 1){

  #Known junctions
  if (length(kji) > 0){
    for(i in kji){
      url_text=paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=", chr, ":", (model[i,"Start"]-50), "-", (model[i,"End"]+50), sep="")
      title_text = "Canonical Junction"
      desc1_text = paste(chr, ":", model[i,"Start"], "-", model[i,"End"], sep="")
      desc2_text = paste("Intron Size: ", model[i,"Size"], " bp", sep="")
      setSVGShapeToolTip(title=title_text, desc1=desc1_text, desc2=desc2_text)
      setSVGShapeURL(url_text, target="_blank") 
      text(x=model[i,"DrawMid"], y=bottom_label_main, labels=model[i,"Name"], cex=model_text_cex, font=2, col=NULL)
    }
  }
  #Alternative junctions
  if (length(aji)){
    for(i in aji){
      url_text=paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=", chr, ":", (model[i,"Start"]-50), "-", (model[i,"End"]+50), sep="")
      title_text = "Alternative Junction"
      desc1_text = paste(chr, ":", model[i,"Start"], "-", model[i,"End"], sep="")
      desc2_text = paste("Intron Size: ", model[i,"Size"], " bp", sep="")
      setSVGShapeToolTip(title=title_text, desc1=desc1_text, desc2=desc2_text)
      setSVGShapeURL(url_text, target="_blank") 
      text(x=model[i,"DrawMid"], y=top_label_main, labels=model[i,"Name"], cex=model_text_cex, font=2, col=NULL)
    }
  }
}

#Note the chromosomal end points being displayed
grand_start_text = paste(chr, ":", prettyNum(grand_start, big.mark=","), " (", strand, ")", sep="")
grand_end_text = paste(chr, ":", prettyNum(grand_end, big.mark=","), " (", strand, ")", sep="")
text(x=0, y=top_label_upper, labels=grand_start_text, cex=model_text_cex, font=2, pos=4)
text(x=grand_length, y=top_label_upper, labels=grand_end_text, cex=model_text_cex, font=2, pos=1)

#Add simple legends
legend1_pos = "topleft"
legend1_text = paste("Mutation (", mut_lib_name_list, ")", sep="")
legend1_cols = c(mutant_color)
legend1_ltys = c(NA)
legend1_lwds = c(NA)
legend1_pchs = c(8)
legend(legend1_pos, legend1_text, col=legend1_cols, lty=legend1_ltys, lwd=legend1_lwds, pch=legend1_pchs, cex=model_text_cex)

legend2_pos = "topright"
legend2_text = c("Canonical junction", "Alternative junction")
legend2_cols = c(kj_color, aj_color)
legend2_ltys = c(1, 2)
legend2_lwds = c(lwd3, lwd3)
legend2_pchs = c(NA, NA)
legend(legend2_pos, legend2_text, col=legend2_cols, lty=legend2_ltys, lwd=legend2_lwds, pch=legend2_pchs, cex=model_text_cex)

#Reset margins to defaults
par(mar=c(10,4,4,2)+0.1)

return(1);
}
#############################################################################################################################################################







