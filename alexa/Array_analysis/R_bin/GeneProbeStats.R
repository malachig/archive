#Gene Level analysis

library(hexbin)
library(lattice)
library(ellipse)
library(affy)

1-A.) CONDUCTING BACKGROUND CORRECTIONS - USING RATIOS AND DIFFS
#Try a background correction where junction probes are corrected using reverse-junction probes
#and exon probes are corrected using an estimate based on negative control probes and a quadratic 
#model fit, as described above.
#Compare background correction methods where the estimate of background is subtracted versus using
#a ratio of the experimental probe over the control probe

#Import the normalized data
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/Processed_Data"
dir(datadir)
setwd(datadir)
data_complete_norm = read.table("summaryDataNimbleGenPilot_Normalize_WithinHybeCondition.txt", header=T, quote="", sep="\t", comment.char="", as.is=3:5, na.strings='na')

#1-B.)  FIRST REMOVE EXTREME VALUES (SATURATED SPOTS)
#If any probe is greater than 60,000 in either of the two experiments for a particular hybe condition
#then remove that probe from consideration

data_complete_norm_24 = data_complete_norm[which(data_complete_norm[,"LnCAP_24"] < 60000),]
data_complete_norm_24 = data_complete_norm_24[which(data_complete_norm_24[,"Brain_24"] < 60000),]
x = c(which(data_complete_norm_24[,"LnCAP_24_revJunct"] < 60000), which(is.na(data_complete_norm_24[,"LnCAP_24_revJunct"])))
data_complete_norm_24 = data_complete_norm_24[x,]
x = c(which(data_complete_norm_24[,"Brain_24_revJunct"] < 60000), which(is.na(data_complete_norm_24[,"Brain_24_revJunct"])))
data_complete_norm_24 = data_complete_norm_24[x,]

data_complete_norm_60 = data_complete_norm[which(data_complete_norm[,"LnCAP_60"] < 60000),]
data_complete_norm_60 = data_complete_norm_60[which(data_complete_norm_60[,"Brain_60"] < 60000),]
x = c(which(data_complete_norm_60[,"LnCAP_60_revJunct"] < 60000), which(is.na(data_complete_norm_60[,"LnCAP_60_revJunct"])))
data_complete_norm_60 = data_complete_norm_60[x,]
x = c(which(data_complete_norm_60[,"Brain_60_revJunct"] < 60000), which(is.na(data_complete_norm_60[,"Brain_60_revJunct"])))
data_complete_norm_60 = data_complete_norm_60[x,]

#1-C.) JUNCTION PROBES - BACKGROUND CORRECTION 

#Background correction of junction probes.  Try Junction/ReverseJunction and Junction-ReverseJunction
#Also try Affy's method: R = (PM-MM) / (PM+MM)
#Note that values of R close to 1 indicate presence, values close to 0 or negative indicate absence

#LnCAP24
data_complete_norm_24[,"LnCAP_24_bgc_ratio"] = data_complete_norm_24[,"LnCAP_24"] / data_complete_norm_24[,"LnCAP_24_revJunct"]
data_complete_norm_24[,"LnCAP_24_bgc_diff"] = data_complete_norm_24[,"LnCAP_24"] - data_complete_norm_24[,"LnCAP_24_revJunct"]
data_complete_norm_24[,"LnCAP_24_bgc_affyR"] = ((data_complete_norm_24[,"LnCAP_24"] - data_complete_norm_24[,"LnCAP_24_revJunct"]) / (data_complete_norm_24[,"LnCAP_24"] + data_complete_norm_24[,"LnCAP_24_revJunct"]))

#Brain24
data_complete_norm_24[,"Brain_24_bgc_ratio"] = data_complete_norm_24[,"Brain_24"] / data_complete_norm_24[,"Brain_24_revJunct"]
data_complete_norm_24[,"Brain_24_bgc_diff"] = data_complete_norm_24[,"Brain_24"] - data_complete_norm_24[,"Brain_24_revJunct"]
data_complete_norm_24[,"Brain_24_bgc_affyR"] = ((data_complete_norm_24[,"Brain_24"] - data_complete_norm_24[,"Brain_24_revJunct"]) / (data_complete_norm_24[,"Brain_24"] + data_complete_norm_24[,"Brain_24_revJunct"]))

#LnCAP60
data_complete_norm_60[,"LnCAP_60_bgc_ratio"] = data_complete_norm_60[,"LnCAP_60"] / data_complete_norm_60[,"LnCAP_60_revJunct"]
data_complete_norm_60[,"LnCAP_60_bgc_diff"] = data_complete_norm_60[,"LnCAP_60"] - data_complete_norm_60[,"LnCAP_60_revJunct"]
data_complete_norm_60[,"LnCAP_60_bgc_affyR"] = ((data_complete_norm_60[,"LnCAP_60"] - data_complete_norm_60[,"LnCAP_60_revJunct"]) / (data_complete_norm_60[,"LnCAP_60"] + data_complete_norm_60[,"LnCAP_60_revJunct"]))

#Brain60
data_complete_norm_60[,"Brain_60_bgc_ratio"] = data_complete_norm_60[,"Brain_60"] / data_complete_norm_60[,"Brain_60_revJunct"]
data_complete_norm_60[,"Brain_60_bgc_diff"] = data_complete_norm_60[,"Brain_60"] - data_complete_norm_60[,"Brain_60_revJunct"]
data_complete_norm_60[,"Brain_60_bgc_affyR"] = ((data_complete_norm_60[,"Brain_60"] - data_complete_norm_60[,"Brain_60_revJunct"]) / (data_complete_norm_60[,"Brain_60"] + data_complete_norm_60[,"Brain_60_revJunct"]))


#1-D.) EXON PROBES - BACKGROUND CORRECTION
#BACKGROUND CORRECTION BASED ON INTENSITIES OBSERVED FOR NEGATIVE CONTROL PROBES (TM BASED)
#Get the control probes and their Tm for all ~4400 negative control probes (normalized values)
dir(datadir)
setwd(datadir)
data_complete_nc = read.table("summaryDataNimbleGenPilot_AllNegativeControlProbes.txt", header=T, quote="", sep="\t", comment.char="", as.is=2:3, na.strings='na')

#Create an array of the Tm and I values for each experiment for all negative control probes
nc_intensities = matrix(0, nrow=length(data_complete_nc[,"Tm"]), ncol=5, dimnames=list(data_complete_nc[,"Probe_ID"], c("Tm","LnCAP_24","LnCAP_60","Brain_24","Brain_60")))
nc_intensities[,"Tm"] = data_complete_nc[,"Tm"]
nc_intensities[,"LnCAP_24"] = data_complete_nc[,"LnCAP_24"]
nc_intensities[,"LnCAP_60"] = data_complete_nc[,"LnCAP_60"]
nc_intensities[,"Brain_24"] = data_complete_nc[,"Brain_24"]
nc_intensities[,"Brain_60"] = data_complete_nc[,"Brain_60"]

#Divide all negative control probes into bins, and remove the outliers from each bin
Tm_bins_ln24 = getTmBinMedians(nc_intensities, "LnCAP_24", 59.0, 75.0, 0.1)
Tm_bins_ln60 = getTmBinMedians(nc_intensities, "LnCAP_60", 59.0, 75.0, 0.1)
Tm_bins_br24 = getTmBinMedians(nc_intensities, "Brain_24", 59.0, 75.0, 0.1)
Tm_bins_br60 = getTmBinMedians(nc_intensities, "Brain_60", 59.0, 75.0, 0.1)

theta0_ln24 = summary(lm(Tm_bins_ln24[,"median"]~Tm_bins_ln24[,"Tm_bin_centre"] + I((Tm_bins_ln24[,"Tm_bin_centre"])^2)))$coefficients
theta0_ln60 = summary(lm(Tm_bins_ln60[,"median"]~Tm_bins_ln60[,"Tm_bin_centre"] + I((Tm_bins_ln60[,"Tm_bin_centre"])^2)))$coefficients
theta0_br24 = summary(lm(Tm_bins_br24[,"median"]~Tm_bins_br24[,"Tm_bin_centre"] + I((Tm_bins_br24[,"Tm_bin_centre"])^2)))$coefficients
theta0_br60 = summary(lm(Tm_bins_br60[,"median"]~Tm_bins_br60[,"Tm_bin_centre"] + I((Tm_bins_br60[,"Tm_bin_centre"])^2)))$coefficients

#Now for each Exon probe, get the Tm of the probe, use the model to calculate the bg I and use that for correction
#Calculate the Affy Present/Absent call using the BG estimate from the negative control probes.  
#  Note:  Not sure if it is valid to do this

#LnCAP24
exon_probe_24_loc = which(data_complete_norm_24[,"Probe_type"]=="Exon")
exon_probe_tms = data_complete_norm_24[exon_probe_24_loc,"Tm"]
exon_probe_bg = theta0_ln24[1] + theta0_ln24[2]*exon_probe_tms + theta0_ln24[3] * (exon_probe_tms)^2
data_complete_norm_24[exon_probe_24_loc,"LnCAP_24_bgc_ratio"] = (data_complete_norm_24[exon_probe_24_loc, "LnCAP_24"]) / exon_probe_bg
data_complete_norm_24[exon_probe_24_loc,"LnCAP_24_bgc_diff"] = (data_complete_norm_24[exon_probe_24_loc, "LnCAP_24"]) - exon_probe_bg
data_complete_norm_24[exon_probe_24_loc,"LnCAP_24_bgc_affyR"] = (data_complete_norm_24[exon_probe_24_loc, "LnCAP_24"] - exon_probe_bg) / (data_complete_norm_24[exon_probe_24_loc, "LnCAP_24"] + exon_probe_bg) 

#Brain24
exon_probe_bg = theta0_br24[1] + theta0_br24[2]*exon_probe_tms + theta0_br24[3] * (exon_probe_tms)^2
data_complete_norm_24[exon_probe_24_loc,"Brain_24_bgc_ratio"] = (data_complete_norm_24[exon_probe_24_loc, "Brain_24"]) / exon_probe_bg
data_complete_norm_24[exon_probe_24_loc,"Brain_24_bgc_diff"] = (data_complete_norm_24[exon_probe_24_loc, "Brain_24"]) - exon_probe_bg
data_complete_norm_24[exon_probe_24_loc,"Brain_24_bgc_affyR"] = (data_complete_norm_24[exon_probe_24_loc, "Brain_24"] - exon_probe_bg) / (data_complete_norm_24[exon_probe_24_loc, "Brain_24"] + exon_probe_bg) 

#LnCAP60
exon_probe_60_loc = which(data_complete_norm_60[,"Probe_type"]=="Exon")
exon_probe_tms = data_complete_norm_60[exon_probe_60_loc,"Tm"]
exon_probe_bg = theta0_ln60[1] + theta0_ln60[2]*exon_probe_tms + theta0_ln60[3] * (exon_probe_tms)^2
data_complete_norm_60[exon_probe_60_loc,"LnCAP_60_bgc_ratio"] = (data_complete_norm_60[exon_probe_60_loc, "LnCAP_60"]) / exon_probe_bg
data_complete_norm_60[exon_probe_60_loc,"LnCAP_60_bgc_diff"] = (data_complete_norm_60[exon_probe_60_loc, "LnCAP_60"]) - exon_probe_bg
data_complete_norm_60[exon_probe_60_loc,"LnCAP_60_bgc_affyR"] = (data_complete_norm_60[exon_probe_60_loc, "LnCAP_60"] - exon_probe_bg) / (data_complete_norm_60[exon_probe_60_loc, "LnCAP_60"] + exon_probe_bg) 

#Brain60
exon_probe_bg = theta0_br60[1] + theta0_br60[2]*exon_probe_tms + theta0_br60[3] * (exon_probe_tms)^2
data_complete_norm_60[exon_probe_60_loc,"Brain_60_bgc_ratio"] = (data_complete_norm_60[exon_probe_60_loc, "Brain_60"]) / exon_probe_bg
data_complete_norm_60[exon_probe_60_loc,"Brain_60_bgc_diff"] = (data_complete_norm_60[exon_probe_60_loc, "Brain_60"]) - exon_probe_bg
data_complete_norm_60[exon_probe_60_loc,"Brain_60_bgc_affyR"] = (data_complete_norm_60[exon_probe_60_loc, "Brain_60"] - exon_probe_bg) / (data_complete_norm_60[exon_probe_60_loc, "Brain_60"] + exon_probe_bg) 

#1-E.) BGC2 - BACKGROUND CORRECTION OF ALL PROBES USING INTENSITIES OBSERVED FOR NEGATIVE CONTROL PROBES (TM BASED) 
#LnCAP24
probe_tms = data_complete_norm_24[,"Tm"]
probe_bg = theta0_ln24[1] + theta0_ln24[2]*probe_tms + theta0_ln24[3] * (probe_tms)^2
data_complete_norm_24[,"LnCAP_24_bgc2_ratio"] = (data_complete_norm_24[, "LnCAP_24"]) / probe_bg
data_complete_norm_24[,"LnCAP_24_bgc2_diff"] = (data_complete_norm_24[, "LnCAP_24"]) - probe_bg
data_complete_norm_24[,"LnCAP_24_bgc2_affyR"] = (data_complete_norm_24[, "LnCAP_24"] - probe_bg) / (data_complete_norm_24[, "LnCAP_24"] + probe_bg) 

#Brain24
probe_bg = theta0_br24[1] + theta0_br24[2]*probe_tms + theta0_br24[3] * (probe_tms)^2
data_complete_norm_24[,"Brain_24_bgc2_ratio"] = (data_complete_norm_24[, "Brain_24"]) / probe_bg
data_complete_norm_24[,"Brain_24_bgc2_diff"] = (data_complete_norm_24[, "Brain_24"]) - probe_bg
data_complete_norm_24[,"Brain_24_bgc2_affyR"] = (data_complete_norm_24[, "Brain_24"] - probe_bg) / (data_complete_norm_24[, "Brain_24"] + probe_bg) 

#LnCAP60
probe_tms = data_complete_norm_60[,"Tm"]
probe_bg = theta0_ln60[1] + theta0_ln60[2]*probe_tms + theta0_ln60[3] * (probe_tms)^2
data_complete_norm_60[,"LnCAP_60_bgc2_ratio"] = (data_complete_norm_60[, "LnCAP_60"]) / probe_bg
data_complete_norm_60[,"LnCAP_60_bgc2_diff"] = (data_complete_norm_60[, "LnCAP_60"]) - probe_bg
data_complete_norm_60[,"LnCAP_60_bgc2_affyR"] = (data_complete_norm_60[, "LnCAP_60"] - probe_bg) / (data_complete_norm_60[, "LnCAP_60"] + probe_bg) 

#Brain60
probe_bg = theta0_br60[1] + theta0_br60[2]*probe_tms + theta0_br60[3] * (probe_tms)^2
data_complete_norm_60[,"Brain_60_bgc2_ratio"] = (data_complete_norm_60[, "Brain_60"]) / probe_bg
data_complete_norm_60[,"Brain_60_bgc2_diff"] = (data_complete_norm_60[, "Brain_60"]) - probe_bg
data_complete_norm_60[,"Brain_60_bgc2_affyR"] = (data_complete_norm_60[, "Brain_60"] - probe_bg) / (data_complete_norm_60[, "Brain_60"] + probe_bg) 


#1-F.) DEAL WITH NEGATIVES OR RATIOS LESS THAN 1
#If a background corrected probe (ratio) is less than 1, reset it to 1
#If this ratio is smaller than 1 (i.e. the intensity for RJ was greater than that for J) reset value to 1
#Do the same thing for background corrected exon probes
#Similarly, for bgc values calculated as a difference, change negative values to 1.

#Since for this analysis we will be borrowing strength across many probes for each gene, try examining data
#WITHOUT reseting background corrected values.
#This means that a gene with many MM > PM will be penalized 
#(considered to have lower expression, less chance of present call) 

#data_complete_norm_24[which(data_complete_norm_24[,"LnCAP_24_bgc_ratio"] < 1),"LnCAP_24_bgc_ratio"] = 1 #Reset ratio < 1 to be 1
#data_complete_norm_24[which(data_complete_norm_24[,"Brain_24_bgc_ratio"] < 1),"Brain_24_bgc_ratio"] = 1 #Reset ratio < 1 to be 1
#data_complete_norm_60[which(data_complete_norm_60[,"LnCAP_60_bgc_ratio"] < 1),"LnCAP_60_bgc_ratio"] = 1 #Reset ratio < 1 to be 1
#data_complete_norm_60[which(data_complete_norm_60[,"Brain_60_bgc_ratio"] < 1),"Brain_60_bgc_ratio"] = 1 #Reset ratio < 1 to be 1

#data_complete_norm_24[which(data_complete_norm_24[,"LnCAP_24_bgc_diff"] <= 0),"LnCAP_24_bgc_diff"] = 1 #Reset diff < 0 to be 1
#data_complete_norm_24[which(data_complete_norm_24[,"Brain_24_bgc_diff"] <= 0),"Brain_24_bgc_diff"] = 1 #Reset diff < 0 to be 1
#data_complete_norm_60[which(data_complete_norm_60[,"LnCAP_60_bgc_diff"] <= 0),"LnCAP_60_bgc_diff"] = 1 #Reset diff < 0 to be 1
#data_complete_norm_60[which(data_complete_norm_60[,"Brain_60_bgc_diff"] <= 0),"Brain_60_bgc_diff"] = 1 #Reset diff < 0 to be 1


#2.) OBSERVING THE EFFECT OF PROBE POSITION WITHIN THE TRANSCRIPT ON OBSERVED INTENSITY

#2-A.) Hexbin plots of probe intensities versus relative position within gene 

#24-mer data
x = c(which(data_complete_norm_24[,"Probe_type"]=="Exon"), which(data_complete_norm_24[,"Exons_skipped"] == 0))
exon_canonical_24 = data_complete_norm_24[x,]
exon_canonical_24_large = exon_canonical_24[which(exon_canonical_24[,"TranscriptSize"] > 3000),]
percent_from_5prime = (exon_canonical_24_large[,"Dist_5prime"]/exon_canonical_24_large[,"TranscriptSize"])*100
hbin = hexbin(percent_from_5prime, exon_canonical_24_large[,"LnCAP_24"], xbins=50)
plot.hexbin(hbin,style="nested.centroids",lcex=0.4, xlab="Percent Distance from 5' End", ylab="Raw Intensity",main="LnCAP_24")

hbin = hexbin(percent_from_5prime, exon_canonical_24_large[,"Brain_24"], xbins=50)
plot.hexbin(hbin,style="nested.centroids",lcex=0.4, xlab="Percent Distance from 5' End", ylab="Raw Intensity",main="Brain_24")

#60-mer data
x = c(which(data_complete_norm_60[,"Probe_type"]=="Exon"), which(data_complete_norm_60[,"Exons_skipped"] == 0))
exon_canonical_60 = data_complete_norm_60[x,]
exon_canonical_60_large = exon_canonical_60[which(exon_canonical_60[,"TranscriptSize"] > 3000),]
percent_from_5prime = (exon_canonical_60_large[,"Dist_5prime"]/exon_canonical_60_large[,"TranscriptSize"])*100
hbin = hexbin(percent_from_5prime, exon_canonical_60_large[,"LnCAP_60"], xbins=50)
plot.hexbin(hbin,style="nested.centroids",lcex=0.4, xlab="Percent Distance from 5' End", ylab="Raw Intensity",main="LnCAP_60")

hbin = hexbin(percent_from_5prime, exon_canonical_60_large[,"Brain_60"], xbins=50)
plot.hexbin(hbin,style="nested.centroids",lcex=0.4, xlab="Percent Distance from 5' End", ylab="Raw Intensity",main="Brain_60")


#2-B.) Box plot comparisons of probe intensities for probes at either end
#What if we consider only probes corresponding to genes above some size (say 2000bp)
#Then plot the probe intensity against distance from the end of the transcript
#Select only the values that are within the first or last 500 bp of the transcript
#compare the intensities of these two groups, plot as two seperate data series against distance from end 
#Are these the intensities of these two populations significantly different?

#Get negative control probes for reference
nc_24 = data_complete_norm_24[which(data_complete_norm_24[,"Probe_type"]=="Control-Negative"),]
nc_60 = data_complete_norm_60[which(data_complete_norm_60[,"Probe_type"]=="Control-Negative"),]

#Get exon and canonical probes only
x = c(which(data_complete_norm_24[,"Probe_type"]=="Exon"), which(data_complete_norm_24[,"Exons_skipped"] == 0))
exon_canonical_24 = data_complete_norm_24[x,]
x = c(which(data_complete_norm_60[,"Probe_type"]=="Exon"), which(data_complete_norm_60[,"Exons_skipped"] == 0))
exon_canonical_60 = data_complete_norm_60[x,]

#Now select those probes corresponding to genes larger than 3000bp
exon_canonical_24_large = exon_canonical_24[which(exon_canonical_24[,"TranscriptSize"] > 3000),]
exon_canonical_60_large = exon_canonical_60[which(exon_canonical_60[,"TranscriptSize"] > 3000),]

#Now select only those probes that are within 500 bp of either end
x = which(exon_canonical_24_large[,"Dist_5prime"] <= 1000 | exon_canonical_24_large[,"Dist_3prime"] <= 1000)
exon_canonical_24_ends = exon_canonical_24_large[x,]
x = which(exon_canonical_60_large[,"Dist_5prime"] <= 1000 | exon_canonical_60_large[,"Dist_3prime"] <= 1000)
exon_canonical_60_ends = exon_canonical_60_large[x,]

par(mfrow=c(2,2))
plot(exon_canonical_24_ends[,"Dist_3prime"], exon_canonical_24_ends[,"LnCAP_24"])
plot(exon_canonical_24_ends[,"Dist_3prime"], exon_canonical_24_ends[,"Brain_24"])
plot(exon_canonical_60_ends[,"Dist_3prime"], exon_canonical_60_ends[,"LnCAP_60"])
plot(exon_canonical_60_ends[,"Dist_3prime"], exon_canonical_60_ends[,"Brain_60"])

#Box plots - Probes in 1st 500bp vs last 500bp - For transcripts larger than 3000bp only
first_500_24 = exon_canonical_24_large[which(exon_canonical_24_large[,"Dist_5prime"] <= 500),]
last_500_24 = exon_canonical_24_large[which(exon_canonical_24_large[,"Dist_3prime"] <= 500),]
first_500_60 = exon_canonical_60_large[which(exon_canonical_60_large[,"Dist_5prime"] <= 500),]
last_500_60 = exon_canonical_60_large[which(exon_canonical_60_large[,"Dist_3prime"] <= 500),]

z = list(first_500_24["LnCAP_24"],last_500_24["LnCAP_24"],first_500_24["Brain_24"],last_500_24["Brain_24"],
         first_500_60["LnCAP_60"],last_500_60["LnCAP_60"],first_500_60["Brain_60"],last_500_60["Brain_60"],
	   nc_24[,"LnCAP_24"],nc_24[,"Brain_24"],nc_60[,"LnCAP_60"],nc_60[,"Brain_60"])
names(z) = c("Ln24 5'end", "Ln24 3'end", "Br24 5'end","Br24 3'end", "Ln60 5'end", "Ln60 3'end", 
             "Br60 5'end","Br60 3'end", "NC_Ln24", "NC_Br24", "NC_Ln60", "NC_Br60")

par(mfrow=c(2,2), ps=10)
boxplot(z[c("Ln24 5'end","Ln24 3'end", "NC_Ln24")], main="LnCAP 24", ylab="Raw Intensity", ylim=c(0,25000))
boxplot(z[c("Br24 5'end","Br24 3'end", "NC_Br24")], main="Brain 24", ylab="Raw Intensity", ylim=c(0,25000))
boxplot(z[c("Ln60 5'end","Ln60 3'end", "NC_Ln60")], main="LnCAP 60", ylab="Raw Intensity", ylim=c(0,15000))
boxplot(z[c("Br60 5'end","Br60 3'end", "NC_Br60")], main="Brain 60", ylab="Raw Intensity", ylim=c(0,15000))





#2-C.) VISUALIZING PROBE INTENSITIES FOR A SINGLE GENE
#test on gene 6782
display_gene = function(gene_id, exp){
  data_24_gene = data_complete_norm_24[which(data_complete_norm_24[,"AlexaGene_ID"]==gene_id),]

  exon_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"]=="Exon"),]
  canonical_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] == 0),]
  skip_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] > 0),]
  ei_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"] == "Exon-Intron" | data_24_gene[,"Probe_type"] == "Intron-Exon") ,]

  data_24 = paste(exp,"_24_bgc2_ratio",sep="")
  data_60 = paste(exp,"_60_bgc2_ratio",sep="")

  min_x = min(data_24_gene[,"Tm"])
  max_x = max(data_24_gene[,"Tm"])
  min_y = min(data_24_gene[,data_24])
  max_y = max(data_24_gene[,data_24])
  title = paste(exp,"_24 - ALEXA GENE:", gene_id)
  plot(exon_24_gene[,"Tm"], exon_24_gene[,data_24], col="green", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title )
  points(canonical_24_gene[,"Tm"], canonical_24_gene[,data_24], col="blue")
  points(skip_24_gene[,"Tm"], skip_24_gene[,data_24], col="red")
  points(ei_24_gene[,"Tm"], ei_24_gene[,data_24], col="black")

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","Exon-Intron"), col=c("green","blue","red","black"), pch=c(1,1,1,1))

  #60-mer
  data_60_gene = data_complete_norm_60[which(data_complete_norm_60[,"AlexaGene_ID"]==gene_id),]

  exon_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"]=="Exon"),]
  canonical_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] == 0),]
  skip_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] > 0),]
  ei_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"] == "Exon-Intron" | data_60_gene[,"Probe_type"] == "Intron-Exon"),]

  min_x = min(data_60_gene[,"Tm"])
  max_x = max(data_60_gene[,"Tm"])
  min_y = min(data_60_gene[,data_60])
  max_y = max(data_60_gene[,data_60])
  title = paste(exp,"_60 - ALEXA GENE:", gene_id)
  windows()
  plot(exon_60_gene[,"Tm"], exon_60_gene[,data_60], col="green", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title)
  points(canonical_60_gene[,"Tm"], canonical_60_gene[,data_60], col="blue")
  points(skip_60_gene[,"Tm"], skip_60_gene[,data_60], col="red")
  points(ei_60_gene[,"Tm"], ei_60_gene[,data_60], col="black")

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","Exon-Intron"), col=c("green","blue","red","black"), pch=c(1,1,1,1))
}



#3-A.) CALCULATE PRESENCE/ABSENCE CALLS
#Try to assess the presence or absence of each gene, by combining all PM/MM values for a gene
#Consider only the exon and canonical probes

#24-MER DATA
x = c(which(data_complete_norm_24[,"Probe_type"]=="Exon"), which(data_complete_norm_24[,"Exons_skipped"] == 0))
exon_canonical_24 = data_complete_norm_24[x,]

#Get a unique list of all gene IDs
genes = unique(exon_canonical_24[,"AlexaGene_ID"])
genes = genes[(which(genes > 0))]

gene_stats_24 = matrix(0, nrow=length(genes), ncol=34, 
                dimnames=list(genes, c("probe_count","LnCAP_454_hits",
                  "mean_bgc_ratio_ln","sd_bgc_ratio_ln","pval_tt_bgc_ratio_ln","pval_wt_bgc_ratio_ln",
                  "mean_bgc_affyR_ln","sd_bgc_affyR_ln","pval_tt_bgc_affyR_ln","pval_wt_bgc_affyR_ln",
                  "mean_bgc2_ratio_ln","sd_bgc2_ratio_ln","pval_tt_bgc2_ratio_ln","pval_wt_bgc2_ratio_ln",
                  "mean_bgc2_affyR_ln","sd_bgc2_affyR_ln","pval_tt_bgc2_affyR_ln","pval_wt_bgc2_affyR_ln",
                  "mean_bgc_ratio_br","sd_bgc_ratio_br","pval_tt_bgc_ratio_br","pval_wt_bgc_ratio_br",
                  "mean_bgc_affyR_br","sd_bgc_affyR_br","pval_tt_bgc_affyR_br","pval_wt_bgc_affyR_br",
                  "mean_bgc2_ratio_br","sd_bgc2_ratio_br","pval_tt_bgc2_ratio_br","pval_wt_bgc2_ratio_br",
                  "mean_bgc2_affyR_br","sd_bgc2_affyR_br","pval_tt_bgc2_affyR_br","pval_wt_bgc2_affyR_br")))

for (i in 1:length(genes)){
#for(i in 1:10){
  #First count the number of exon and canonical probes for each gene
  gene_probes = exon_canonical_24[which(exon_canonical_24[,"AlexaGene_ID"]==genes[i]),]
  probe_count = length(gene_probes[,"Probe_ID"])
  gene_stats_24[i,"probe_count"] = probe_count

  if (probe_count == 1){
    next()
  }

  #BGC1 - Use reverse-junctions to correct junction probes and NC probes to correct exon probes

  #Now calculate the mean and SD of the background corrected ratio for these probes for LnCAP and Brain
  #Also try a T-test and Wilcoxon's test (as used by Affy) to evaluate the significance of the ratios
  #i.e. test whether they are significantly greater than 1
  gene_stats_24[i,"mean_bgc_ratio_ln"] = mean(gene_probes[,"LnCAP_24_bgc_ratio"])
  gene_stats_24[i,"sd_bgc_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc_ratio"]))
  gene_stats_24[i,"mean_bgc_ratio_br"] = mean(gene_probes[,"Brain_24_bgc_ratio"])
  gene_stats_24[i,"sd_bgc_ratio_br"] = sqrt(var(gene_probes[,"Brain_24_bgc_ratio"]))

  gene_stats_24[i,"pval_tt_bgc_ratio_ln"] = t.test(gene_probes[,"LnCAP_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_tt_bgc_ratio_br"] = t.test(gene_probes[,"Brain_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  gene_stats_24[i,"pval_wt_bgc_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_wt_bgc_ratio_br"] = wilcox.test(gene_probes[,"Brain_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now do the same calculations for the Affy Present/Absent metric values ((PM-MM)/(PM+MM))
  gene_stats_24[i,"mean_bgc_affyR_ln"] = mean(gene_probes[,"LnCAP_24_bgc_affyR"])
  gene_stats_24[i,"sd_bgc_affyR_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc_affyR"]))
  gene_stats_24[i,"mean_bgc_affyR_br"] = mean(gene_probes[,"Brain_24_bgc_affyR"])
  gene_stats_24[i,"sd_bgc_affyR_br"] = sqrt(var(gene_probes[,"Brain_24_bgc_affyR"]))

  gene_stats_24[i,"pval_tt_bgc_affyR_ln"] = t.test(gene_probes[,"LnCAP_24_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_24[i,"pval_tt_bgc_affyR_br"] = t.test(gene_probes[,"Brain_24_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  gene_stats_24[i,"pval_wt_bgc_affyR_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_24[i,"pval_wt_bgc_affyR_br"] = wilcox.test(gene_probes[,"Brain_24_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  #BGC2 - Use NC probes to correct ALL probes
  gene_stats_24[i,"mean_bgc2_ratio_ln"] = mean(gene_probes[,"LnCAP_24_bgc2_ratio"])
  gene_stats_24[i,"sd_bgc2_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc2_ratio"]))
  gene_stats_24[i,"mean_bgc2_ratio_br"] = mean(gene_probes[,"Brain_24_bgc2_ratio"])
  gene_stats_24[i,"sd_bgc2_ratio_br"] = sqrt(var(gene_probes[,"Brain_24_bgc2_ratio"]))

  gene_stats_24[i,"pval_tt_bgc2_ratio_ln"] = t.test(gene_probes[,"LnCAP_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_tt_bgc2_ratio_br"] = t.test(gene_probes[,"Brain_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  gene_stats_24[i,"pval_wt_bgc2_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_wt_bgc2_ratio_br"] = wilcox.test(gene_probes[,"Brain_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now do the same calculations for the Affy Present/Absent metric values ((PM-MM)/(PM+MM))
  gene_stats_24[i,"mean_bgc2_affyR_ln"] = mean(gene_probes[,"LnCAP_24_bgc2_affyR"])
  gene_stats_24[i,"sd_bgc2_affyR_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc2_affyR"]))
  gene_stats_24[i,"mean_bgc2_affyR_br"] = mean(gene_probes[,"Brain_24_bgc2_affyR"])
  gene_stats_24[i,"sd_bgc2_affyR_br"] = sqrt(var(gene_probes[,"Brain_24_bgc2_affyR"]))

  gene_stats_24[i,"pval_tt_bgc2_affyR_ln"] = t.test(gene_probes[,"LnCAP_24_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_24[i,"pval_tt_bgc2_affyR_br"] = t.test(gene_probes[,"Brain_24_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  gene_stats_24[i,"pval_wt_bgc2_affyR_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_24[i,"pval_wt_bgc2_affyR_br"] = wilcox.test(gene_probes[,"Brain_24_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
}




#60-MER DATA
x = c(which(data_complete_norm_60[,"Probe_type"]=="Exon"), which(data_complete_norm_60[,"Exons_skipped"] == 0))
exon_canonical_60 = data_complete_norm_60[x,]

#Get a unique list of all gene IDs
genes = unique(exon_canonical_60[,"AlexaGene_ID"])
genes = genes[(which(genes > 0))]

gene_stats_60 = matrix(0, nrow=length(genes), ncol=34, 
                dimnames=list(genes, c("probe_count","LnCAP_454_hits",
                  "mean_bgc_ratio_ln","sd_bgc_ratio_ln","pval_tt_bgc_ratio_ln","pval_wt_bgc_ratio_ln",
                  "mean_bgc_affyR_ln","sd_bgc_affyR_ln","pval_tt_bgc_affyR_ln","pval_wt_bgc_affyR_ln",
                  "mean_bgc2_ratio_ln","sd_bgc2_ratio_ln","pval_tt_bgc2_ratio_ln","pval_wt_bgc2_ratio_ln",
                  "mean_bgc2_affyR_ln","sd_bgc2_affyR_ln","pval_tt_bgc2_affyR_ln","pval_wt_bgc2_affyR_ln",
                  "mean_bgc_ratio_br","sd_bgc_ratio_br","pval_tt_bgc_ratio_br","pval_wt_bgc_ratio_br",
                  "mean_bgc_affyR_br","sd_bgc_affyR_br","pval_tt_bgc_affyR_br","pval_wt_bgc_affyR_br",
                  "mean_bgc2_ratio_br","sd_bgc2_ratio_br","pval_tt_bgc2_ratio_br","pval_wt_bgc2_ratio_br",
                  "mean_bgc2_affyR_br","sd_bgc2_affyR_br","pval_tt_bgc2_affyR_br","pval_wt_bgc2_affyR_br")))

for (i in 1:length(genes)){
  #First count the number of exon and canonical probes for each gene
  gene_probes = exon_canonical_60[which(exon_canonical_60[,"AlexaGene_ID"]==genes[i]),]
  probe_count = length(gene_probes[,"Probe_ID"])
  gene_stats_60[i,"probe_count"] = probe_count

  if (probe_count == 1){
    next()
  }

  #BGC1 - Use reverse-junctions to correct junction probes and NC probes to correct exon probes

  #Now calculate the mean and SD of the background corrected ratio for these probes for LnCAP and Brain
  #Also try a T-test and Wilcoxon's test (as used by Affy) to evaluate the significance of the ratios
  #i.e. test whether they are significantly greater than 1
  gene_stats_60[i,"mean_bgc_ratio_ln"] = mean(gene_probes[,"LnCAP_60_bgc_ratio"])
  gene_stats_60[i,"sd_bgc_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc_ratio"]))
  gene_stats_60[i,"mean_bgc_ratio_br"] = mean(gene_probes[,"Brain_60_bgc_ratio"])
  gene_stats_60[i,"sd_bgc_ratio_br"] = sqrt(var(gene_probes[,"Brain_60_bgc_ratio"]))

  gene_stats_60[i,"pval_tt_bgc_ratio_ln"] = t.test(gene_probes[,"LnCAP_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_tt_bgc_ratio_br"] = t.test(gene_probes[,"Brain_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  gene_stats_60[i,"pval_wt_bgc_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_wt_bgc_ratio_br"] = wilcox.test(gene_probes[,"Brain_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now do the same calculations for the Affy Present/Absent metric values ((PM-MM)/(PM+MM))
  gene_stats_60[i,"mean_bgc_affyR_ln"] = mean(gene_probes[,"LnCAP_60_bgc_affyR"])
  gene_stats_60[i,"sd_bgc_affyR_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc_affyR"]))
  gene_stats_60[i,"mean_bgc_affyR_br"] = mean(gene_probes[,"Brain_60_bgc_affyR"])
  gene_stats_60[i,"sd_bgc_affyR_br"] = sqrt(var(gene_probes[,"Brain_60_bgc_affyR"]))

  gene_stats_60[i,"pval_tt_bgc_affyR_ln"] = t.test(gene_probes[,"LnCAP_60_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_60[i,"pval_tt_bgc_affyR_br"] = t.test(gene_probes[,"Brain_60_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  gene_stats_60[i,"pval_wt_bgc_affyR_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_60[i,"pval_wt_bgc_affyR_br"] = wilcox.test(gene_probes[,"Brain_60_bgc_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  #BGC2 - Use NC probes to correct ALL probes
  gene_stats_60[i,"mean_bgc2_ratio_ln"] = mean(gene_probes[,"LnCAP_60_bgc2_ratio"])
  gene_stats_60[i,"sd_bgc2_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc2_ratio"]))
  gene_stats_60[i,"mean_bgc2_ratio_br"] = mean(gene_probes[,"Brain_60_bgc2_ratio"])
  gene_stats_60[i,"sd_bgc2_ratio_br"] = sqrt(var(gene_probes[,"Brain_60_bgc2_ratio"]))

  gene_stats_60[i,"pval_tt_bgc2_ratio_ln"] = t.test(gene_probes[,"LnCAP_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_tt_bgc2_ratio_br"] = t.test(gene_probes[,"Brain_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  gene_stats_60[i,"pval_wt_bgc2_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_wt_bgc2_ratio_br"] = wilcox.test(gene_probes[,"Brain_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now do the same calculations for the Affy Present/Absent metric values ((PM-MM)/(PM+MM))
  gene_stats_60[i,"mean_bgc2_affyR_ln"] = mean(gene_probes[,"LnCAP_60_bgc2_affyR"])
  gene_stats_60[i,"sd_bgc2_affyR_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc2_affyR"]))
  gene_stats_60[i,"mean_bgc2_affyR_br"] = mean(gene_probes[,"Brain_60_bgc2_affyR"])
  gene_stats_60[i,"sd_bgc2_affyR_br"] = sqrt(var(gene_probes[,"Brain_60_bgc2_affyR"]))

  gene_stats_60[i,"pval_tt_bgc2_affyR_ln"] = t.test(gene_probes[,"LnCAP_60_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_60[i,"pval_tt_bgc2_affyR_br"] = t.test(gene_probes[,"Brain_60_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value

  gene_stats_60[i,"pval_wt_bgc2_affyR_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
  gene_stats_60[i,"pval_wt_bgc2_affyR_br"] = wilcox.test(gene_probes[,"Brain_60_bgc2_affyR"], paired=FALSE, alternative = "greater", mu = 0)$p.value
}


#3-B.) IMPORT LNCAP 454 DATA AND ASSOCIATE WITH EACH GENE
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/ExpressionValidation"
dir(datadir)
setwd(datadir)

data_LnCAP_454 = read.table("prototypeGenes_454hits_joined.txt", header=T, quote="", sep="\t", comment.char="")
dimnames(data_LnCAP_454) = list(data_LnCAP_454[,"AlexaGene_ID"], c("AlexaGene_ID","LnCAP_454_hits"))

gene_names = dimnames(data_LnCAP_454)[[1]]
gene_stats_24[gene_names,"LnCAP_454_hits"] = data_LnCAP_454[gene_names,"LnCAP_454_hits"]
gene_stats_60[gene_names,"LnCAP_454_hits"] = data_LnCAP_454[gene_names,"LnCAP_454_hits"]


#Compare the present/absent calls calculated above to the hit counts for LnCAP 454 data
#Try comparing only the genes with at least one 454 hit 
#P-values
genes_454 = which(gene_stats_24[,"LnCAP_454_hits"] > 0)
par(mfrow=c(2,2))
plot(gene_stats_24[genes_454,"pval_wt_bgc_ratio_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_24[genes_454,"pval_wt_bgc_ratio_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"pval_wt_bgc_affyR_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_24[genes_454,"pval_wt_bgc_affyR_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"pval_wt_bgc2_ratio_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_24[genes_454,"pval_wt_bgc2_ratio_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"pval_wt_bgc2_affyR_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_24[genes_454,"pval_wt_bgc2_affyR_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")

genes_454 = which(gene_stats_60[,"LnCAP_454_hits"] > 0)
par(mfrow=c(2,2))
plot(gene_stats_60[genes_454,"pval_wt_bgc_ratio_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_60[genes_454,"pval_wt_bgc_ratio_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"pval_wt_bgc_affyR_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_60[genes_454,"pval_wt_bgc_affyR_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"pval_wt_bgc2_ratio_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_60[genes_454,"pval_wt_bgc2_ratio_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"pval_wt_bgc2_affyR_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"])
cor(x=gene_stats_60[genes_454,"pval_wt_bgc2_affyR_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")

#BGC mean values
genes_454 = which(gene_stats_24[,"LnCAP_454_hits"] > 0)
par(mfrow=c(2,2))
plot(gene_stats_24[genes_454,"mean_bgc_ratio_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC Ratio", ylab="LnCAP 454 Hits", main="Gene Expression - 24mer")
cor(x=gene_stats_24[genes_454,"mean_bgc_ratio_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"mean_bgc_affyR_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC AffyR", ylab="LnCAP 454 Hits", main="Gene Expression - 24mer")
cor(x=gene_stats_24[genes_454,"mean_bgc_affyR_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"mean_bgc2_ratio_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC2 Ratio", ylab="LnCAP 454 Hits", main="Gene Expression - 24mer")
cor(x=gene_stats_24[genes_454,"mean_bgc2_ratio_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_24[genes_454,"mean_bgc2_affyR_ln"], gene_stats_24[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC2 AffyR", ylab="LnCAP 454 Hits", main="Gene Expression - 24mer")
cor(x=gene_stats_24[genes_454,"mean_bgc2_affyR_ln"], y=gene_stats_24[genes_454,"LnCAP_454_hits"], method="pearson")

genes_454 = which(gene_stats_60[,"LnCAP_454_hits"] > 0)
par(mfrow=c(2,2))
plot(gene_stats_60[genes_454,"mean_bgc_ratio_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC Ratio", ylab="LnCAP 454 Hits", main="Gene Expression - 60mer")
cor(x=gene_stats_60[genes_454,"mean_bgc_ratio_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"mean_bgc_affyR_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC AffyR", ylab="LnCAP 454 Hits", main="Gene Expression - 60mer")
cor(x=gene_stats_60[genes_454,"mean_bgc_affyR_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"mean_bgc2_ratio_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC2 Ratio", ylab="LnCAP 454 Hits", main="Gene Expression - 60mer")
cor(x=gene_stats_60[genes_454,"mean_bgc2_ratio_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")
plot(gene_stats_60[genes_454,"mean_bgc2_affyR_ln"], gene_stats_60[genes_454,"LnCAP_454_hits"], xlab="LnCAP Mean BGC2 AffyR", ylab="LnCAP 454 Hits", main="Gene Expression - 60mer")
cor(x=gene_stats_60[genes_454,"mean_bgc2_affyR_ln"], y=gene_stats_60[genes_454,"LnCAP_454_hits"], method="pearson")

#Compare the BGC values of 24-mer versus 60-mer (for Brain and LnCAP)
par(mfrow=c(2,2))
plot (gene_stats_24[,"mean_bgc_ratio_ln"], gene_stats_60[,"mean_bgc_ratio_ln"], xlab="LnCAP_24", ylab="LnCAP_60", main="Mean BGC Ratios")
plot (gene_stats_24[,"mean_bgc_affyR_ln"], gene_stats_60[,"mean_bgc_affyR_ln"], xlab="LnCAP_24", ylab="LnCAP_60", main="Mean BGC AffyR")
plot (gene_stats_24[,"mean_bgc2_ratio_ln"], gene_stats_60[,"mean_bgc2_ratio_ln"], xlab="LnCAP_24", ylab="LnCAP_60", main="Mean BGC2 Ratios")
plot (gene_stats_24[,"mean_bgc2_affyR_ln"], gene_stats_60[,"mean_bgc2_affyR_ln"], xlab="LnCAP_24", ylab="LnCAP_60", main="Mean BGC2 AffyR")

par(mfrow=c(2,2))
plot (gene_stats_24[,"mean_bgc_ratio_br"], gene_stats_60[,"mean_bgc_ratio_br"], xlab="Brain_24", ylab="Brain_60", main="Mean BGC Ratios")
plot (gene_stats_24[,"mean_bgc_affyR_br"], gene_stats_60[,"mean_bgc_affyR_br"], xlab="Brain_24", ylab="Brain_60", main="Mean BGC AffyR")
plot (gene_stats_24[,"mean_bgc2_ratio_br"], gene_stats_60[,"mean_bgc2_ratio_br"], xlab="Brain_24", ylab="Brain_60", main="Mean BGC2 Ratios")
plot (gene_stats_24[,"mean_bgc2_affyR_br"], gene_stats_60[,"mean_bgc2_affyR_br"], xlab="Brain_24", ylab="Brain_60", main="Mean BGC2 AffyR")

#The same sample in different hybe conditions seems to correlate pretty well
#BUT, how well do different samples compare with the same hybe conditions
par(mfrow=c(2,2))
plot (gene_stats_24[,"mean_bgc_ratio_ln"], gene_stats_24[,"mean_bgc_ratio_br"])
plot (gene_stats_24[,"mean_bgc_affyR_ln"], gene_stats_24[,"mean_bgc_affyR_br"])
plot (gene_stats_60[,"mean_bgc2_ratio_ln"], gene_stats_60[,"mean_bgc2_ratio_br"])
plot (gene_stats_60[,"mean_bgc2_affyR_ln"], gene_stats_60[,"mean_bgc2_affyR_br"])

#Try a correlation matrix for the mean gene ratios and AffyR gene values
#BGC1 - Ratios
ratios = matrix(0,nrow=length(gene_stats_24[,"probe_count"]), ncol=4, dimnames=list(dimnames(gene_stats_24)[[1]], c("LnCAP_24","Brain_24","LnCAP_60","Brain_60"))) 
ratios[,"LnCAP_24"] = gene_stats_24[,"mean_bgc_ratio_ln"]
ratios[,"Brain_24"] = gene_stats_24[,"mean_bgc_ratio_br"]
ratios[,"LnCAP_60"] = gene_stats_60[,"mean_bgc_ratio_ln"]
ratios[,"Brain_60"] = gene_stats_60[,"mean_bgc_ratio_br"]
r_pear_mean_gene_ratios = cor(ratios, method="pearson")
plotcorr(r_pear_mean_gene_ratios, col="dark green", main="Pearson Corr - Mean BGC Ratios")
r_pear_mean_gene_ratios

#BGC1 - AffyR
affyR = matrix(0,nrow=length(gene_stats_24[,"probe_count"]), ncol=4, dimnames=list(dimnames(gene_stats_24)[[1]], c("LnCAP_24","Brain_24","LnCAP_60","Brain_60"))) 
affyR[,"LnCAP_24"] = gene_stats_24[,"mean_bgc_affyR_ln"]
affyR[,"Brain_24"] = gene_stats_24[,"mean_bgc_affyR_br"]
affyR[,"LnCAP_60"] = gene_stats_60[,"mean_bgc_affyR_ln"]
affyR[,"Brain_60"] = gene_stats_60[,"mean_bgc_affyR_br"]
r_pear_mean_gene_affyR = cor(affyR, method="pearson")
plotcorr(r_pear_mean_gene_affyR, col="dark green", main="Pearson Corr - Mean BGC AffyR")
r_pear_mean_gene_affyR

#BGC2 - Ratios
ratios = matrix(0,nrow=length(gene_stats_24[,"probe_count"]), ncol=4, dimnames=list(dimnames(gene_stats_24)[[1]], c("LnCAP_24","Brain_24","LnCAP_60","Brain_60"))) 
ratios[,"LnCAP_24"] = gene_stats_24[,"mean_bgc2_ratio_ln"]
ratios[,"Brain_24"] = gene_stats_24[,"mean_bgc2_ratio_br"]
ratios[,"LnCAP_60"] = gene_stats_60[,"mean_bgc2_ratio_ln"]
ratios[,"Brain_60"] = gene_stats_60[,"mean_bgc2_ratio_br"]
r_pear_mean_gene_ratios = cor(ratios, method="pearson")
plotcorr(r_pear_mean_gene_ratios, col="dark green", main="Pearson Corr - Mean BGC2 Ratios")
r_pear_mean_gene_ratios

#BGC2 - AffyR
affyR = matrix(0,nrow=length(gene_stats_24[,"probe_count"]), ncol=4, dimnames=list(dimnames(gene_stats_24)[[1]], c("LnCAP_24","Brain_24","LnCAP_60","Brain_60"))) 
affyR[,"LnCAP_24"] = gene_stats_24[,"mean_bgc2_affyR_ln"]
affyR[,"Brain_24"] = gene_stats_24[,"mean_bgc2_affyR_br"]
affyR[,"LnCAP_60"] = gene_stats_60[,"mean_bgc2_affyR_ln"]
affyR[,"Brain_60"] = gene_stats_60[,"mean_bgc2_affyR_br"]
r_pear_mean_gene_affyR = cor(affyR, method="pearson")
plotcorr(r_pear_mean_gene_affyR, col="dark green", main="Pearson Corr - Mean BGC2 AffyR")
r_pear_mean_gene_affyR



#Try a corelation of the BGC log2 DE values between 24-mer and 60-mer hybe conditions
r_pear_de = cor(x=gene_stats_24[,"mean_log2_DE"], y=gene_stats_60[,"mean_log2_DE"], method="pearson")

#Try a correlation matrix between mean_bgc_ratio for Ln and Br, 24-mer and 60-mer
#For all genes, then for only the top DE genes 



#Write out the gene level summary statistics
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/summary_out"
dir(datadir)
setwd(datadir)
write.table(gene_stats_24[,], file="GeneStats_24.txt", sep="\t", eol="\n", row.names=TRUE, quote=FALSE)
write.table(gene_stats_60[,], file="GeneStats_60.txt", sep="\t", eol="\n", row.names=TRUE, quote=FALSE)



###################################################################################################
#Create a function to do this for any set of negative control probes with a range of Tm
#Provide an array containing Tm in the first column, followed by labeled columns of experimental I values
#Specify the experiment name to consider, specify the min_tm and max_tm possible
#Specify the bin_size to use
#Return an array of median intensity values for each bin along with the max and min Tm of each bin
##################################################################################################
getTmBinMedians = function(tm_I_array,exp_name,min_tm, max_tm, bin_size){
  tm_diff=(abs(max_tm-min_tm))
  bin_count = round((tm_diff/bin_size), digits=0)

  Tm_start = min_tm
  nc_median_tm = matrix(0,nrow=bin_count,ncol=5, dimnames=list(1:bin_count, c("Tm_start","Tm_end","Tm_bin_centre","median","mean")))

  for(i in 1:bin_count){
    x=which((tm_I_array[,"Tm"] > Tm_start) & (tm_I_array[,"Tm"] < (Tm_start+bin_size)))
    y = tm_I_array[x,exp_name]

    #Identify all the I values in this Tm bin that are NOT outliers
    #Where outliers are defined as values that are less than Q1-(IQR*.5) or more than Q3+(IQR*.5)
    q1 = quantile(y, 0.25)
    q3 = quantile(y, 0.75)
    iqr = (IQR(y))*0.5
    lower_limit = q1 - iqr
    upper_limit = q3 + iqr

    non_outliers = which((y > lower_limit) & (y < upper_limit))
    bin_I = y[non_outliers]

    med_I = median(bin_I)
    mean_I = mean(bin_I)

    Tm_bin_centre = Tm_start + ((abs((Tm_start+bin_size)-Tm_start))/2)

    nc_median_tm[i,"median"] = med_I
    nc_median_tm[i,"mean"] = mean_I
    nc_median_tm[i,"Tm_start"] = Tm_start
    nc_median_tm[i,"Tm_end"] = Tm_start+bin_size
    nc_median_tm[i,"Tm_bin_centre"] = Tm_bin_centre
    Tm_start = Tm_start+bin_size
  }
  return(nc_median_tm)
}
####################################################################################################



