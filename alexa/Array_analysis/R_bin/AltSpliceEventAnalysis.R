#Experiment with different methods for identifying Alt. Splice Events
#Both Absolute (single sample) and differential (sample A vs. sample B)

#1-A.) Import the normalized data
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/Processed_Data"
dir(datadir)
setwd(datadir)
data_complete_norm = read.table("allProbes_Normalize_WithinHybeCondition.txt", header=T, quote="", sep="\t", comment.char="", as.is=3:5, na.strings='na')

#1-B.) FIRST REMOVE EXTREME VALUES (SATURATED SPOTS)
#If any probe is greater than 60,000 in either of the two experiments for a particular hybe condition
#then remove that probe from consideration

data_complete_norm_24 = data_complete_norm[which(data_complete_norm[,"LnCAP_24"] < 60000),]
data_complete_norm_24 = data_complete_norm_24[which(data_complete_norm_24[,"Brain_24"] < 60000),]

data_complete_norm_60 = data_complete_norm[which(data_complete_norm[,"LnCAP_60"] < 60000),]
data_complete_norm_60 = data_complete_norm_60[which(data_complete_norm_60[,"Brain_60"] < 60000),]

#1-C.) CONDUCTING BACKGROUND CORRECTIONS - Using Tm method for ALL probes
#BACKGROUND CORRECTION BASED ON INTENSITIES OBSERVED FOR NEGATIVE CONTROL PROBES (TM BASED)
#Get the control probes and their Tm for all ~4400 negative control probes (normalized values)
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

#Now for each probe, get the Tm of the probe, use the model to calculate the bg I and use that for correction

#LnCAP24
probe_tms_24 = data_complete_norm_24[,"Tm"]
probe_bg = theta0_ln24[1] + theta0_ln24[2]*probe_tms_24 + theta0_ln24[3] * (probe_tms_24)^2
data_complete_norm_24[,"LnCAP_24_bgc2_ratio"] = (data_complete_norm_24[, "LnCAP_24"]) / probe_bg

#Brain24
probe_bg = theta0_br24[1] + theta0_br24[2]*probe_tms_24 + theta0_br24[3] * (probe_tms_24)^2
data_complete_norm_24[,"Brain_24_bgc2_ratio"] = (data_complete_norm_24[, "Brain_24"]) / probe_bg

#LnCAP60
probe_tms_60 = data_complete_norm_60[,"Tm"]
probe_bg = theta0_ln60[1] + theta0_ln60[2]*probe_tms_60 + theta0_ln60[3] * (probe_tms_60)^2
data_complete_norm_60[,"LnCAP_60_bgc2_ratio"] = (data_complete_norm_60[, "LnCAP_60"]) / probe_bg

#Brain60
probe_bg = theta0_br60[1] + theta0_br60[2]*probe_tms_60 + theta0_br60[3] * (probe_tms_60)^2
data_complete_norm_60[,"Brain_60_bgc2_ratio"] = (data_complete_norm_60[, "Brain_60"]) / probe_bg

#Convert BGC values <1 to equal 1.
ln24_bgc2 = data_complete_norm_24[,"LnCAP_24_bgc2_ratio"]
br24_bgc2 = data_complete_norm_24[,"Brain_24_bgc2_ratio"]

ln24_bgc2[which(ln24_bgc2 < 1)] = 1
br24_bgc2[which(br24_bgc2 < 1)] = 1
ln60_bgc2 = data_complete_norm_60[,"LnCAP_60_bgc2_ratio"]
br60_bgc2 = data_complete_norm_60[,"Brain_60_bgc2_ratio"]

ln60_bgc2[which(ln60_bgc2 < 1)] = 1
br60_bgc2[which(br60_bgc2 < 1)] = 1

#1-D.) Calculate differential expression for all probes
DE_24_bgc2 = log2(ln24_bgc2) - log2(br24_bgc2)
data_complete_norm_24[,"DE_24_bgc2"] = DE_24_bgc2

DE_60_bgc2 = log2(ln60_bgc2) - log2(br60_bgc2)
data_complete_norm_60[,"DE_60_bgc2"] = DE_60_bgc2



#VISUALIZING PROBE INTENSITIES FOR A SINGLE GENE
#test on gene 6782
#display_gene(8703,"LnCAP")
display_gene = function(gene_id, exp){
  data_24_gene = data_complete_norm_24[which(data_complete_norm_24[,"AlexaGene_ID"]==gene_id),]

  exon_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"]=="Exon"),]
  canonical_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] == 0),]
  skip_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] > 0),]
  ei_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"] == "Exon-Intron" | data_24_gene[,"Probe_type"] == "Intron-Exon") ,]
  rj_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"] == "Rev-Junct_EE" | data_24_gene[,"Probe_type"] == "Rev-Junct_EI" | data_24_gene[,"Probe_type"] == "Rev-Junct_IE") ,]

  data_24 = paste(exp,"_24_bgc2_ratio",sep="")
  data_60 = paste(exp,"_60_bgc2_ratio",sep="")

  #Use two times mean plus SD of the RJ probes to draw an upper limit on the graph
  mux = mean(rj_24_gene[,data_24])
  sdx = sqrt(var(rj_24_gene[,data_24]))
  upper_lim = mux + (2*sdx)  

  min_x = min(data_24_gene[,"Tm"])
  max_x = max(data_24_gene[,"Tm"])
  min_y = min(data_24_gene[,data_24])
  max_y = max(data_24_gene[,data_24])
  title = paste(exp,"_24 - ALEXA GENE:", gene_id)
  plot(rj_24_gene[,"Tm"], rj_24_gene[,data_24], col="purple", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title, pch=18)
  points(skip_24_gene[,"Tm"], skip_24_gene[,data_24], col="red", pch=17)
  points(ei_24_gene[,"Tm"], ei_24_gene[,data_24], col="black", pch=17)
  points(canonical_24_gene[,"Tm"], canonical_24_gene[,data_24], col="blue", pch=16)
  points(exon_24_gene[,"Tm"], exon_24_gene[,data_24], col="dark green", pch=16)
  abline (h=upper_lim, lty=2)

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","ExonIntron","RJ"), col=c("dark green","blue","red","black","purple"), pch=c(16,16,17,17,18))

  #60-mer
  data_60_gene = data_complete_norm_60[which(data_complete_norm_60[,"AlexaGene_ID"]==gene_id),]

  exon_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"]=="Exon"),]
  canonical_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] == 0),]
  skip_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] > 0),]
  ei_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"] == "Exon-Intron" | data_60_gene[,"Probe_type"] == "Intron-Exon"),]
  rj_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"] == "Rev-Junct_EE" | data_60_gene[,"Probe_type"] == "Rev-Junct_EI" | data_60_gene[,"Probe_type"] == "Rev-Junct_IE") ,]

  #Use two times mean plus SD of the RJ probes to draw an upper limit on the graph
  mux = mean(rj_60_gene[,data_60])
  sdx = sqrt(var(rj_24_gene[,data_24]))
  upper_lim = mux + (2*sdx)

  min_x = min(data_60_gene[,"Tm"])
  max_x = max(data_60_gene[,"Tm"])
  min_y = min(data_60_gene[,data_60])
  max_y = max(data_60_gene[,data_60])
  title = paste(exp,"_60 - ALEXA GENE:", gene_id)
  windows()
  plot(rj_60_gene[,"Tm"], rj_60_gene[,data_60], col="purple", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title, pch=18)
  points(skip_60_gene[,"Tm"], skip_60_gene[,data_60], col="red", pch=17)
  points(ei_60_gene[,"Tm"], ei_60_gene[,data_60], col="black", pch=17)
  points(canonical_60_gene[,"Tm"], canonical_60_gene[,data_60], col="blue", pch=16)
  points(exon_60_gene[,"Tm"], exon_60_gene[,data_60], col="dark green", pch=16)
  abline (h=upper_lim, lty=2)

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","ExonIntron","RJ"), col=c("dark green","blue","red","black","purple"), pch=c(16,16,17,17,18))
}


#VISUALIZING DIFFERENTIAL PROBE INTENSITIES FOR A SINGLE GENE
#test on gene 6782
#display_gene_de(6782,"LnCAP","Brain")
display_gene_de = function(gene_id, exp1, exp2){
  data_24_gene = data_complete_norm_24[which(data_complete_norm_24[,"AlexaGene_ID"]==gene_id),]

  data_24_exp1 = paste(exp1,"_24_bgc2_ratio",sep="")
  data_24_exp2 = paste(exp2,"_24_bgc2_ratio",sep="")

  exon_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"]=="Exon"),]
  canonical_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] == 0),]
  skip_24_gene = data_24_gene[which(data_24_gene[,"Exons_skipped"] > 0),]
  ei_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"] == "Exon-Intron" | data_24_gene[,"Probe_type"] == "Intron-Exon") ,]
  rj_24_gene = data_24_gene[which(data_24_gene[,"Probe_type"] == "Rev-Junct_EE" | data_24_gene[,"Probe_type"] == "Rev-Junct_EI" | data_24_gene[,"Probe_type"] == "Rev-Junct_IE") ,]

  #Use two times mean plus/minus SD of the RJ probes to draw an upper/lower limit on the graph
  mux = mean(rj_24_gene[,"DE_24_bgc2"])
  sdx = sqrt(var(rj_24_gene[,"DE_24_bgc2"]))
  lower_lim = mux - (2*sdx)
  upper_lim = mux + (2*sdx)

  min_x = min(data_24_gene[,"Tm"])
  max_x = max(data_24_gene[,"Tm"])
  min_y = min(data_24_gene[,"DE_24_bgc2"])
  max_y = max(data_24_gene[,"DE_24_bgc2"])
  title = paste(exp1,"_vs_",exp2,"DE_24 - ALEXA GENE:", gene_id)
  plot(rj_24_gene[,"Tm"], rj_24_gene[,"DE_24_bgc2"], col="purple", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title, pch=18)
  points(skip_24_gene[,"Tm"], skip_24_gene[,"DE_24_bgc2"], col="red", pch=17)
  points(ei_24_gene[,"Tm"], ei_24_gene[,"DE_24_bgc2"], col="black", pch=17)
  points(canonical_24_gene[,"Tm"], canonical_24_gene[,"DE_24_bgc2"], col="blue", pch=16)
  points(exon_24_gene[,"Tm"], exon_24_gene[,"DE_24_bgc2"], col="dark green", pch=16)
  abline (h=lower_lim, lty=2)
  abline (h=upper_lim, lty=2)

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","ExonIntron","RJ"), col=c("dark green","blue","red","black","purple"), pch=c(16,16,17,17,18))

  #60-mer
  data_60_gene = data_complete_norm_60[which(data_complete_norm_60[,"AlexaGene_ID"]==gene_id),]

  data_60_exp1 = paste(exp1,"_60_bgc2_ratio",sep="")
  data_60_exp2 = paste(exp2,"_60_bgc2_ratio",sep="")

  exon_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"]=="Exon"),]
  skip_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] > 0),]
  canonical_60_gene = data_60_gene[which(data_60_gene[,"Exons_skipped"] == 0),]
  ei_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"] == "Exon-Intron" | data_60_gene[,"Probe_type"] == "Intron-Exon") ,]
  rj_60_gene = data_60_gene[which(data_60_gene[,"Probe_type"] == "Rev-Junct_EE" | data_60_gene[,"Probe_type"] == "Rev-Junct_EI" | data_60_gene[,"Probe_type"] == "Rev-Junct_IE") ,]

  #Use two times mean plus/minus SD of the RJ probes to draw an upper/lower limit on the graph
  mux = mean(rj_60_gene[,"DE_60_bgc2"])
  sdx = sqrt(var(rj_60_gene[,"DE_60_bgc2"]))
  lower_lim = mux - (2*sdx)
  upper_lim = mux + (2*sdx)

  min_x = min(data_60_gene[,"Tm"])
  max_x = max(data_60_gene[,"Tm"])
  min_y = min(data_60_gene[,"DE_60_bgc2"])
  max_y = max(data_60_gene[,"DE_60_bgc2"])
  title = paste(exp1,"_vs_",exp2,"DE_60 - ALEXA GENE:", gene_id)
  plot(rj_60_gene[,"Tm"], rj_60_gene[,"DE_60_bgc2"], col="purple", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="LOG2 Differential Expression Ratio", main=title, pch=18)
  points(ei_60_gene[,"Tm"], ei_60_gene[,"DE_60_bgc2"], col="black", pch=17)
  points(skip_60_gene[,"Tm"], skip_60_gene[,"DE_60_bgc2"], col="red", pch=17)
  points(canonical_60_gene[,"Tm"], canonical_60_gene[,"DE_60_bgc2"], col="blue", pch=16)
  points(exon_60_gene[,"Tm"], exon_60_gene[,"DE_60_bgc2"], col="dark green", pch=16)
  abline (h=lower_lim, lty=2)
  abline (h=upper_lim, lty=2)

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","ExonIntron","RJ"), col=c("dark green","blue","red","black","purple"), pch=c(16,16,17,17,18))
}



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





