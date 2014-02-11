#USE NORMALIZED PROBE DATA AS INPUT
#CONDUCT BACKGROUND CORRECTION AND CALCULATE DIFFERENTIAL EXPRESSION FOR PROBES - OUTPUT TO FILE
#ALSO CALCULATE DIFFERENTIAL EXPRESSION AT THE GENE LEVEL

#1-A.) Import the normalized data
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/Processed_Data"
dir(datadir)
setwd(datadir)
data_complete_norm = read.table("summaryDataNimbleGenPilot_Normalize_WithinHybeCondition.txt", header=T, quote="", sep="\t", comment.char="", as.is=3:5, na.strings='na')

#1-B.) FIRST REMOVE EXTREME VALUES (SATURATED SPOTS)
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


#1-C.) CONDUCTING BACKGROUND CORRECTIONS - USING RATIOS
#BGC1 - JUNCTION PROBES taken as ratio over reverse-junction
data_complete_norm_24[,"LnCAP_24_bgc_ratio"] = data_complete_norm_24[,"LnCAP_24"] / data_complete_norm_24[,"LnCAP_24_revJunct"]
data_complete_norm_24[,"Brain_24_bgc_ratio"] = data_complete_norm_24[,"Brain_24"] / data_complete_norm_24[,"Brain_24_revJunct"]
data_complete_norm_60[,"LnCAP_60_bgc_ratio"] = data_complete_norm_60[,"LnCAP_60"] / data_complete_norm_60[,"LnCAP_60_revJunct"]
data_complete_norm_60[,"Brain_60_bgc_ratio"] = data_complete_norm_60[,"Brain_60"] / data_complete_norm_60[,"Brain_60_revJunct"]

#BGC1 - EXON PROBES taken as ratio over an estimate calculated using a bin of probes with the same Tm
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

#Brain24
exon_probe_bg = theta0_br24[1] + theta0_br24[2]*exon_probe_tms + theta0_br24[3] * (exon_probe_tms)^2
data_complete_norm_24[exon_probe_24_loc,"Brain_24_bgc_ratio"] = (data_complete_norm_24[exon_probe_24_loc, "Brain_24"]) / exon_probe_bg

#LnCAP60
exon_probe_60_loc = which(data_complete_norm_60[,"Probe_type"]=="Exon")
exon_probe_tms = data_complete_norm_60[exon_probe_60_loc,"Tm"]
exon_probe_bg = theta0_ln60[1] + theta0_ln60[2]*exon_probe_tms + theta0_ln60[3] * (exon_probe_tms)^2
data_complete_norm_60[exon_probe_60_loc,"LnCAP_60_bgc_ratio"] = (data_complete_norm_60[exon_probe_60_loc, "LnCAP_60"]) / exon_probe_bg

#Brain60
exon_probe_bg = theta0_br60[1] + theta0_br60[2]*exon_probe_tms + theta0_br60[3] * (exon_probe_tms)^2
data_complete_norm_60[exon_probe_60_loc,"Brain_60_bgc_ratio"] = (data_complete_norm_60[exon_probe_60_loc, "Brain_60"]) / exon_probe_bg

#1-D.) BGC2 - BACKGROUND CORRECTION OF ALL PROBES USING INTENSITIES OBSERVED FOR NEGATIVE CONTROL PROBES (TM BASED) 
#LnCAP24
probe_tms = data_complete_norm_24[,"Tm"]
probe_bg = theta0_ln24[1] + theta0_ln24[2]*probe_tms + theta0_ln24[3] * (probe_tms)^2
data_complete_norm_24[,"LnCAP_24_bgc2_ratio"] = (data_complete_norm_24[, "LnCAP_24"]) / probe_bg

#Brain24
probe_bg = theta0_br24[1] + theta0_br24[2]*probe_tms + theta0_br24[3] * (probe_tms)^2
data_complete_norm_24[,"Brain_24_bgc2_ratio"] = (data_complete_norm_24[, "Brain_24"]) / probe_bg

#LnCAP60
probe_tms = data_complete_norm_60[,"Tm"]
probe_bg = theta0_ln60[1] + theta0_ln60[2]*probe_tms + theta0_ln60[3] * (probe_tms)^2
data_complete_norm_60[,"LnCAP_60_bgc2_ratio"] = (data_complete_norm_60[, "LnCAP_60"]) / probe_bg

#Brain60
probe_bg = theta0_br60[1] + theta0_br60[2]*probe_tms + theta0_br60[3] * (probe_tms)^2
data_complete_norm_60[,"Brain_60_bgc2_ratio"] = (data_complete_norm_60[, "Brain_60"]) / probe_bg

#1-E.) CALCULATE DIFFERENTIAL EXPRESSION ESTIMATES (LnCAP vs Brain) FOR ALL PROBES
#BGC1
ln24_bgc = data_complete_norm_24[,"LnCAP_24_bgc_ratio"]
br24_bgc = data_complete_norm_24[,"Brain_24_bgc_ratio"]

ln60_bgc = data_complete_norm_60[,"LnCAP_60_bgc_ratio"]
br60_bgc = data_complete_norm_60[,"Brain_60_bgc_ratio"]

#RESET VALUES < 1 TO 1
ln24_bgc[which(ln24_bgc < 1)] = 1
br24_bgc[which(br24_bgc < 1)] = 1
DE_24_bgc = log2(ln24_bgc) - log2(br24_bgc)
data_complete_norm_24[,"DE_24_bgc"] = DE_24_bgc

ln60_bgc[which(ln60_bgc < 1)] = 1
br60_bgc[which(br60_bgc < 1)] = 1
DE_60_bgc = log2(ln60_bgc) - log2(br60_bgc)
data_complete_norm_60[,"DE_60_bgc"] = DE_60_bgc

#BGC2
ln24_bgc2 = data_complete_norm_24[,"LnCAP_24_bgc2_ratio"]
br24_bgc2 = data_complete_norm_24[,"Brain_24_bgc2_ratio"]

ln60_bgc2 = data_complete_norm_60[,"LnCAP_60_bgc2_ratio"]
br60_bgc2 = data_complete_norm_60[,"Brain_60_bgc2_ratio"]

#RESET VALUES < 1 TO 1
ln24_bgc2[which(ln24_bgc2 < 1)] = 1
br24_bgc2[which(br24_bgc2 < 1)] = 1
DE_24_bgc2 = log2(ln24_bgc2) - log2(br24_bgc2)
data_complete_norm_24[,"DE_24_bgc2"] = DE_24_bgc2

ln60_bgc2[which(ln60_bgc2 < 1)] = 1
br60_bgc2[which(br60_bgc2 < 1)] = 1
DE_60_bgc2 = log2(ln60_bgc2) - log2(br60_bgc2)
data_complete_norm_60[,"DE_60_bgc2"] = DE_60_bgc2

#1-F.) OUTPUT PROBE VALUES TO FILE
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/summary_out"
dir(datadir)
setwd(datadir)
write.table(data_complete_norm_24[,], file="Probes_LnBr_BGC1-2_DE_24.txt", sep="\t", eol="\n", quote=FALSE, row.names=FALSE)
write.table(data_complete_norm_60[,], file="Probes_LnBr_BGC1-2_DE_60.txt", sep="\t", eol="\n", quote=FALSE, row.names=FALSE)

#1-G.) OUTPUT ALL PROBES WITH DE >= 2 FOLD 
x_24_2fold = data_complete_norm_24[which(abs(data_complete_norm_24[,"DE_24_bgc"]) > 1 ),]
x_60_2fold = data_complete_norm_60[which(abs(data_complete_norm_60[,"DE_60_bgc"]) > 1 ),]

datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/summary_out"
dir(datadir)
setwd(datadir)
write.table(x_24_2fold, file="DE24_2fold_BGC2.txt", sep="\t", eol="\n", quote=FALSE, row.names=FALSE)
write.table(x_60_2fold, file="DE60_2fold_BGC2.txt", sep="\t", eol="\n", quote=FALSE, row.names=FALSE)


#2-A.) GENE LEVEL DIFFERENTIAL EXPRESSION
#24-MER DATA
x = c(which(data_complete_norm_24[,"Probe_type"]=="Exon"), which(data_complete_norm_24[,"Exons_skipped"] == 0))
exon_canonical_24 = data_complete_norm_24[x,]

#Get a unique list of all gene IDs
genes = unique(exon_canonical_24[,"AlexaGene_ID"])
genes = genes[(which(genes > 0))]

gene_stats_24 = matrix(0, nrow=length(genes), ncol=17, 
                dimnames=list(genes, c("exon_canonical_probes",
                  "mean_bgc_ratio_ln","sd_bgc_ratio_ln","pval_wt_bgc_ratio_ln",
                  "mean_bgc2_ratio_ln","sd_bgc2_ratio_ln","pval_wt_bgc2_ratio_ln",
                  "mean_bgc_ratio_br","sd_bgc_ratio_br","pval_wt_bgc_ratio_br",
                  "mean_bgc2_ratio_br","sd_bgc2_ratio_br","pval_wt_bgc2_ratio_br",
			"mean_DE_bgc_ratio","pval_wt_DE_bgc_ratio","mean_DE_bgc2_ratio","pval_wt_DE_bgc2_ratio")))

for (i in 1:length(genes)){
#for (i in 1:10){
  #First count the number of exon and canonical probes for each gene
  gene_probes = exon_canonical_24[which(exon_canonical_24[,"AlexaGene_ID"]==genes[i]),]
  probe_count = length(gene_probes[,"Probe_ID"])
  gene_stats_24[i,"exon_canonical_probes"] = probe_count

  if (probe_count == 1){
    next()
  }

  #BGC1 - Use reverse-junctions to correct junction probes and NC probes to correct exon probes

  #Now calculate the mean and SD of the background corrected ratio for these probes for LnCAP and Brain
  #Try a Wilcoxon's test (as used by Affy) to evaluate the significance of the ratios
  #i.e. test whether they are significantly greater than 1
  gene_stats_24[i,"mean_bgc_ratio_ln"] = mean(gene_probes[,"LnCAP_24_bgc_ratio"])
  gene_stats_24[i,"sd_bgc_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc_ratio"]))
  gene_stats_24[i,"mean_bgc_ratio_br"] = mean(gene_probes[,"Brain_24_bgc_ratio"])
  gene_stats_24[i,"sd_bgc_ratio_br"] = sqrt(var(gene_probes[,"Brain_24_bgc_ratio"]))

  gene_stats_24[i,"pval_wt_bgc_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_wt_bgc_ratio_br"] = wilcox.test(gene_probes[,"Brain_24_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #BGC2 - Use NC probes to correct ALL probes
  gene_stats_24[i,"mean_bgc2_ratio_ln"] = mean(gene_probes[,"LnCAP_24_bgc2_ratio"])
  gene_stats_24[i,"sd_bgc2_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_24_bgc2_ratio"]))
  gene_stats_24[i,"mean_bgc2_ratio_br"] = mean(gene_probes[,"Brain_24_bgc2_ratio"])
  gene_stats_24[i,"sd_bgc2_ratio_br"] = sqrt(var(gene_probes[,"Brain_24_bgc2_ratio"]))

  gene_stats_24[i,"pval_wt_bgc2_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_24[i,"pval_wt_bgc2_ratio_br"] = wilcox.test(gene_probes[,"Brain_24_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now calculate the DE for each gene using BGC1 and BGC2 values
  #Remember to convert ratios < 1 to 1 before calculating the log2 ratio

  #BGC1
  ln24_bgc = gene_probes[,"LnCAP_24_bgc_ratio"]
  ln24_bgc[which(ln24_bgc < 1)] = 1
  br24_bgc = gene_probes[,"Brain_24_bgc_ratio"]
  br24_bgc[which(br24_bgc < 1)] = 1 
  de_24_bgc = log2(ln24_bgc) - log2(br24_bgc)

  gene_stats_24[i,"mean_DE_bgc_ratio"] = mean(de_24_bgc)
  gene_stats_24[i,"pval_wt_DE_bgc_ratio"] = wilcox.test(de_24_bgc, paired=FALSE, alternative = "two.sided", mu = 0)$p.value

  #BGC2
  ln24_bgc2 = gene_probes[,"LnCAP_24_bgc2_ratio"]
  ln24_bgc2[which(ln24_bgc2 < 1)] = 1
  br24_bgc2 = gene_probes[,"Brain_24_bgc2_ratio"]
  br24_bgc2[which(br24_bgc2 < 1)] = 1 
  de_24_bgc2 = log2(ln24_bgc2) - log2(br24_bgc2)

  gene_stats_24[i,"mean_DE_bgc2_ratio"] = mean(de_24_bgc2)
  gene_stats_24[i,"pval_wt_DE_bgc2_ratio"] = wilcox.test(de_24_bgc2, paired=FALSE, alternative = "two.sided", mu = 0)$p.value

}

#60-MER DATA
x = c(which(data_complete_norm_60[,"Probe_type"]=="Exon"), which(data_complete_norm_60[,"Exons_skipped"] == 0))
exon_canonical_60 = data_complete_norm_60[x,]

#Get a unique list of all gene IDs
genes = unique(exon_canonical_60[,"AlexaGene_ID"])
genes = genes[(which(genes > 0))]

gene_stats_60 = matrix(0, nrow=length(genes), ncol=17, 
                dimnames=list(genes, c("exon_canonical_probes",
                  "mean_bgc_ratio_ln","sd_bgc_ratio_ln","pval_wt_bgc_ratio_ln",
                  "mean_bgc2_ratio_ln","sd_bgc2_ratio_ln","pval_wt_bgc2_ratio_ln",
                  "mean_bgc_ratio_br","sd_bgc_ratio_br","pval_wt_bgc_ratio_br",
                  "mean_bgc2_ratio_br","sd_bgc2_ratio_br","pval_wt_bgc2_ratio_br",
			"mean_DE_bgc_ratio","pval_wt_DE_bgc_ratio","mean_DE_bgc2_ratio","pval_wt_DE_bgc2_ratio")))

for (i in 1:length(genes)){
#for (i in 1:10){
  #First count the number of exon and canonical probes for each gene
  gene_probes = exon_canonical_60[which(exon_canonical_60[,"AlexaGene_ID"]==genes[i]),]
  probe_count = length(gene_probes[,"Probe_ID"])
  gene_stats_60[i,"exon_canonical_probes"] = probe_count

  if (probe_count == 1){
    next()
  }

  #BGC1 - Use reverse-junctions to correct junction probes and NC probes to correct exon probes

  #Now calculate the mean and SD of the background corrected ratio for these probes for LnCAP and Brain
  #Try a Wilcoxon's test (as used by Affy) to evaluate the significance of the ratios
  #i.e. test whether they are significantly greater than 1
  gene_stats_60[i,"mean_bgc_ratio_ln"] = mean(gene_probes[,"LnCAP_60_bgc_ratio"])
  gene_stats_60[i,"sd_bgc_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc_ratio"]))
  gene_stats_60[i,"mean_bgc_ratio_br"] = mean(gene_probes[,"Brain_60_bgc_ratio"])
  gene_stats_60[i,"sd_bgc_ratio_br"] = sqrt(var(gene_probes[,"Brain_60_bgc_ratio"]))

  gene_stats_60[i,"pval_wt_bgc_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_wt_bgc_ratio_br"] = wilcox.test(gene_probes[,"Brain_60_bgc_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #BGC2 - Use NC probes to correct ALL probes
  gene_stats_60[i,"mean_bgc2_ratio_ln"] = mean(gene_probes[,"LnCAP_60_bgc2_ratio"])
  gene_stats_60[i,"sd_bgc2_ratio_ln"] = sqrt(var(gene_probes[,"LnCAP_60_bgc2_ratio"]))
  gene_stats_60[i,"mean_bgc2_ratio_br"] = mean(gene_probes[,"Brain_60_bgc2_ratio"])
  gene_stats_60[i,"sd_bgc2_ratio_br"] = sqrt(var(gene_probes[,"Brain_60_bgc2_ratio"]))

  gene_stats_60[i,"pval_wt_bgc2_ratio_ln"] = wilcox.test(gene_probes[,"LnCAP_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value
  gene_stats_60[i,"pval_wt_bgc2_ratio_br"] = wilcox.test(gene_probes[,"Brain_60_bgc2_ratio"], paired=FALSE, alternative = "greater", mu = 1)$p.value

  #Now calculate the DE for each gene using BGC1 and BGC2 values
  #Remember to convert ratios < 1 to 1 before calculating the log2 ratio

  #BGC1
  ln60_bgc = gene_probes[,"LnCAP_60_bgc_ratio"]
  ln60_bgc[which(ln60_bgc < 1)] = 1
  br60_bgc = gene_probes[,"Brain_60_bgc_ratio"]
  br60_bgc[which(br60_bgc < 1)] = 1 
  de_60_bgc = log2(ln60_bgc) - log2(br60_bgc)

  gene_stats_60[i,"mean_DE_bgc_ratio"] = mean(de_60_bgc)
  gene_stats_60[i,"pval_wt_DE_bgc_ratio"] = wilcox.test(de_60_bgc, paired=FALSE, alternative = "two.sided", mu = 0)$p.value

  #BGC2
  ln60_bgc2 = gene_probes[,"LnCAP_60_bgc2_ratio"]
  ln60_bgc2[which(ln60_bgc2 < 1)] = 1
  br60_bgc2 = gene_probes[,"Brain_60_bgc2_ratio"]
  br60_bgc2[which(br60_bgc2 < 1)] = 1 
  de_60_bgc2 = log2(ln60_bgc2) - log2(br60_bgc2)

  gene_stats_60[i,"mean_DE_bgc2_ratio"] = mean(de_60_bgc2)
  gene_stats_60[i,"pval_wt_DE_bgc2_ratio"] = wilcox.test(de_60_bgc2, paired=FALSE, alternative = "two.sided", mu = 0)$p.value

}

#2-B.) WRITE OUT THE GENE LEVEL SUMMARY STATISTICS
datadir = "C:/Documents and Settings/MG/My Documents/Grad Studies/Project Documents/ArrayAnalysis/Malachi_Nimblegen_pilot_010606/summary_out"
dir(datadir)
setwd(datadir)
write.table(gene_stats_24[,], file="GeneStats_24.txt", sep="\t", eol="\n", row.names=TRUE, quote=FALSE)
write.table(gene_stats_60[,], file="GeneStats_60.txt", sep="\t", eol="\n", row.names=TRUE, quote=FALSE)

#3.) CREATE VIEW OF THE BGC PROBE INTENSITY RATIOS FOR EACH PROBE OF A GENE

#display_gene(data_complete_norm_24, 6782, "LnCAP", "24")
display_gene = function(dataset, gene_id, sample, hybe){

  gene_data = dataset[which(dataset[,"AlexaGene_ID"]==gene_id),]

  exon_probes = gene_data[which(gene_data[,"Probe_type"]=="Exon"),]
  canonical_probes = gene_data[which(gene_data[,"Exons_skipped"] == 0),]
  skip_probes = gene_data[which(gene_data[,"Exons_skipped"] > 0),]
  ei_probes = gene_data[which(gene_data[,"Probe_type"] == "Exon-Intron" | gene_data[,"Probe_type"] == "Intron-Exon") ,]

  data_name = paste(sample,"_",hybe,"_bgc2_ratio",sep="")

  min_x = min(gene_data[,"Tm"])
  max_x = max(gene_data[,"Tm"])
  min_y = min(gene_data[,data_name])
  max_y = max(gene_data[,data_name])

  title = paste(sample,"_",hybe,"- ALEXA GENE:", gene_id, sep="")

  plot(gene_data[,"Tm"], gene_data[,data_name], col="green", xlim=c(min_x,max_x), ylim=c(min_y,max_y), xlab="Tm", ylab="Background Corrected Intensity Ratio", main=title)
  points(canonical_probes[,"Tm"], canonical_probes[,data_name], col="blue")
  points(skip_probes[,"Tm"], skip_probes[,data_name], col="red")
  points(ei_probes[,"Tm"], ei_probes[,data_name], col="black")

  legend(locator(n=1,type="n"),c("Exon","Canonical","Skip","Exon-Intron"), col=c("green","blue","red","black"), pch=c(1,1,1,1))
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


