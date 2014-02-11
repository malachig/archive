#Basic stats and exploration of raw data

#Specify the working directory 
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted"
dir(datadir)
setwd(datadir)

#Now input the actual intensity measurements measurements 
#File contains data for multiple hybridizations

#IMPORT DATA FILE 
#Note: Convert 'na' values to NA.  Use 'as.is' lines with strings

#Complete file
data_complete = read.table("All_hybes_withProbeInfo.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,12,15), na.strings='na')

#CHANGE TO AN OUTPUT DIRECTORY
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/Figures"
dir(datadir)
setwd(datadir)

#Sanity checks:

#Total probes:
length(data_complete[,"Probe_Type"])

#Count probes of each type
table(data_complete[,"Probe_Type"])

#Total genes:
length(unique(data_complete[,"AlexaGene_ID"]))

#Simple box plots of raw data for all probes
par(mfrow=c(3,2))
boxplot(data_complete[,"MIP_EF_A"], main="MIP_EF_A")
boxplot(data_complete[,"FUR_EF_A"], main="FUR_EF_A")
boxplot(data_complete[,"MIP_EF_B"], main="MIP_EF_B")
boxplot(data_complete[,"FUR_EF_B"], main="FUR_EF_B")
boxplot(data_complete[,"MIP_GH_A"], main="MIP_GH_A")
boxplot(data_complete[,"FUR_GH_A"], main="FUR_GH_A")

#Apply a simple ceiling to see the distibution near the majority of points
par(mfrow=c(3,2))
boxplot(data_complete[,"MIP_EF_A"], main="MIP_EF_A", ylim=c(0,2000))
boxplot(data_complete[,"FUR_EF_A"], main="FUR_EF_A", ylim=c(0,2000))
boxplot(data_complete[,"MIP_EF_B"], main="MIP_EF_B", ylim=c(0,2000))
boxplot(data_complete[,"FUR_EF_B"], main="FUR_EF_B", ylim=c(0,2000))
boxplot(data_complete[,"MIP_GH_A"], main="MIP_GH_A", ylim=c(0,2000))
boxplot(data_complete[,"FUR_GH_A"], main="FUR_GH_A", ylim=c(0,2000))

# CREATE PLOTS OF ALL RAW INTENSITY VALUES FOR EACH EXPERIMENT

#Create plots of all raw intensities - ORDERED
mip_ef_a_sorted = data_complete[sort.list(data_complete[,"MIP_EF_A"]),"MIP_EF_A"]
mip_ef_b_sorted = data_complete[sort.list(data_complete[,"MIP_EF_B"]),"MIP_EF_B"]
mip_gh_a_sorted = data_complete[sort.list(data_complete[,"MIP_GH_A"]),"MIP_GH_A"]
fur_ef_a_sorted = data_complete[sort.list(data_complete[,"FUR_EF_A"]),"FUR_EF_A"]
fur_ef_b_sorted = data_complete[sort.list(data_complete[,"FUR_EF_B"]),"FUR_EF_B"]
fur_gh_a_sorted = data_complete[sort.list(data_complete[,"FUR_GH_A"]),"FUR_GH_A"]

probes = matrix(0,nrow=385000,ncol=6)
probes[,1:6] = 1:385000

intensities = matrix(0,nrow=385000,ncol=6, dimnames=list(1:385000, c("mip_ef_a","mip_ef_b","mip_gh_a","fur_ef_a","fur_ef_b","fur_gh_a")))
intensities[,1:6] = c(mip_ef_a_sorted,mip_ef_b_sorted,mip_gh_a_sorted,fur_ef_a_sorted,fur_ef_b_sorted,fur_gh_a_sorted)

#Plot all of the ordered intensities on the same graph 
matplot(probes, intensities, type="l", lty=1:6, col=1:6, lwd=2, xlab="probe rank", ylab="raw intensity", main="Ordered Raw Intensities")
legend(locator(n=1,type="n"), c("mip_ef_a","mip_ef_b","mip_gh_a","fur_ef_a","fur_ef_b","fur_gh_a"), lty=1:6, col=1:6, lwd=2)

matplot(probes, log2(intensities), type="l", lty=1:6, col=1:6, lwd=2, xlab="probe rank", ylab="log2 raw intensity", main="Ordered Raw Intensities")
legend(locator(n=1,type="n"), c("mip_ef_a","mip_ef_b","mip_gh_a","fur_ef_a","fur_ef_b","fur_gh_a"), lty=1:6, col=1:6, lwd=2)


# SUMMARIZE INTENSITY VALUES FOR EACH OF THE VARIOUS PROBE TYPES
exon_probes = data_complete[which(data_complete[,"Probe_Type"]=="Exon"),]
canonical_probes = data_complete[which(data_complete[,"Exons_Skipped"]==0),]
skip1_probes = data_complete[which(data_complete[,"Exons_Skipped"]==1),]
skip2_probes = data_complete[which(data_complete[,"Exons_Skipped"]==2),]
skip3_probes = data_complete[which(data_complete[,"Exons_Skipped"]==3),]
exonIntron_probes = data_complete[which(data_complete[,"Probe_Type"]=="Exon-Intron"),]
intronExon_probes = data_complete[which(data_complete[,"Probe_Type"]=="Intron-Exon"),]
intron_probes = data_complete[which(data_complete[,"Probe_Type"]=="Intron"),]
control_probes = data_complete[which(data_complete[,"Probe_Type"]=="Control-Negative"),]

createlist = function(exp){
  probe_list = list (exon_probes[,exp], canonical_probes[,exp], skip1_probes[,exp],
                     skip2_probes[,exp], skip3_probes[,exp], exonIntron_probes[,exp], intronExon_probes[,exp],
	               intron_probes[,exp], control_probes[,exp])
  names(probe_list) = c("E","S0","S1","S2","S3","EI","IE","I","NC")
  return(probe_list)
}

#Create a plot of raw intensities, one line for each probe type, probes ordered within each type
x = createlist("MIP_EF_A")
y = lapply(x,order)
ordered_probelist = list (x$NC[y$NC], x$I[y$I], x$EI[y$EI], x$IE[y$IE], x$S3[y$S3],
                          x$S2[y$S2], x$S1[y$S1], x$S0[y$S0], x$E[y$E]) 
names(ordered_probelist) = c("NC","I","EI","IE","S3","S2","S1","S0","E")
levels(ordered_probelist) = c("NC","I","EI","IE","S3","S2","S1","S0","E")

plot (1:length(ordered_probelist$E), ordered_probelist$E, type="n", xlab="probe rank", ylab="raw intensity", main="MIP_EF_A - Ordered Raw Intensities (by probe type)") 
points(1:length(ordered_probelist$NC), ordered_probelist$NC, col=1, pch=1)
points(1:length(ordered_probelist$I), ordered_probelist$I, col=2, pch=2)
points(1:length(ordered_probelist$EI), ordered_probelist$EI, col=3, pch=3)
points(1:length(ordered_probelist$IE), ordered_probelist$IE, col=4, pch=4)
points(1:length(ordered_probelist$S3), ordered_probelist$S3, col=5, pch=5)
points(1:length(ordered_probelist$S2), ordered_probelist$S2, col=6, pch=6)
points(1:length(ordered_probelist$S1), ordered_probelist$S1, col=7, pch=7)
points(1:length(ordered_probelist$S0), ordered_probelist$S0, col=8, pch=8)
points(1:length(ordered_probelist$E), ordered_probelist$E, col=9, pch=9)
legend(locator(n=1,type="n"),c("NC","I","EI","IE","S3","S2","S1","S0","E"), col=1:9, pch=1:9)

#Boxplots by probe type
boxplot(createlist("MIP_EF_A"), main="MIP_EF_A - Raw intensities by probe type")
boxplot(createlist("MIP_EF_B"), main="MIP_EF_B - Raw intensities by probe type")
boxplot(createlist("MIP_GH_A"), main="MIP_GH_A - Raw intensities by probe type")
boxplot(createlist("FUR_EF_A"), main="FUR_EF_A - Raw intensities by probe type")
boxplot(createlist("FUR_EF_B"), main="FUR_EF_B - Raw intensities by probe type")
boxplot(createlist("FUR_GH_A"), main="FUR_GH_A - Raw intensities by probe type")

#Try zooming into the area of interest by cutting off the outliers on the graph
boxplot(createlist("MIP_EF_A"), main="MIP_EF_A - Raw intensities by probe type", ylim=c(0,10000))
boxplot(createlist("MIP_EF_B"), main="MIP_EF_B - Raw intensities by probe type", ylim=c(0,10000))
boxplot(createlist("MIP_GH_A"), main="MIP_GH_A - Raw intensities by probe type", ylim=c(0,10000))
boxplot(createlist("FUR_EF_A"), main="FUR_EF_A - Raw intensities by probe type", ylim=c(0,10000))
boxplot(createlist("FUR_EF_B"), main="FUR_EF_B - Raw intensities by probe type", ylim=c(0,10000))
boxplot(createlist("FUR_GH_A"), main="FUR_GH_A - Raw intensities by probe type", ylim=c(0,10000))


#Summarize some basic stats for each probe type for each experiment
summarize_probe_means = function(exp){
  probe_means = c(mean(exon_probes[,exp]), mean(canonical_probes[,exp]),
                  mean(skip1_probes[,exp]), mean(skip2_probes[,exp]), mean(skip3_probes[,exp]),
		      mean(exonIntron_probes[,exp]), mean(intronExon_probes[,exp]), mean(intron_probes[,exp]),
	            mean(control_probes[,exp]))
  names(probe_means) = c("E","S0","S1","S2","S3","EI","IE","I","NC")
  return(probe_means)
}
summarize_probe_means("MIP_EF_A")

summarize_probe_medians = function(exp){
  probe_medians = c(median(exon_probes[,exp]), median(canonical_probes[,exp]),
                  median(skip1_probes[,exp]), median(skip2_probes[,exp]), median(skip3_probes[,exp]),
		      median(exonIntron_probes[,exp]), median(intronExon_probes[,exp]), median(intron_probes[,exp]),
	            median(control_probes[,exp]))
  names(probe_medians) = c("E","S0","S1","S2","S3","EI","IE","I","NC")
  return(probe_medians)
}
summarize_probe_medians("MIP_EF_A")

#COMPARE THE VARIABILITY WITHIN PROBESETS COMPARED TO THE VARIABILITY ACROSS PROBESETS FOR ALL EXON PROBES

#COMPARE THE RAW INTENSITIES FOR NEGATIVE CONTROL PROBES TO THOSE OF EXON PROBES








