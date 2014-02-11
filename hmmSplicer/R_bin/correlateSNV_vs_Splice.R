
#Load neccessary libraries
library("gcrma")
library("preprocessCore")
library("Cairo")
library("RSVGTipsDevice")


#Specify input data - Triple Negative Breast Libraries
id_file="/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/SSMutations_Junctions_IDs.tsv"
exp_file="/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/JunctionExpressionMatrix.tsv"
exon_content_file="/projects/alexa2/hmmSplicer/ReferenceAnnotations/hg18/ALL.ExonContent"
junc_coord_file="/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/SSMutations_JunctionCoords.tsv"
outdir="/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/figures/"

#Specify input data - TCGA Libraries
id_file="/projects/alexa2/hmmSplicer/TCGA/summary/spliceSiteMutations/SSMutations_Junctions_IDs.tsv"
exp_file="/projects/alexa2/hmmSplicer/TCGA/summary/unfiltered/JunctionExpressionMatrix.tsv"
exon_content_file="/projects/alexa2/hmmSplicer/ReferenceAnnotations/hg18/ALL.ExonContent"
junc_coord_file="/projects/alexa2/hmmSplicer/TCGA/summary/spliceSiteMutations/SSMutations_JunctionCoords.tsv"
outdir="/projects/alexa2/hmmSplicer/TCGA/summary/spliceSiteMutations/figures/"


#Load the junction expression matrix
exp_data = read.table(file=exp_file, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings="na")
id_data = read.table(file=id_file, header=TRUE, sep="\t", as.is=c(1:6,9:12), na.strings=c("na","NA"))
ec_data = read.table(file=exon_content_file, header=FALSE, sep="\t", as.is=c(1,5), na.strings=c("na","NA"))
names(ec_data)=c("chr","start","end","tcount","strand")
jc_data = read.table(file=junc_coord_file, header=TRUE, sep="\t", as.is=c(1:3), na.strings=c("na","NA"))

#Set rownames to be the JIDs
rownames(exp_data)=exp_data[,"JID"]
rownames(jc_data) = jc_data[,"jid"]

#Change to output dir
setwd(outdir)

#Define global colors
kj_color = "#31A354"
aj_color = "#756BB1"
intron_color = "black"
exon_color = "gray90"
mutant_color = "red"


#Normalize the junction expression matrix
x_norm = as.data.frame(normalize.quantiles(as.matrix(exp_data[,13:(ncol(exp_data)-2)])))
lib_ids = colnames(exp_data[,13:(ncol(exp_data)-2)])
colnames(x_norm)=lib_ids

#Replace the original read counts with the normalized values - then remove the norm array to reduce memory usage
l=length(exp_data[,1])
exp_data[1:l,lib_ids]=x_norm[1:l,]
remove(x_norm)

#Load the splicing model drawing function
source(file="/home/malachig/svn/hmmSplicer/R_bin/splicingModelPlot.R")  

#Run the loop that produced SVG figures for each mutation and stores calculated values in the id_data object
#Generate 'active' version using devSVRTips
active = 1
source(file="/home/malachig/svn/hmmSplicer/R_bin/splicingFigureLoop.R")  

#Generate 'inactive' version with prettier formating using Cairo
active = 0
source(file="/home/malachig/svn/hmmSplicer/R_bin/splicingFigureLoop.R")  


#Assign classes to each mutation event in this order - within each class, rank according to average median diff
#1.) Gain of alternative junction                          -> 'GAIN_ALTER'
#2.) Reciprocal gain of alternative and loss of canonical  -> 'RECIPROCAL'
#3.) Loss of canonical only                                -> 'LOSS_CANON'
#4.) Gain of canonical only                                -> 'GAIN_CANON'
#5.) NA                                                    -> 'UNKNOWN'
id_data[,"class"] = "UNCHANGED"
id_data[which(id_data[,"max_gain_alter"] > 0), "class"] = "GAIN_ALTER"
id_data[which(id_data[,"max_recip_diff"] > 0), "class"] = "RECIPROCAL"
id_data[which(id_data[,"max_loss_canon"] < 0 & id_data[,"class"] == "UNCHANGED"), "class"] = "LOSS_CANON"
id_data[which(id_data[,"max_gain_canon"] > 0 & id_data[,"class"] == "UNCHANGED"), "class"] = "GAIN_CANON"


#write out a new summary results file ranked by cumulative median difference
#o=order(id_data[,"cum_median_diff"], decreasing=TRUE)
o=order(id_data[,"avg_median_diff"], decreasing=TRUE)
write.table(file="RankedList.tsv", id_data[o,], sep="\t", quote=FALSE, row.names=FALSE)

