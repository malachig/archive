#Draft analysis of compendium junction expressoin matrix created by hmmSplicer pipeline

#Load neccessary libraries
library("gcrma")
library("preprocessCore")
library("Cairo")
library("RSVGTipsDevice")


#Specify input data - Triple Negative Breast Libraries
exp_file="/projects/alexa2/hmmSplicer/Summary/Compendium/JunctionExpressionMatrix.tsv"
lib_file="/projects/alexa2/hmmSplicer/Summary/Compendium/LibraryClasses.tsv"
sum_file="/projects/alexa2/hmmSplicer/Summary/Compendium/LibrarySummary.tsv"
outdir="/projects/alexa2/hmmSplicer/Summary/Compendium/figures/"

#Load the junction expression matrix
exp_data = read.table(file=exp_file, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings=c("na", "NA", "n/a", "N/A"))
lib_data = read.table(file=lib_file, header=TRUE, sep="\t", as.is=c(1:8), na.strings=c("na", "NA", "n/a", "N/A"))
sum_data = read.table(file=sum_file, header=TRUE, sep="\t", as.is=c(1,27), na.strings=c("na", "NA", "n/a", "N/A"))

#Set column names
names(lib_data) = c("library_id", "library_name", "group_name", "tissue", "disease_status", "tissue_type", "subtype", "description")

#Set rownames to be the JIDs
rownames(exp_data)=exp_data[,"JID"]

#Change to output dir
setwd(outdir)

#Normalize the junction expression matrix
x_norm = as.data.frame(normalize.quantiles(as.matrix(exp_data[,13:(ncol(exp_data)-2)])))
lib_ids = colnames(exp_data[,13:(ncol(exp_data)-2)])
colnames(x_norm)=lib_ids

#Replace the original read counts with the normalized values - then remove the norm array to reduce memory usage
l=length(exp_data[,1])
exp_data[1:l,lib_ids]=x_norm[1:l,]
remove(x_norm)

#Perform a garbage collection
gc(verbose=TRUE)

#Create figures and stats to summarize the distribution of libraries (tumour/normal, tissue, sample type, etc.)

#number of libraries total - 599
length(lib_ids)

#Piechart of tissue
pdf(file="TissueTypes_Pie.pdf")
pie(table(lib_data[,"tissue"]), main="Tissues")
dev.off()

#Piechart of disease states
pdf(file="DiseaseStatus_Pie.pdf")
pie(table(lib_data[,"disease_status"]), main="Disease Status")
dev.off()

#Piechart of tissue types
pdf(file="TissueTypes_Pie.pdf")
pie(table(lib_data[,"tissue_type"]), main="Tissue Type")
dev.off()


##############################################################################################################
#Try processing the input data with DESeq
library(DESeq)

#define conditions to be used
group_name_conds = as.factor(lib_data[,"group_name"])
tissue_conds = as.factor(lib_data[,"tissue"])
disease_status_conds = lib_data[,"disease_status"]
tissue_type_conds = as.factor(lib_data[,"tissue_type"])
subtype_conds = as.factor(lib_data[,"subtype"])

#Load the raw counts and do NOT normalize them
exp_data = read.table(file=exp_file, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings=c("na", "NA", "n/a", "N/A"))
rownames(exp_data)=exp_data[,"JID"]
gc(verbose=TRUE)

#Get some simple stats

#Total junction mapping reads = 4,319,439,488
sum(as.numeric(exp_data[,"Total"]))

#Total known junctions detected ('DA') = 223,714
table(exp_data[,"Anchored"])

#Total known junction possible = 314,935

#Percent of all known junctions detected = 71.0%
(223714/314935)*100

#How many high quality novel junctions were discovered of each type (NDA, D, A, or N)
min_score=1000
min_lib_count=5
min_read_count=100
t=length(which(exp_data[,"Score"] > min_score & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
tda=length(which(exp_data[,"Score"] > min_score & exp_data[,"Anchored"] == "DA" & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
tnda=length(which(exp_data[,"Score"] > min_score & exp_data[,"Anchored"] == "NDA" & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
td=length(which(exp_data[,"Score"] > min_score & exp_data[,"Anchored"] == "D" & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
ta=length(which(exp_data[,"Score"] > min_score & exp_data[,"Anchored"] == "A" & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
tn=length(which(exp_data[,"Score"] > min_score & exp_data[,"Anchored"] == "N" & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count))
(tda/t)*100
(tnda/t)*100
(td/t)*100
(ta/t)*100
(tn/t)*100


#Create a sample counts table
countsTable=exp_data[1:10000,13:(ncol(exp_data)-2)]

#Try just 'high' quality values
i=which(exp_data[,"Score"] > min_score & exp_data[,"Library_Count"] > min_lib_count & exp_data[,"Read_Count"] > min_read_count)

countsTable=exp_data[i,13:(ncol(exp_data)-2)]


#Example DE analysis
cds = newCountDataSet(countsTable, disease_status_conds)
cds = estimateSizeFactors(cds)
cds = estimateVarianceFunctions(cds)
res = nbinomTest(cds, "Normal", "Cancer")
gc(verbose=TRUE)
res_backup=res
gc(verbose=TRUE)


#Pull out the IDs for all significant events
sig = which(res$padj < 0.05 & res$log2FoldChange > 5 | res$log2FoldChange < -5)
resSig = res[sig, ]
resSigId=resSig$id
exp_data_sig = exp_data[resSigId,]

#Set the rownames to allow easier lookup
rownames(resSig)=(resSig$id)



#Subset that are NDA, D and A
nda_sig_i=which(exp_data_sig[,"Anchored"] == "NDA")
length(nda_sig_i) #318
nda_sig_jid=exp_data_sig[nda_sig_i,"JID"]
resSig[nda_sig_jid,]$sig = exp_data_sig[nda_sig_i,"Gene_Name"]
resSigNDA = resSig[nda_sig_jid,]
resSigNDA=resSigNDA[which(resSigNDA$baseMean > 1 & resSigNDA$log2FoldChange > 3),]
head(resSigNDA[order(resSigNDA$pval),], n=20)

d_sig_i=which(exp_data_sig[,"Anchored"] == "D")
length(d_sig_i) #326
d_sig_jid=exp_data_sig[d_sig_i, "JID"]
resSig[d_sig_jid,]$sig = exp_data_sig[d_sig_i,"Gene_Name"]
resSigD = resSig[d_sig_jid,]
resSigD=resSigD[which(resSigD$baseMean > 1 & resSigD$log2FoldChange > 3),]
head(resSigD[order(resSigD$pval),], n=20)


a_sig_i=which(exp_data_sig[,"Anchored"] == "A")
length(a_sig_i) #538
a_sig_jid=exp_data_sig[a_sig_i,"JID"]
resSig[a_sig_jid,]$sig = exp_data_sig[a_sig_i,"Gene_Name"]
resSigA = resSig[a_sig_jid,]
resSigA=resSigA[which(resSigA$baseMean > 0.5 & resSigA$log2FoldChange > 10),]
head(resSigA[order(resSigA$pval),], n=20)



checkCountsTable = function(jid){
  normal_i = which(disease_status_conds == "Normal")
  cancer_i = which(disease_status_conds == "Cancer")
  countsTable[jid, normal_i]
  ncount = length(which(countsTable[jid, normal_i] > 0))
  ccount = length(which(countsTable[jid, cancer_i] > 0))
  paste("Normal library count = ", ncount, " | ", "Cancer library count = ", ccount, sep="")
}

#View each list
exp_data_sig[nda_sig_i,"Gene_Name"]


head(resSig[order(abs(resSig$log2FoldChange)),])




#Create a PDF of the differentially expressed features and the subset of those that are novel AS events
#Function to plot DE values
pdf(file="CancerVsNormal_DE_Scatter4.pdf")
sig = which(res$padj < 0.05 & (res$log2FoldChange > 3 | res$log2FoldChange < -3))
sig = which(res$padj < 0.01)

length(sig)

res$sig=0
res[sig,]$sig=1

plotDE = function(res){
	plot (res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.1, col=ifelse(res$sig > 0, "blue", "grey50"), xlab="Junction expression level", ylab="Fold change (Cancer vs. Normal)", cex.lab=1.3, cex.axis=1.2)
}
plotDE(res)

resSig = res[sig, ]

#Grab the significant JIDs 
jids=resSig$id

#Now grab a slice of the original array to tie the DE events back to genes, etc.
temp = exp_data[jids, 1:12]
x=temp[which(temp[,"Exons_Skipped"] > 0 & temp[,"Anchored"] == "DA"),"JID"]
xx = which(res$id %in% x)
length(x)
points (res$baseMean[xx], res$log2FoldChange[xx], pch=20, cex=0.2, col="green")

x=temp[which(temp[,"Anchored"] == "NDA" | temp[,"Anchored"] == "D" | temp[,"Anchored"] == "D"),"JID"]
length(x)
xx = which(res$id %in% x)
points (res$baseMean[xx], res$log2FoldChange[xx], pch=20, cex=0.2, col="red")

legend_text=c("Unchanged", "DE", "DE & known splicing", "DE & novel splicing")
legend_cols=c("grey50","blue","green","red")
legend("topleft", legend_text, col=legend_cols, pch=20, bg="white")
dev.off()




