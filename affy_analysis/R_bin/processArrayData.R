#Written by Malachi Griffith

#This script takes data files from Expression Console (already normalized, background corrected, and summarized to the gene or exon level)
#It calculates the average of biological replicates and the differential expression values for a pairwise comparison

#GENE-LEVEL DATA - MIP101/MIP5FUR only data

gene_datadir = "/projects/malachig/affymetrix_exon_arrays/Data/MIP101_and_MIP5FU_only/EnsEMBL_v49_gene-level"
setwd(gene_datadir)
dir()

#Import gene data
gene_data = read.table("iterplier-gene-other.summary.txt", header = TRUE, sep = "\t", comment.char = "#", nrows=27244)

#Add an arbitrary small value to raw data before taking log to stabalize variance
#var_stab_adjust = 16
#var_stab_adjust = 4
var_stab_adjust = 1

#Rename columns to something more convenient
colnames(gene_data) = c("alexa_gene_id","MIP101_Rep1","MIP5FUR_Rep1","MIP101_Rep2","MIP5FUR_Rep2","MIP101_Rep3","MIP5FUR_Rep3","MIP101_Rep6","MIP5FUR_Rep6","MIP101_Rep2b")

#Take the average of biological replicates - dont use polyA+ samples
gene_data[,"MIP101_mean"] = apply (gene_data[,c("MIP101_Rep1", "MIP101_Rep2", "MIP101_Rep3", "MIP101_Rep2b")], 1, mean)
gene_data[,"MIP5FUR_mean"] = apply (gene_data[,c("MIP5FUR_Rep1", "MIP5FUR_Rep2", "MIP5FUR_Rep3")], 1, mean)

#Convert means to log2 scale
gene_data[,"MIP101_mean_log2"] = log2(gene_data[,"MIP101_mean"] + var_stab_adjust)
gene_data[,"MIP5FUR_mean_log2"] = log2(gene_data[,"MIP5FUR_mean"] + var_stab_adjust)

#Calculate log2 DE (MIP101 - MIP5FUR)
gene_data[,"mean_log2_DE"] = gene_data[,"MIP101_mean_log2"] - gene_data[,"MIP5FUR_mean_log2"]


#write.table(gene_data[,c("alexa_gene_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="GeneMeans_PLIER_plus16_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)

#write.table(gene_data[,c("alexa_gene_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="GeneMeans_PLIER_plus4_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)

write.table(gene_data[,c("alexa_gene_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="GeneMeans_PLIER_plus1_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)



#EXON-LEVEL DATA - MIP101/MIP5FUR only data
exon_datadir = "/projects/malachig/affymetrix_exon_arrays/Data/MIP101_and_MIP5FU_only/exon-level"
setwd(exon_datadir)
dir()

exon_data = read.table("iterplier-exon-all.summary.txt", header = TRUE, sep = "\t", comment.char = "#", nrows=1412447)

#Add an arbitrary small value to raw data before taking log to stabalize variance
#var_stab_adjust = 16
#var_stab_adjust = 4
var_stab_adjust = 1

#Rename columns to something more convenient
colnames(exon_data) = c("probeset_id","MIP101_Rep1","MIP5FUR_Rep1","MIP101_Rep2","MIP5FUR_Rep2","MIP101_Rep3","MIP5FUR_Rep3","MIP101_Rep6","MIP5FUR_Rep6","MIP101_Rep2b")

#Take the average of biological replicates - dont use polyA+ samples
exon_data[,"MIP101_mean"] = apply (exon_data[,c("MIP101_Rep1", "MIP101_Rep2", "MIP101_Rep3", "MIP101_Rep2b")], 1, mean)
exon_data[,"MIP5FUR_mean"] = apply (exon_data[,c("MIP5FUR_Rep1", "MIP5FUR_Rep2", "MIP5FUR_Rep3")], 1, mean)

#Convert means to log2 scale
exon_data[,"MIP101_mean_log2"] = log2(exon_data[,"MIP101_mean"] + var_stab_adjust)
exon_data[,"MIP5FUR_mean_log2"] = log2(exon_data[,"MIP5FUR_mean"] + var_stab_adjust)

#Calculate log2 DE (MIP101 - MIP5FUR)
exon_data[,"mean_log2_DE"] = exon_data[,"MIP101_mean_log2"] - exon_data[,"MIP5FUR_mean_log2"]


#write.table(exon_data[,c("probeset_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="ExonMeans_PLIER_plus16_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)

#write.table(exon_data[,c("probeset_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="ExonMeans_PLIER_plus4_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)

write.table(exon_data[,c("probeset_id","MIP101_mean_log2","MIP5FUR_mean_log2","mean_log2_DE")], file="ExonMeans_PLIER_plus1_Log2_DE.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)




#EXON-LEVEL DATA - ALL 45 EXPERIMENTS
exon_datadir = "/projects/malachig/affymetrix_exon_arrays/Data/Probeset-level_45Experiments"
setwd(exon_datadir)
dir()

exon_data = read.table("iterplier_all_exons_summary.txt", header = TRUE, sep = "\t", comment.char = "#", nrows=1412482)

#Rename columns to something more convenient
colnames(exon_data) = c("probeset_id","MIP101_Rep1","MIP5FUR_Rep1","MIP101_Rep2","MIP5FUR_Rep2","MIP101_Rep3","MIP5FUR_Rep3",
			"HEK_CrkRS","HEK_CrkRS_DN","HEK_CrkRS_Vec","LnCAP_AR_Rep1","LnCAP_EtOH_Rep1","LnCAP_AR_Rep2","LnCAP_EtOH_Rep2",
			"Lymphoma_1","Lymphoma_2","Lymphoma_3","Lymphoma_4","Lymphoma_5","Lymphoma_6","HL60_Rep1","HL60_Rep2","HL60_Rep3",
		        "MIP101_Rep6","MIP5FUR_Rep6","RKO_Rep1","RKO5FUR_Rep1","MIP101_Rep2b","RKO_Rep1","RKO_Rep2","RKO_Rep3",
		        "RKO5FUR_Rep1","RKO5FUR_Rep2","RKO5FUR_Rep3","HCT_Rep1","HCT_Rep2","HCT_Rep3","HCT5FUR_Rep1","HCT5FUR_Rep2","HCT5FUR_Rep3",
			"MIP101_plus5FU_Rep1","MIP101_plus5FU_Rep2","MIP101_plus5FU_Rep3","MIP5FUR_plus5FU_Rep1","MIP5FUR_plus5FU_Rep2",
			"MIP5FUR_plus5FU_Rep3")

exon_data[,"MIP101_mean"] = log2((apply (exon_data[,c("MIP101_Rep1", "MIP101_Rep2", "MIP101_Rep3", "MIP101_Rep2b")], 1, mean))+16)

exon_data_means = exon_data[,c("probeset_id","MIP101_mean")]
exon_data_means[,"MIP101_polyA_Rep6"] = log2(exon_data[,"MIP101_Rep6"] + 16)
exon_data_means[,"MIP5FUR_mean"] = log2((apply (exon_data[,c("MIP5FUR_Rep1", "MIP5FUR_Rep2", "MIP5FUR_Rep3")], 1, mean))+16)
exon_data_means[,"MIP5FUR_polyA_Rep6"] = log2(exon_data[,"MIP5FUR_Rep6"] + 16)
exon_data_means[,"RKO_mean"] = log2((apply (exon_data[,c("RKO_Rep1", "RKO_Rep2", "RKO_Rep3")], 1, mean))+16)
exon_data_means[,"RKO5FUR_mean"] = log2((apply (exon_data[,c("RKO5FUR_Rep1", "RKO5FUR_Rep2", "RKO5FUR_Rep3")], 1, mean))+16)
exon_data_means[,"HCT_mean"] = log2((apply (exon_data[,c("HCT_Rep1", "HCT_Rep2", "HCT_Rep3")], 1, mean))+16)
exon_data_means[,"HCT5FUR_mean"] = log2((apply (exon_data[,c("HCT5FUR_Rep1", "HCT5FUR_Rep2", "HCT5FUR_Rep3")], 1, mean))+16)
exon_data_means[,"MIP101_plus5FU"] = log2((apply (exon_data[,c("MIP101_plus5FU_Rep1", "MIP101_plus5FU_Rep2", "MIP101_plus5FU_Rep3")], 1, mean))+16)
exon_data_means[,"MIP5FUR_plus5FU"] = log2((apply (exon_data[,c("MIP5FUR_plus5FU_Rep1", "MIP5FUR_plus5FU_Rep2", "MIP5FUR_plus5FU_Rep3")], 1, mean))+16)
exon_data_means[,"LnCAP_AR_mean"] = log2((apply (exon_data[,c("LnCAP_AR_Rep1", "LnCAP_AR_Rep2")], 1, mean))+16)
exon_data_means[,"LnCAP_EtOH_mean"] = log2((apply (exon_data[,c("LnCAP_EtOH_Rep1", "LnCAP_EtOH_Rep2")], 1, mean))+16)
exon_data_means[,"HL60_mean"] =  log2((apply (exon_data[,c("HL60_Rep1", "HL60_Rep2","HL60_Rep3")], 1, mean))+16)
exon_data_means[,"Lymphoma_1"] = log2(exon_data[,"Lymphoma_1"] + 16)
exon_data_means[,"Lymphoma_2"] = log2(exon_data[,"Lymphoma_2"] + 16)
exon_data_means[,"Lymphoma_3"] = log2(exon_data[,"Lymphoma_3"] + 16)
exon_data_means[,"Lymphoma_4"] = log2(exon_data[,"Lymphoma_4"] + 16)
exon_data_means[,"Lymphoma_5"] = log2(exon_data[,"Lymphoma_5"] + 16)
exon_data_means[,"Lymphoma_6"] = log2(exon_data[,"Lymphoma_6"] + 16)
exon_data_means[,"HEK_CrkRS"] = log2(exon_data[,"HEK_CrkRS"] + 16)
exon_data_means[,"HEK_CrkRS_DN"] = log2(exon_data[,"HEK_CrkRS_DN"] + 16)
exon_data_means[,"HEK_CrkRS_Vec"] = log2(exon_data[,"HEK_CrkRS_Vec"] + 16)

write.table(exon_data_means[,], file="ExonMeans.txt", sep="\t", eol="\n", row.names=FALSE, quote=FALSE)




 
