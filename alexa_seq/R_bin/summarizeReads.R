#!/home/malachig/R64/R-2.9.0/bin/Rscript
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes a tab delimited file of read record data and summarizes basic stats
#These data represent paired end reads that have been mapped to an EnsEMBL transcriptome by BLAST with a word size of 11 and pre-filtered on bit score

#0.) The # of reads
#1.) The # of 'Top_Hit' reads >= length_cutoff bp in length.  For R1 and R2
#2.) The # of 'Top_Hit' reads >= length_cutoff bp in length where BOTH R1 and R2 meet this criteria
#3.) The # of 'Ambiguous' reads, and # of 'None' reads (no hits at all to the transcriptome)
#4.) The distribution of alignment lengths (for hits >= length_cutoff bp in length)
#5.) The distribution of percent position in the transcript for all quality reads (Top_Hit and >= length_cutoff bp)
#    - Do this for R1, R2 and R1 & R2 combined
#    - Do this for small (<1000 bp), medium (1000-5000 bp), large (5000-10000 bp) and very large (>10000) transcripts
#    - Make note of the # of reads occuring in each of these categories
#6.) The distribution of distances between R1 and R2 (again only for top_hits of length > length_cutoff)
#7.) The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to the SAME GENE
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to DIFFERENT GENES
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to the SAME CHROMOSOME
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to DIFFERENT CHROMOSOMES

args = (commandArgs(TRUE))
datafile = args[1]
results_dir = args[2]
record_count = args[3]


#Open input file
datafile = "/projects/malachig/solexa/read_records/HS04391/ENST_v53/HS04391_Lanes1-23_ENST_v53_SELECTED_COLUMNS.txt"
results_dir = ""

dir(datadir)
setwd(datadir)
read_data_hs04391 = read.table("HS04391_Lanes1-23_ENST_v49_SELECTED_COLUMNS.txt", header=T, quote="",sep="\t", comment.char="", as.is=c(1,2,3,4,9), na.strings='NA', nrows=127186727)

#Change to output directory for storage of figures
outdir = "/projects/malachig/solexa/read_records/HS04391/ENST_v49/figures"
setwd(outdir)


#Consider any unambiguous 'Top_Hit' of length 21 or greater to be 'real'.  
bitscore_cutoff=48.1

##############################################################################################################
#0.) The # of reads.  Total mapped paired reads = 127186726 (where one or both reads was mapped)
length(read_data_hs04391[,1])
total_reads = length(read_data_hs04391[,1])


##############################################################################################################
#1a.) The # of 'Top_Hit' reads (unambiguous mappings, any length of hit)  For R1, R2

#R1 = 2329467 (73.2%)
length(which(read_data_hs04391[,"R1_HitType"]=="Top_Hit"))
(length(which(read_data_hs04391[,"R1_HitType"]=="Top_Hit")) / total_reads)*100

#R2 = 2257741 (70.9%)
length(which(read_data_hs04391[,"R2_HitType"]=="Top_Hit"))
(length(which(read_data_hs04391[,"R2_HitType"]=="Top_Hit")) / total_reads)*100


##############################################################################################################
#1b.) The # of 'Top_Hit' reads >= length_cutoff bp in length.  For R1, R2

#R1 = 1824349 (57.3%)
length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff))
(length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff)) / total_reads)*100

#R2 = 1714088 (53.9%)
length(which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff))
(length(which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff)) / total_reads)*100


##############################################################################################################
#2b.) The # of 'Top_Hit' reads where BOTH R1 and R2 meet this criteria

#BOTH R1 and R2 = 1940836 (61.0%)
length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R2_HitType"]=="Top_Hit"))
(length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R2_HitType"]=="Top_Hit")) / total_reads)*100


##############################################################################################################
#2b.) The # of 'Top_Hit' reads >= length_cutoff bp in length where BOTH R1 and R2 meet this criteria

#BOTH R1 and R2 = 1547322 (48.6%)
length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff))
(length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff)) / total_reads)*100


##############################################################################################################
#3.) The # of 'Ambiguous' reads, and # of 'None' reads (no hits at all to the transcriptome)

#R1.  654097 (20.6%) Ambiguous.  199039 (6.3%) No hits.  505118 (15.9%) Top_Hit but less than length_cutoff bp.
#R2.  707866 (22.2%) Ambiguous.  216996 (6.8%) No hits . 543653 (17.1%) Top_Hit but less than length_cutoff bp.

length(which(read_data[,"R1_HitType"]=="Ambiguous"))
(length(which(read_data[,"R1_HitType"]=="Ambiguous")) / total_reads)*100

length(which(read_data[,"R1_HitType"]=="None"))
(length(which(read_data[,"R1_HitType"]=="None")) / total_reads)*100

length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]<length_cutoff))
(length(which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]<length_cutoff))/total_reads)*100


length(which(read_data[,"R2_HitType"]=="Ambiguous"))
(length(which(read_data[,"R2_HitType"]=="Ambiguous")) / total_reads)*100

length(which(read_data[,"R2_HitType"]=="None"))
(length(which(read_data[,"R2_HitType"]=="None")) / total_reads)*100

length(which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]<length_cutoff))
(length(which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]<length_cutoff))/total_reads)*100


##############################################################################################################
#4b.) The distribution of alignment lengths (for top_hits, any length)
temp_I_R1 = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=1)
temp_val_R1 = read_data[temp_I_R1, "R1_AlignmentLength"]

temp_I_R2 = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=1)
temp_val_R2 = read_data[temp_I_R2, "R2_AlignmentLength"]

temp_val_R12 = c(temp_val_R1,temp_val_R2);
count = length(temp_val_R12)
y_title = paste("Read Count (Total = ", count, ")", sep="") 

pdf(file = "ReadLengthsHist_TopHits_anyLength.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
hist(temp_val_R12, col="dark green", xlab="Length of BLAST hit", ylab=y_title, main="Distribution of read hits lengths (any length, R1+R2)", breaks=100)
dev.off()

#4b.) The distribution of alignment lengths (for top_hits >= length_cutoff bp in length)
temp_I_R1 = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff)
temp_val_R1 = read_data[temp_I_R1, "R1_AlignmentLength"]

temp_I_R2 = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff)
temp_val_R2 = read_data[temp_I_R2, "R2_AlignmentLength"]

temp_val_R12 = c(temp_val_R1,temp_val_R2);
count = length(temp_val_R12)
y_title = paste("Read Count (Total = ", count, ")", sep="") 

pdf(file = "ReadLengthsHist_TopHits_21plus.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
hist(temp_val_R12, col="dark green", xlab="Length of BLAST hit", ylab=y_title, main="Distribution of read hits lengths (>= 21 bp, R1+R2)", breaks=100)
dev.off()


#5.) The distribution of percent position in the transcript for all quality reads (Top_Hit and >= length_cutoff bp)
#    - Do this for R1, R2 and R1 & R2 combined?
#    - Do this for small (<1000 bp), medium (1000-5000 bp), large (5000-10000 bp), and very large (>10000 bp) transcripts
#    - Make note of the # of reads occuring in each of these categories
R1_smallTrans_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R1_TranscriptSize"]<1000)
R1_smallTrans_val = read_data[R1_smallTrans_I,"R1_RelativePosition"]

R1_mediumTrans_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R1_TranscriptSize"]>=1000 & read_data[,"R1_TranscriptSize"]<=5000)
R1_mediumTrans_val = read_data[R1_mediumTrans_I,"R1_RelativePosition"]

R1_largeTrans_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R1_TranscriptSize"]>5000 & read_data[,"R1_TranscriptSize"]<=10000)
R1_largeTrans_val = read_data[R1_largeTrans_I,"R1_RelativePosition"]

R1_veryLargeTrans_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R1_TranscriptSize"]>10000)
R1_veryLargeTrans_val = read_data[R1_veryLargeTrans_I,"R1_RelativePosition"]


R2_smallTrans_I = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff & read_data[,"R2_TranscriptSize"]<1000)
R2_smallTrans_val = read_data[R2_smallTrans_I,"R2_RelativePosition"]

R2_mediumTrans_I = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff & read_data[,"R2_TranscriptSize"]>=1000 & read_data[,"R2_TranscriptSize"]<=5000)
R2_mediumTrans_val = read_data[R2_mediumTrans_I,"R2_RelativePosition"]

R2_largeTrans_I = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff & read_data[,"R2_TranscriptSize"]>5000 & read_data[,"R2_TranscriptSize"]<=10000)
R2_largeTrans_val = read_data[R2_largeTrans_I,"R2_RelativePosition"]

R2_veryLargeTrans_I = which(read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff & read_data[,"R2_TranscriptSize"]>10000)
R2_veryLargeTrans_val = read_data[R2_veryLargeTrans_I,"R2_RelativePosition"]


R12_small_positions = c(R1_smallTrans_val,R2_smallTrans_val)
R12_medium_positions = c(R1_mediumTrans_val,R2_mediumTrans_val)
R12_large_positions = c(R1_largeTrans_val,R2_largeTrans_val)
R12_veryLarge_positions = c(R1_veryLargeTrans_val,R2_veryLargeTrans_val)


pdf(file = "ReadPositions_TopHits_21plus_SmallTranscripts.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title_small = paste("Read Count (Total = ", length(R12_small_positions), ")", sep="") 
hist(R12_small_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_small, main="Distribution of read positions (Transcripts < 1000 bp)")
dev.off()

pdf(file = "ReadPositions_TopHits_21plus_MediumTranscripts.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title_medium = paste("Read Count (Total = ", length(R12_medium_positions), ")", sep="") 
hist(R12_medium_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_medium, main="Distribution of read positions (Transcripts 1000-5000 bp)")
dev.off()

pdf(file = "ReadPositions_TopHits_21plus_LargeTranscripts.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title_large = paste("Read Count (Total = ", length(R12_large_positions), ")", sep="") 
hist(R12_large_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_large, main="Distribution of read positions (Transcripts 5000-10000 bp)")
dev.off()

pdf(file = "ReadPositions_TopHits_21plus_veryLargeTranscripts.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title_veryLarge = paste("Read Count (Total = ", length(R12_veryLarge_positions), ")", sep="") 
hist(R12_veryLarge_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_veryLarge, main="Distribution of read positions (Transcripts > 10000 bp)")
dev.off()

#Make a combo figure with all 4 together
pdf(file = "ReadPositions_TopHits_21plus_Combo.pdf", width = 12, height = 12, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")

par(mfrow=c(2,2))
hist(R12_small_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_small, main="Distribution of read positions (Transcripts < 1000 bp)")
hist(R12_medium_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_medium, main="Distribution of read positions (Transcripts 1000-5000 bp)")
hist(R12_large_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_large, main="Distribution of read positions (Transcripts 5000-10000 bp)")
hist(R12_veryLarge_positions, breaks=100, col="blue", xlab="% Read Position in Transcript", ylab=y_title_veryLarge, main="Distribution of read positions (Transcripts > 10000 bp)")
dev.off()


#6.) The distribution of distances between R1 and R2 (again only for top_hits of length > length_cutoff)
R12_quality_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff & read_data[,"DistanceBetweenReads_Transcript"]!="Gene_Mismatch")

R12_val = as.numeric(read_data[R12_quality_I,"DistanceBetweenReads_Transcript"])

#Median distance between all reads where both reads of a pair map unambiguously to a single gene
# = 177
median(R12_val)
R12_median = median(R12_val)

#Percent that are with +/- 100bp of the median
(length(which(R12_val > (R12_median-100) & R12_val < (R12_median+100))) / length(R12_val))*100


R12_distance_val_small = R12_val[which(R12_val<=500)]
R12_distance_val_large = R12_val[which(R12_val>500)]

pdf(file = "DistanceBetweenReads_SMALL.pdf", width = 8, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title = paste("Read Count (Total = ", length(R12_distance_val_small), ")", sep="") 
hist(as.numeric(R12_distance_val_small), col="red", xlab="Distance between paired reads in Transcript (bp)", ylab=y_title, main="Distribution of distances between paired reads SMALL (<=500 bp apart)", breaks=100)
dev.off()

pdf(file = "DistanceBetweenReads_LARGE.pdf", width = 8, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title = paste("Read Count (Total = ", length(R12_distance_val_large), ")", sep="") 
hist(as.numeric(R12_distance_val_large), col="red", xlab="Distance between paired reads in Transcript (bp)", ylab=y_title, main="Distribution of distances between paired reads LARGE (>500 bp apart)", breaks=100)
dev.off()

pdf(file = "DistanceBetweenReads_ALL.pdf", width = 8, height = 6, onefile = TRUE, family = "Helvetica", title = "R Graphics Output", fonts = NULL, version = "1.1", paper = "special")
y_title = paste("Read Count (Total = ", length(R12_val), ")", sep="") 
hist(as.numeric(R12_val), col="red", xlab="Distance between paired reads in Transcript (bp)", ylab=y_title, main="Distribution of distances between paired reads", breaks=100)
dev.off()



#7.) The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to the SAME GENE - 1520330 (98.3%)
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to DIFFERENT GENES - 26992 (1.7%)
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to the SAME CHROMOSOME - 1522996 (98.4%)
#    The # of Read Pairs, where both reads map unamgiguously (and both are >= length_cutoff) to DIFFERENT CHROMOSOMES - 24326 (1.6%)


R12_both_quality_I = which(read_data[,"R1_HitType"]=="Top_Hit" & read_data[,"R1_AlignmentLength"]>=length_cutoff & read_data[,"R2_HitType"]=="Top_Hit" & read_data[,"R2_AlignmentLength"]>=length_cutoff)
R12_both_quality = read_data[R12_both_quality_I,]

total = length(R12_both_quality_I)
total_same_gene = total - length(which(R12_both_quality[,"DistanceBetweenReads_Transcript"]=="Gene_Mismatch"))
total_diff_genes = length(which(R12_both_quality[,"DistanceBetweenReads_Transcript"]=="Gene_Mismatch"))
total_same_chromosome = total - length(which(R12_both_quality[,"DistanceBetweenReads_Genomic"]=="Chromosome_Mismatch"))
total_diff_chromosomes = length(which(R12_both_quality[,"DistanceBetweenReads_Genomic"]=="Chromosome_Mismatch"))

#Total
total
(total/total)*100

#Total same gene
total_same_gene
(total_same_gene/total)*100

#Total different genes
total_diff_genes
(total_diff_genes/total)*100

#Total same chromosome
total_same_chromosome
(total_same_chromosome/total)*100

#Total different chromosomes
total_diff_chromosomes
(total_diff_chromosomes/total)*100
