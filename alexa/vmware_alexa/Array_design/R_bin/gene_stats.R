#Written by Malachi Griffith
#Generate simple statistics from files containing gene,transcripts and exon information

#This simple script assumes that it is being passed to 'R' FROM the directory containing the following files: 
transcripts_per_gene = read.table("transcript_counts.txt", header=F, quote="", sep="\t", comment.char="")
transcript_lengths = read.table("transcript_lengths.txt", header=F, quote="", sep="\t", comment.char="")
exons_per_transcript = read.table("exon_counts.txt", header=F, quote="", sep="\t", comment.char="")
exon_lengths = read.table("exon_lengths.txt", header=F, quote="", sep="\t", comment.char="")

#Make sure file has some data in it. Some simple species may have 0 known exon skipping events
test = file.info("exon_skip_counts.txt")
if (test["size"] > 0){
  exon_skip_counts = read.table("exon_skip_counts.txt", header=F, quote="", sep="\t", comment.char="")
}


#If you wish to implicitly state their location you can do something like the following
#datadir = "/home/user/alexa/ALEXA_version/stats/genes"
#setwd(datadir)

summary(transcripts_per_gene)
summary(transcript_lengths)
summary(exons_per_transcript)
summary(exon_lengths[,1])

#Create histograms to show the distributions of basic gene statistics 
jpeg(file = "f01_KnownNonpseudo_GeneStats_hist.jpg", quality=100, res=600)
par(mfrow=c(2,2), adj=0.5)
hist(transcripts_per_gene[,1], main="Fig. 1: Transcripts Per Gene", ylab="Gene Count", xlab="Transcripts Per Gene", col="red")
hist(transcript_lengths[,1], main="Fig. 2: Transcript Lengths", ylab="Transcript Count", xlab="Length (bp)", col="blue")
hist(exons_per_transcript[,1], main="Fig. 3: Exons Per Transcript", ylab="Transcript Count", xlab="Exons Per Transcript", col="blue")
hist(exon_lengths[,1], main="Fig. 4: Exon Lengths", ylab="Exon Count", xlab="Length (bp)", col="dark green")
dev.off()

#Transcripts Per Gene - Exclude genes with only a single transcript
genes_multiple_transcripts = transcripts_per_gene[which(transcripts_per_gene[,1]!=1),1]
if (length(genes_multiple_transcripts) > 5){
  jpeg(file = "f05_KnownNonpseudo_TransPerGene_hist.jpg", quality=100, res=600)
  hist(genes_multiple_transcripts, main="Fig. 5: Transcripts Per Gene - excluding singleton genes",
        ylab="Gene Count", xlab="Transcripts Per Gene", col="red")
  dev.off()
}

#Transcript lengths - Divide into large and small transcripts 
transcripts_lt_10000 = transcript_lengths[which(transcript_lengths[,1] <= 10000),1]
transcripts_gt_10000 = transcript_lengths[which(transcript_lengths[,1] > 10000),1]
jpeg(file = "f06_KnownNonpseudo_Transcripts_hist.jpg", quality=100, res=600)
par(mfrow=c(2,1), adj=0.5)
hist(transcripts_lt_10000, main="Fig. 6: Transcript Lengths - Less than 10000 bp", ylab="Transcript Count", xlab="Length (bp)", col="blue")
hist(transcripts_gt_10000, main="Fig. 7: Transcript Lengths - Greater than 10000 bp", ylab="Transcript Count", xlab="Length (bp)", col="blue")
dev.off()

#Exons Per transcript - Divide into those with a reasonable number of exons and those with a large number
exons_lt_50 = exons_per_transcript[which(exons_per_transcript[,1] <= 50),1]
exons_gt_50 = exons_per_transcript[which(exons_per_transcript[,1] > 50),1]
jpeg(file = "f08_KnownNonpseudo_Exons_hist.jpg", quality=100, res=600)
par(mfrow=c(2,1), adj=0.5)
hist(exons_lt_50, main="Fig. 8: Exons Per Transcript", ylab="Transcript Count",
     xlab="Exons Per Transcript - for those with less than 50 exons", col="blue")
if (length(exons_gt_50) > 5){
  hist(exons_gt_50, main="Fig. 9: Exons Per Transcript", ylab="Transcript Count",
       xlab="Exons Per Transcript - for those with more than 50 exons", col="blue")
}
dev.off()

#Exon lengths - Divide into large and small exons
exons_small = exon_lengths[which(exon_lengths[,1] <= 500),1]
exons_large = exon_lengths[which(exon_lengths[,1] > 500),1]
jpeg(file = "f10_KnownNonpseudo_ExonSizes_hist.jpg", quality=100, res=600)
par(mfrow=c(2,1), adj=0.5)
hist(exons_small, main="Fig. 10: Exon Lengths - Exons smaller than 500 bp", ylab="Exon Count", xlab="Length (bp)", col="dark green")
hist(exons_large, main="Fig. 11: Exon Lengths - Exons larger than 500 bp", ylab="Exon Count", xlab="Length (bp)", col="dark green")
dev.off()

#Create boxplots to show the distributions of basic gene statistics 
jpeg(file = "f12_KnownNonpseudo_GeneStats_boxplot.jpg", quality=100, res=600)
par(mfrow=c(2,2), adj=0.5)
boxplot(transcripts_per_gene[,1], main="Fig. 12: Transcripts Per Gene", ylab="Transcript Count", col="blue")
boxplot(transcript_lengths[,1], main="Fig. 13: Transcript Lengths", ylab="Length (bp)", col="blue")
boxplot(exons_per_transcript[,1], main="Fig. 14: Exons Per Transcript", ylab="Exon Count", col="dark green")
boxplot(exon_lengths[,1], main="Fig. 15: Exon Lengths", ylab="Length (bp)", col="dark green")
dev.off()

#Exon skipping - Summarize the number exons that get skipped in every observable exon skipping event in Ensembl
#An analysis of all multi-transcript genes in Ensembl 
if (test["size"] > 0){
  table(exon_skip_counts[,1])
  quantile(exon_skip_counts[,1], probs= seq(0, 1, 0.05))
  jpeg(file = "f13_KnownNonpseudo_ExonsSkipped_hist.jpg", quality=100, res=600)
  hist(exon_skip_counts[,1], main="Fig. 13: Exons skipped per skip event", ylab="Event Count", xlab="Exons Skipped", col="dark green")
  dev.off()
}


