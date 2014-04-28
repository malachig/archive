#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create plots of genes identified and exonic bases covered as a function of increasing read depth
#The input data for this script was generated from read record files using a perl script

#MIP101
#10X coverage files
gene_base_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Gene-ExonicBases_Coverage_10xCoverage_v53_RAND2.txt"
gene_base_multi_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Gene-ExonicBases_Coverage_MULTI.txt"
exon_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Exon_Coverage_10xCoverage_v53.txt"
junction_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Junction_Coverage_10xCoverage_v53.txt"
intron_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Intron_Coverage_10xCoverage_v53.txt"
intergenic_file = "/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Intergenic_Coverage_10xCoverage_v53_RAND.txt"
results_dir = "/projects/malachig/solexa/figures_and_stats/HS04391/coverageVersusDepth_v53"
x_label = "% Library Depth (of 334 million reads)"

#MIP5FU
gene_base_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Gene-ExonicBases_Coverage_10xCoverage_v53_RAND2.txt"
gene_base_multi_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Gene-ExonicBases_Coverage_MULTI.txt"
exon_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Exon_Coverage_10xCoverage_v53_RAND.txt"
junction_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Junction_Coverage_10xCoverage_v53_RAND.txt"
intron_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Intron_Coverage_10xCoverage_v53_RAND.txt"
intergenic_file = "/projects/malachig/solexa/read_records/HS04401/Summary/coverageVsDepth/HS04401_Lanes1-16_Intergenic_Coverage_10xCoverage_v53_RAND.txt"
results_dir = "/projects/malachig/solexa/figures_and_stats/HS04401/coverageVersusDepth_v53"
x_label = "% Library Depth (of 191 million reads)"


setwd(results_dir)
dir()
bg="white"

#Theoretical maximums being targeted:  Genes, Exons, Junctions, Exonic Bases, Intronic Bases, Intergenic Bases
#These need to be updated if the EnsEMBL version used for mapping changes!!
max_genes_t = 36953
max_exons_t = 242709
max_known_junctions_t = 218463
max_novel_junctions_t = 1992798
max_exonic_bases_t = 69793264
max_intronic_bases_t = 594103021
max_intergenic_bases_t = 783298713


#Import Gene and Exonic Base Count data - then convert to % values
data_gene_exonic_base = read.table(gene_base_file, header = FALSE, sep = "", comment.char = "#")
names(data_gene_exonic_base) = c("LibraryDepth", "ExpressedGenes", "ExonicBaseCoverage")
data_gene_exonic_base[1,]
data_length = length(data_gene_exonic_base[,1])
max_library_depth = data_gene_exonic_base[data_length, "LibraryDepth"]
data_gene_exonic_base[,"LibraryDepth_p"] = (data_gene_exonic_base[,"LibraryDepth"]/max_library_depth)*100
max_expressed_genes = data_gene_exonic_base[data_length, "ExpressedGenes"]
data_gene_exonic_base[,"ExpressedGenes_p"] = (data_gene_exonic_base[,"ExpressedGenes"]/max_expressed_genes)*100
data_gene_exonic_base[,"ExpressedGenes_pt"] = (data_gene_exonic_base[,"ExpressedGenes"]/max_genes_t)*100
max_exon_coverage = data_gene_exonic_base[data_length, "ExonicBaseCoverage"]
data_gene_exonic_base[,"ExonicBaseCoverage_p"] = (data_gene_exonic_base[,"ExonicBaseCoverage"]/max_exon_coverage)*100
data_gene_exonic_base[,"ExonicBaseCoverage_pt"] = (data_gene_exonic_base[,"ExonicBaseCoverage"]/max_exonic_bases_t)*100

i = 2:length(data_gene_exonic_base[,1])
data_gene_exonic_base[,"discovery_rate_exonic_p"] = c(data_gene_exonic_base[1,"ExonicBaseCoverage_p"], data_gene_exonic_base[i,"ExonicBaseCoverage_p"] - data_gene_exonic_base[i-1,"ExonicBaseCoverage_p"])
data_gene_exonic_base[,"discovery_rate_gene_p"] = c(data_gene_exonic_base[1,"ExpressedGenes_p"], data_gene_exonic_base[i,"ExpressedGenes_p"] - data_gene_exonic_base[i-1,"ExpressedGenes_p"])
data_gene_exonic_base[,"discovery_rate_exonic_pt"] = c(data_gene_exonic_base[1,"ExonicBaseCoverage_pt"], data_gene_exonic_base[i,"ExonicBaseCoverage_pt"] - data_gene_exonic_base[i-1,"ExonicBaseCoverage_pt"])
data_gene_exonic_base[,"discovery_rate_gene_pt"] = c(data_gene_exonic_base[1,"ExpressedGenes_pt"], data_gene_exonic_base[i,"ExpressedGenes_pt"] - data_gene_exonic_base[i-1,"ExpressedGenes_pt"])

data_gene_exonic_base[,"discovery_rate_exonic_r"] = c((data_gene_exonic_base[1,"ExonicBaseCoverage"]/data_gene_exonic_base[1,"ExonicBaseCoverage"])*100, ((data_gene_exonic_base[i,"ExonicBaseCoverage"] - data_gene_exonic_base[i-1,"ExonicBaseCoverage"])/data_gene_exonic_base[i-1,"ExonicBaseCoverage"])*100)


#Import Exonic base count data for multiple cutoffs - then convert to % values
data_gene_exonic_base_multi = read.table(gene_base_multi_file, header = FALSE, sep = "", comment.char = "#")
names(data_gene_exonic_base_multi) = c("LibraryDepth","Genes_1X","Genes_5X","Genes_10X","Genes_50X","Genes_100X","Genes_500X","ExonicBases_1X","ExonicBases_5X","ExonicBases_10X","ExonicBases_50X","ExonicBases_100X","ExonicBases_500X")
data_length = length(data_gene_exonic_base_multi[,1])
max_library_depth = data_gene_exonic_base_multi[data_length, "LibraryDepth"]
data_gene_exonic_base_multi[,"LibraryDepth_p"] = (data_gene_exonic_base_multi[,"LibraryDepth"]/max_library_depth)*100

data_gene_exonic_base_multi[,"ExonicBaseCoverage_1x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_1X"]/max_exonic_bases_t)*100
data_gene_exonic_base_multi[,"ExonicBaseCoverage_5x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_5X"]/max_exonic_bases_t)*100
data_gene_exonic_base_multi[,"ExonicBaseCoverage_10x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_10X"]/max_exonic_bases_t)*100
data_gene_exonic_base_multi[,"ExonicBaseCoverage_50x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_50X"]/max_exonic_bases_t)*100
data_gene_exonic_base_multi[,"ExonicBaseCoverage_100x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_100X"]/max_exonic_bases_t)*100
data_gene_exonic_base_multi[,"ExonicBaseCoverage_500x_pt"] = (data_gene_exonic_base_multi[,"ExonicBases_500X"]/max_exonic_bases_t)*100

#Get the discovery rate at each iteration for each cutoff
i = 2:length(data_gene_exonic_base_multi[,1])
data_gene_exonic_base_multi[,"discovery_rate_1x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_1x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_1x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_1x_pt"])
data_gene_exonic_base_multi[,"discovery_rate_5x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_5x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_5x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_5x_pt"])
data_gene_exonic_base_multi[,"discovery_rate_10x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_10x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_10x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_10x_pt"])
data_gene_exonic_base_multi[,"discovery_rate_50x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_50x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_50x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_50x_pt"])
data_gene_exonic_base_multi[,"discovery_rate_100x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_100x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_100x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_100x_pt"])
data_gene_exonic_base_multi[,"discovery_rate_500x"] = c(data_gene_exonic_base_multi[1,"ExonicBaseCoverage_500x_pt"], data_gene_exonic_base_multi[i,"ExonicBaseCoverage_500x_pt"] - data_gene_exonic_base_multi[i-1,"ExonicBaseCoverage_500x_pt"])


#Import Exon data - then convert to % values
data_exon = read.table(exon_file, header = FALSE, sep = "", comment.char = "#")
names(data_exon) = c("LibraryDepth","ExpressedExons")
data_exon[1,]
data_length = length(data_exon[,1])
max_library_depth = data_exon[data_length, "LibraryDepth"]
data_exon[,"LibraryDepth_p"] = (data_exon[,"LibraryDepth"]/max_library_depth)*100

max_expressed_exons = data_exon[data_length, "ExpressedExons"]
data_exon[,"ExpressedExons_p"] = (data_exon[,"ExpressedExons"]/max_expressed_exons)*100
data_exon[,"ExpressedExons_pt"] = (data_exon[,"ExpressedExons"]/max_exons_t)*100
i = 2:length(data_exon[,1])
data_exon[,"discovery_rate_exon_p"] = c(data_exon[1,"ExpressedExons_p"], data_exon[i,"ExpressedExons_p"] - data_exon[i-1,"ExpressedExons_p"])
data_exon[,"discovery_rate_exon_pt"] = c(data_exon[1,"ExpressedExons_pt"], data_exon[i,"ExpressedExons_pt"] - data_exon[i-1,"ExpressedExons_pt"])

#Import Junction data - then convert to % values
data_junction = read.table(junction_file, header = FALSE, sep = "", comment.char = "#")
names(data_junction) = c("LibraryDepth","ExpressedKnownJunctions","ExpressedNovelJunctions")
data_junction[1,]
data_length = length(data_junction[,1])
max_library_depth = data_junction[data_length, "LibraryDepth"]
data_junction[,"LibraryDepth_p"] = (data_junction[,"LibraryDepth"]/max_library_depth)*100

max_expressed_known_junctions = data_junction[data_length, "ExpressedKnownJunctions"]
max_expressed_novel_junctions = data_junction[data_length, "ExpressedNovelJunctions"]

data_junction[,"ExpressedKnownJunctions_p"] = (data_junction[,"ExpressedKnownJunctions"]/max_expressed_known_junctions)*100
data_junction[,"ExpressedKnownJunctions_pt"] = (data_junction[,"ExpressedKnownJunctions"]/max_known_junctions_t)*100
data_junction[,"ExpressedNovelJunctions_p"] = (data_junction[,"ExpressedNovelJunctions"]/max_expressed_novel_junctions)*100
data_junction[,"ExpressedNovelJunctions_pt"] = (data_junction[,"ExpressedNovelJunctions"]/max_novel_junctions_t)*100

i = 2:length(data_junction[,1])
data_junction[,"discovery_rate_known_junction_p"] = c(data_junction[1,"ExpressedKnownJunctions_p"], data_junction[i,"ExpressedKnownJunctions_p"] - data_junction[i-1,"ExpressedKnownJunctions_p"])
data_junction[,"discovery_rate_novel_junction_p"] = c(data_junction[1,"ExpressedNovelJunctions_p"], data_junction[i,"ExpressedNovelJunctions_p"] - data_junction[i-1,"ExpressedNovelJunctions_p"])
data_junction[,"discovery_rate_known_junction_pt"] = c(data_junction[1,"ExpressedKnownJunctions_pt"], data_junction[i,"ExpressedKnownJunctions_pt"] - data_junction[i-1,"ExpressedKnownJunctions_pt"])
data_junction[,"discovery_rate_novel_junction_pt"] = c(data_junction[1,"ExpressedNovelJunctions_pt"], data_junction[i,"ExpressedNovelJunctions_pt"] - data_junction[i-1,"ExpressedNovelJunctions_pt"])

data_junction[,"discovery_rate_known_junction_r"] = c((data_junction[1,"ExpressedKnownJunctions"]/data_junction[1,"ExpressedKnownJunctions"])*100, ((data_junction[i,"ExpressedKnownJunctions"] - data_junction[i-1,"ExpressedKnownJunctions"])/data_junction[i-1,"ExpressedKnownJunctions"])*100)
data_junction[,"discovery_rate_novel_junction_r"] = c((data_junction[1,"ExpressedNovelJunctions"]/data_junction[1,"ExpressedNovelJunctions"])*100, ((data_junction[i,"ExpressedNovelJunctions"] - data_junction[i-1,"ExpressedNovelJunctions"])/data_junction[i-1,"ExpressedNovelJunctions"])*100)


#Import Intronic Base Count data
data_intronic_base = read.table(intron_file, header = FALSE, sep = "", comment.char = "#")
names(data_intronic_base) = c("LibraryDepth", "IntronicBaseCoverage")
data_intronic_base[1,]
data_length = length(data_intronic_base[,1])
max_library_depth = data_intronic_base[data_length, "LibraryDepth"]
data_intronic_base[,"LibraryDepth_p"] = (data_intronic_base[,"LibraryDepth"]/max_library_depth)*100
max_intron_coverage = data_intronic_base[data_length, "IntronicBaseCoverage"]
data_intronic_base[,"IntronicBaseCoverage_p"] = (data_intronic_base[,"IntronicBaseCoverage"]/max_intron_coverage)*100
data_intronic_base[,"IntronicBaseCoverage_pt"] = (data_intronic_base[,"IntronicBaseCoverage"]/max_intronic_bases_t)*100
i = 2:length(data_intronic_base[,1])
data_intronic_base[,"discovery_rate_intronic_p"] = c(data_intronic_base[1,"IntronicBaseCoverage_p"], data_intronic_base[i,"IntronicBaseCoverage_p"] - data_intronic_base[i-1,"IntronicBaseCoverage_p"])
data_intronic_base[,"discovery_rate_intronic_pt"] = c(data_intronic_base[1,"IntronicBaseCoverage_pt"], data_intronic_base[i,"IntronicBaseCoverage_pt"] - data_intronic_base[i-1,"IntronicBaseCoverage_pt"])
data_intronic_base[,"discovery_rate_intronic_r"] = c((data_intronic_base[1,"IntronicBaseCoverage"]/data_intronic_base[1,"IntronicBaseCoverage"])*100, ((data_intronic_base[i,"IntronicBaseCoverage"] - data_intronic_base[i-1,"IntronicBaseCoverage"])/data_intronic_base[i-1,"IntronicBaseCoverage"])*100)


#Import Intergenic Base Count data
data_intergenic_base = read.table(intergenic_file, header = FALSE, sep = "", comment.char = "#")
names(data_intergenic_base) = c("LibraryDepth", "IntergenicBaseCoverage")
data_intergenic_base[1,]
data_length = length(data_intergenic_base[,1])
max_library_depth = data_intergenic_base[data_length, "LibraryDepth"]
data_intergenic_base[,"LibraryDepth_p"] = (data_intergenic_base[,"LibraryDepth"]/max_library_depth)*100
max_intergenic_coverage = data_intergenic_base[data_length, "IntergenicBaseCoverage"]
data_intergenic_base[,"IntergenicBaseCoverage_p"] = (data_intergenic_base[,"IntergenicBaseCoverage"]/max_intergenic_coverage)*100
data_intergenic_base[,"IntergenicBaseCoverage_pt"] = (data_intergenic_base[,"IntergenicBaseCoverage"]/max_intergenic_bases_t)*100
i = 2:length(data_intergenic_base[,1])
data_intergenic_base[,"discovery_rate_intergenic_p"] = c(data_intergenic_base[1,"IntergenicBaseCoverage_p"], data_intergenic_base[i,"IntergenicBaseCoverage_p"] - data_intergenic_base[i-1,"IntergenicBaseCoverage_p"])
data_intergenic_base[,"discovery_rate_intergenic_pt"] = c(data_intergenic_base[1,"IntergenicBaseCoverage_pt"], data_intergenic_base[i,"IntergenicBaseCoverage_pt"] - data_intergenic_base[i-1,"IntergenicBaseCoverage_pt"])
data_intergenic_base[,"discovery_rate_intergenic_r"] = c((data_intergenic_base[1,"IntergenicBaseCoverage"]/data_intergenic_base[1,"IntergenicBaseCoverage"])*100, ((data_intergenic_base[i,"IntergenicBaseCoverage"] - data_intergenic_base[i-1,"IntergenicBaseCoverage"])/data_intergenic_base[i-1,"IntergenicBaseCoverage"])*100)



#First plot as the percent of all possible theoretical events
cols = c("blue", "dark green", "red", "grey10", "grey50")
tiff(file="CoverageVersusDepth_Plot_PercentOfMax.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Percent of events detected with increasing library depth"
y_label = "% of possible events detected at >= 10x"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"ExpressedGenes_pt"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,65), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

#lines(x=data_exon[,"LibraryDepth_p"], y=data_exon[,"ExpressedExons_pt"], lty=1, col=cols[2], lwd=2)
lines(x=data_junction[,"LibraryDepth_p"], y=data_junction[,"ExpressedKnownJunctions_pt"], lty=1, col=cols[2], lwd=2)

lines(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"ExonicBaseCoverage_pt"], lty=2, col=cols[3], lwd=2)
lines(x=data_intronic_base[,"LibraryDepth_p"], y=data_intronic_base[,"IntronicBaseCoverage_pt"], lty=2, col=cols[4], lwd=2)
lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"IntergenicBaseCoverage_pt"], lty=2, col=cols[5], lwd=2)

genes_text = paste ("Genes (n = ",  prettyNum(max_expressed_genes, big.mark=","), ")", sep="")
#exons_text = paste ("Exons (n = ",  prettyNum(max_expressed_exons, big.mark=","), ")", sep="")
junctions_text = paste ("Junctions (n = ",  prettyNum(max_expressed_known_junctions, big.mark=","), ")", sep="")
exonic_bases_text = paste ("Exonic Bases (n = ",  prettyNum(max_exon_coverage, big.mark=","), ")", sep="")
intronic_bases_text = paste ("Intronic Bases (n = ",  prettyNum(max_intron_coverage, big.mark=","), ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  prettyNum(max_intergenic_coverage, big.mark=","), ")", sep="")

legend_text = c(genes_text, junctions_text, exonic_bases_text, intronic_bases_text, intergenic_bases_text)
legend_position = "topleft"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1,2,2,2), lwd=2, cex=1.4)
dev.off();

#Now plot as the percent of all events actually observed
cols = c("blue", "dark green", "red", "grey10", "grey50")
tiff(file="CoverageVersusDepth_Plot_PercentOfObserved.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Events identified with increasing library depth"
y_label = "% of events detected at >= 10x"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"ExpressedGenes_p"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,100), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

#lines(x=data_exon[,"LibraryDepth_p"], y=data_exon[,"ExpressedExons_pt"], lty=1, col=cols[2], lwd=2)
lines(x=data_junction[,"LibraryDepth_p"], y=data_junction[,"ExpressedKnownJunctions_p"], lty=1, col=cols[2], lwd=2)

lines(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"ExonicBaseCoverage_p"], lty=2, col=cols[3], lwd=2)
lines(x=data_intronic_base[,"LibraryDepth_p"], y=data_intronic_base[,"IntronicBaseCoverage_p"], lty=2, col=cols[4], lwd=2)
lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"IntergenicBaseCoverage_p"], lty=2, col=cols[5], lwd=2)

genes_text = paste ("Genes (n = ",  max_expressed_genes, ")", sep="")
#exons_text = paste ("Exons (n = ",  max_expressed_exons, ")", sep="")
junctions_text = paste ("Junctions (n = ",  prettyNum(max_expressed_known_junctions, big.mark=","), ")", sep="")
exonic_bases_text = paste ("Exonic Bases (n = ",  prettyNum(max_exon_coverage, big.mark=","), ")", sep="")
intronic_bases_text = paste ("Intronic Bases (n = ",  prettyNum(max_intron_coverage, big.mark=","), ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  prettyNum(max_intergenic_coverage, big.mark=","), ")", sep="")

legend_text = c(genes_text, junctions_text, exonic_bases_text, intronic_bases_text, intergenic_bases_text)
legend_position = "bottomright"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1,2,2,2), lwd=2)
dev.off();


#Now plot the discovery rates (change in percent of ALL POSSIBLE events of each type observed from one interation to the next)
cols = c("blue", "dark green", "orange", "red")
tiff(file="CoverageVersusDepth_Plot_Discovery_pt.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate with increasing library depth"
y_label = "% discovery rate of events detected at >= 10x (% of theoretical events)"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"discovery_rate_exonic_pt"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,0.1), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(x=data_junction[,"LibraryDepth_p"], y=data_junction[,"discovery_rate_known_junction_pt"], lty=1, col=cols[2], lwd=2)
lines(x=data_intronic_base[,"LibraryDepth_p"], y=data_intronic_base[,"discovery_rate_intronic_pt"], lty=2, col=cols[3], lwd=2)
lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"discovery_rate_intergenic_pt"], lty=2, col=cols[4], lwd=2)
abline(h=0, col="grey50", lwd=1, lty=3)

exonic_bases_text = paste ("Exonic Bases (n = ",  prettyNum(max_exon_coverage, big.mark=","), ")", sep="")
junctions_text = paste ("Junctions (n = ",  prettyNum(max_expressed_known_junctions, big.mark=","), ")", sep="")
intronic_bases_text = paste ("Intronic Bases (n = ",  prettyNum(max_intron_coverage, big.mark=","), ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  prettyNum(max_intergenic_coverage, big.mark=","), ")", sep="")

legend_text = c(exonic_bases_text, junctions_text, intronic_bases_text, intergenic_bases_text)
legend_position = "topright"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1,2,2), lwd=2)
dev.off();

#At what point would the two discover rate lines (exonic and intergenic or intronic) intersect?  Fit data to tail of distribution only - LINEAR MODEL FIT and POWER MODEL FIT
tail = 0.9
exonic_lower = floor((length(data_gene_exonic_base[,1]))*tail)
exonic_upper = length(data_gene_exonic_base[,1])
intergenic_lower = floor((length(data_intergenic_base[,1]))*tail)
intergenic_upper = length(data_intergenic_base[,1])
intronic_lower = floor((length(data_intronic_base[,1]))*tail)
intronic_upper = length(data_intronic_base[,1])

#LINEAR MODEL
exonic_fit = lm(data_gene_exonic_base[exonic_lower:exonic_upper,"discovery_rate_exonic_pt"]~data_gene_exonic_base[exonic_lower:exonic_upper,"LibraryDepth_p"])
intergenic_fit = lm(data_intergenic_base[intergenic_lower:intergenic_upper,"discovery_rate_intergenic_pt"]~data_intergenic_base[intergenic_lower:intergenic_upper,"LibraryDepth_p"])
intronic_fit = lm(data_intronic_base[intronic_lower:intronic_upper,"discovery_rate_intronic_pt"]~data_intronic_base[intronic_lower:intronic_upper,"LibraryDepth_p"])

#Find intersection between two straight lines
#y1 = mx2+b and y2 = mx2+b. set y1=y2.  solve for x. Then plug into one of the two equations.
#exonic vs. intergenic
slope_diff = exonic_fit$coeff[2]-intergenic_fit$coeff[2]
int_diff = intergenic_fit$coeff[1]-exonic_fit$coeff[1]
library_depth_needed_intergenic = int_diff/slope_diff
corresponding_discovery_rate_intergenic = (exonic_fit$coeff[2]*x)+exonic_fit$coeff[1]
library_depth_needed_intergenic

#exonic vs. intronic
slope_diff = exonic_fit$coeff[2]-intronic_fit$coeff[2]
int_diff = intronic_fit$coeff[1]-exonic_fit$coeff[1]
library_depth_needed_intronic = int_diff/slope_diff
corresponding_discovery_rate_intronic = (exonic_fit$coeff[2]*x)+exonic_fit$coeff[1]
library_depth_needed_intronic

#Redraw the curve above to depict the prejected point of true saturation
cols = c("blue", "red")
tiff(file="CoverageVersusDepth_Plot_Discovery_pt_projection_linear.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate with increasing library depth"
y_label = "% discovery rate of events detected at >= 10x (% of theoretical events)"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"discovery_rate_exonic_pt"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,150), ylim=c(0,0.025), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"discovery_rate_intergenic_pt"], lty=1, col=cols[2], lwd=2)
abline(h=0, col="grey50", lwd=1, lty=3)
abline(exonic_fit, col=cols[1], lty=2)
abline(intergenic_fit, col=cols[2], lty=2)
exonic_bases_text = paste ("Exonic Bases (n = ",  prettyNum(max_exon_coverage, big.mark=","), ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  prettyNum(max_intergenic_coverage, big.mark=","), ")", sep="")

legend_text = c(exonic_bases_text, intergenic_bases_text)
legend_position = "topright"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1), lwd=2)
dev.off();


#POWER MODEL FIT
#Y = a *X^b

#Determine the best amount of tail data to use for the following model fit...
tails=seq(0.05,0.5,0.01)
for (i in tails){
  tail=i
  exonic_lower = floor((length(data_gene_exonic_base[,1]))*tail)
  exonic_upper = length(data_gene_exonic_base[,1])
  x = data_gene_exonic_base[exonic_lower:exonic_upper,"LibraryDepth_p"]
  y = data_gene_exonic_base[exonic_lower:exonic_upper,"discovery_rate_exonic_pt"]
  exonic_fit = nls(y~a*x^b, start = list(a = 0.12345, b = 0.54321))
  r = round(((cor(data_gene_exonic_base[exonic_lower:exonic_upper,"discovery_rate_exonic_pt"], predict(exonic_fit)))^2), digits=3)
  print(c(tail, r))
}

#MIP101 setting
tail = 0.13
x_limit=275
x_limit_test = 257

#MIP5FU setting
tail = 0.23
x_limit=475 
x_limit_test = 440

exonic_lower = floor((length(data_gene_exonic_base[,1]))*tail)
exonic_upper = length(data_gene_exonic_base[,1])
intergenic_lower = floor((length(data_intergenic_base[,1]))*tail)
intergenic_upper = length(data_intergenic_base[,1])

x = data_gene_exonic_base[exonic_lower:exonic_upper,"LibraryDepth_p"]
y = data_gene_exonic_base[exonic_lower:exonic_upper,"discovery_rate_exonic_pt"]
exonic_fit = nls(y~a*x^b, start = list(a = 0.12345, b = 0.54321), trace = TRUE)

x = data_intergenic_base[intergenic_lower:intergenic_upper,"LibraryDepth_p"]
y = data_intergenic_base[intergenic_lower:intergenic_upper,"discovery_rate_intergenic_pt"]
intergenic_fit = nls(y~a*x^b, start = list(a = 0.002910054, b = -0.137870576), trace = TRUE)

#Get the R^2 values for each fit (use the equation: R^2 = (cor(fitted.y, y))^2
R_squared_exonic = round(((cor(data_gene_exonic_base[exonic_lower:exonic_upper,"discovery_rate_exonic_pt"], predict(exonic_fit)))^2), digits=3)
R_squared_intergenic = round(((cor(data_intergenic_base[intergenic_lower:intergenic_upper,"discovery_rate_intergenic_pt"], predict(intergenic_fit)))^2), digits=3)

tiff(file="CoverageVersusDepth_Plot_Discovery_pt_projection_power.tiff", width=700, height=700, compression="none")
cols = c("light blue", "blue", "orange", "red", "black")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate with increasing library depth"
y_label = "% discovery rate of events detected at >= 10x (% of theoretical events)"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"discovery_rate_exonic_pt"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,x_limit), ylim=c(0,0.025), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"discovery_rate_intergenic_pt"], lty=1, col=cols[3], lwd=2)
abline(h=0, col="grey50", lwd=1, lty=3)
exonic_bases_text = paste ("Exonic Bases", sep="")
exonic_fit_text = paste ("Power curve fit (R = ", R_squared_exonic, ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases", sep="")
intergenic_fit_text = paste ("Power curve fit (R = ", R_squared_intergenic, ")", sep="")

legend_text = c(exonic_bases_text, exonic_fit_text, intergenic_bases_text, intergenic_fit_text, "point of convergence")
legend_position = "topright"
legend(legend_position, legend=legend_text, col=c(cols[1],cols[2],cols[3],cols[4],cols[5]), lty=c(1,2,1,2,0), pch=c(NA,NA,NA,NA,16), lwd=c(3,3,3,3,NA), cex=1.4)

lines(x=seq(0,x_limit,0.1), y=predict(exonic_fit, list(x=seq(0,x_limit,0.1))), lty=2, lwd=3, col=cols[2])
lines(x=seq(0,x_limit,0.1), y=predict(intergenic_fit, list(x=seq(0,x_limit,0.1))), lty=2, lwd=3, col=cols[4])

#NOTE.  not sure how to calculate the intersection between two power curves, but you can observe the convergence in these data by the following
#the two curves and determining where they cross by subtracting the y-value from one from the y-value of the other
predict(exonic_fit, list(x=seq(0,x_limit_test,1)))-predict(intergenic_fit, list(x=seq(0,x_limit_test,1)))

#Once you have determined the correct x_limit_test value you can grab the corresponding y-value
y_result = predict(exonic_fit, list(x=x_limit_test))

points(x=x_limit_test, y=y_result, pch=16, col=cols[5])
dev.off();


#Now plot the discovery rates (change in percent of ALL OBSERVED events of each type observed from one interation to the next)
cols = c("blue", "dark green", "orange", "red")
tiff(file="CoverageVersusDepth_Plot_Discovery_p.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate with increasing library depth"
y_label = "% discovery rate of events detected at >= 10x (% of observed events)"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"discovery_rate_exonic_p"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(x=data_junction[,"LibraryDepth_p"], y=data_junction[,"discovery_rate_known_junction_p"], lty=1, col=cols[2], lwd=2)
lines(x=data_intronic_base[,"LibraryDepth_p"], y=data_intronic_base[,"discovery_rate_intronic_p"], lty=2, col=cols[3], lwd=2)
lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"discovery_rate_intergenic_p"], lty=2, col=cols[4], lwd=2)
abline(h=0, col="grey50", lwd=1, lty=3)

exonic_bases_text = paste ("Exonic Bases (n = ",  max_exon_coverage, ")", sep="")
junctions_text = paste ("Junctions (n = ",  max_expressed_known_junctions, ")", sep="")
intronic_bases_text = paste ("Intronic Bases (n = ",  max_intron_coverage, ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  max_intergenic_coverage, ")", sep="")

legend_text = c(exonic_bases_text, junctions_text, intronic_bases_text, intergenic_bases_text)
legend_position = "topright"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1,2,2), lwd=2)
dev.off();

#Now plot the discovery rates (change in percent of RELATIVE only to the previous interation)
cols = c("blue", "dark green", "orange", "red")
tiff(file="CoverageVersusDepth_Plot_Discovery_r.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate with increasing library depth"
y_label = "% discovery rate of events detected at >= 10x (% of observed events)"
plot(x=data_gene_exonic_base[,"LibraryDepth_p"], y=data_gene_exonic_base[,"discovery_rate_exonic_r"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,1), type="l", lty=1, col=cols[1], lwd=2, 
     col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)

lines(x=data_junction[,"LibraryDepth_p"], y=data_junction[,"discovery_rate_known_junction_r"], lty=1, col=cols[2], lwd=2)
lines(x=data_intronic_base[,"LibraryDepth_p"], y=data_intronic_base[,"discovery_rate_intronic_r"], lty=2, col=cols[3], lwd=2)
lines(x=data_intergenic_base[,"LibraryDepth_p"], y=data_intergenic_base[,"discovery_rate_intergenic_r"], lty=2, col=cols[4], lwd=2)
abline(h=0, col="grey50", lwd=1, lty=3)

exonic_bases_text = paste ("Exonic Bases (n = ",  max_exon_coverage, ")", sep="")
junctions_text = paste ("Junctions (n = ",  max_expressed_known_junctions, ")", sep="")
intronic_bases_text = paste ("Intronic Bases (n = ",  max_intron_coverage, ")", sep="")
intergenic_bases_text = paste ("Intergenic Bases (n = ",  max_intergenic_coverage, ")", sep="")

legend_text = c(exonic_bases_text, junctions_text, intronic_bases_text, intergenic_bases_text)
legend_position = "topright"
legend(legend_position, legend=legend_text, col=cols, lty=c(1,1,2,2), lwd=2)
dev.off();




###############################################################################################################################
#Now generate a plot to compare the identification of exonic bases at increasing cutoff requirements (1x, 5x, 10x, 50x, 100x, 500x)

#Get the slope of the line fit to the following sets of data (using the last 10% of the line) to estimate the rate of discovery
lower = floor(length(data_gene_exonic_base_multi[,1])*0.9)
upper = length(data_gene_exonic_base_multi[,1])
slope_1x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_1x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)
slope_5x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_5x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)
slope_10x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_10x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)
slope_50x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_50x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)
slope_100x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_100x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)
slope_500x = round(as.numeric(lm(data_gene_exonic_base_multi[lower:upper,"ExonicBaseCoverage_500x_pt"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])$coef[2]), digits=3)

legend_text = c(paste("1X cutoff (slope = ", slope_1x, ")", sep=""),
                paste("5X cutoff (slope = ", slope_5x, ")", sep=""),
                paste("10X cutoff (slope = ", slope_10x, ")", sep=""),
                paste("50X cutoff (slope = ", slope_50x, ")", sep=""),
                paste("100X cutoff (slope = ", slope_100x, ")", sep=""),
                paste("500X cutoff (slope = ", slope_500x, ")", sep=""))

cols = rainbow(6)
tiff(file="CoverageVersusDepth_Plot_MultiCutoff.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Exonic bases identified with increasing library depth (at increasing cutoffs)"
y_label = "% of possible events detected at specified cutoff"
plot(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_1x_pt"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,100), type="l", lty=1, col=cols[1], lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_5x_pt"], lty=1, col=cols[2], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_10x_pt"], lty=1, col=cols[3], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_50x_pt"], lty=1, col=cols[4], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_100x_pt"], lty=1, col=cols[5], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"ExonicBaseCoverage_500x_pt"], lty=1, col=cols[6], lwd=2)
legend("topright", legend_text, lty=1, lwd=2, col=cols, cex=1.4)
dev.off()


#Now plot the discovery rates instead of cumulative percent discovered

#Try fitting linear model to tail of the data
#NOTE: A linear model is not ideal here and would underestimate the additional depth needed to achieve the same minimum discovery rate achieved at 1x cutoff...
#Use fit to the last 10% of the data
lower = 1839
upper = length(data_gene_exonic_base_multi[,1])
fit_1x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_1x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])
fit_5x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_5x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])
fit_10x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_10x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])
fit_50x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_50x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])
fit_100x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_100x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])
fit_500x = lm(data_gene_exonic_base_multi[lower:upper,"discovery_rate_500x"]~data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"])

#What does the discovery rate drop to at 100% for the 1x line?  y = mX+b
min_discovery_1x = fit_1x$coeff[2]*100+fit_1x$coeff[1]

#At what %library size would this discovery rate be achieved for each of the other cutoff values?
#x = (y - b)/m
library_size_needed_1x = (min_discovery_1x - fit_1x$coeff[1])/fit_1x$coeff[2] 
library_size_needed_5x = (min_discovery_1x - fit_5x$coeff[1])/fit_5x$coeff[2] 
library_size_needed_10x = (min_discovery_1x - fit_10x$coeff[1])/fit_10x$coeff[2] 
library_size_needed_50x = (min_discovery_1x - fit_50x$coeff[1])/fit_50x$coeff[2] 
library_size_needed_100x = (min_discovery_1x - fit_100x$coeff[1])/fit_100x$coeff[2] 
library_size_needed_500x = (min_discovery_1x - fit_500x$coeff[1])/fit_500x$coeff[2] 

legend_text = c(paste("1X cutoff", sep=""),
                paste("5X cutoff", sep=""),
                paste("10X cutoff", sep=""),
                paste("50X cutoff", sep=""),
                paste("100X cutoff", sep=""),
                paste("500X cutoff", sep=""))

cols = rainbow(6)
tiff(file="CoverageVersusDepth_Plot_MultiCutoff_Discovery.tiff", width=700, height=700, compression="none")
par(bg=bg, font.main = 2, font.lab = 2); 
main_title = "Discovery rate of exonic bases with increasing library depth (at increasing cutoffs)"
y_label = "% discovery rate"
plot(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_1x"], main=main_title, xlab=x_label, ylab=y_label,
     xlim=c(0,100), ylim=c(0,0.05),type="l", lty=1, col=cols[1], lwd=2, col.lab = gray(.1), cex.main = 1.5, cex.lab = 1.4, cex.axis=1.4)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_5x"], lty=1, col=cols[2], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_10x"], lty=1, col=cols[3], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_50x"], lty=1, col=cols[4], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_100x"], lty=1, col=cols[5], lwd=2)
lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=data_gene_exonic_base_multi[,"discovery_rate_500x"], lty=1, col=cols[6], lwd=2)
legend("topright", legend_text, lty=1, lwd=2, col=cols)

#lines(x=data_gene_exonic_base_multi[lower:upper,"LibraryDepth_p"], y=predict(fit_1x), lty=2)
#lines(x=data_gene_exonic_base_multi[,"LibraryDepth_p"], y=predict.lm(fit_1x, newdata=data_gene_exonic_base_multi[,"discovery_rate_1x"]), lty=2)
abline(fit_1x, lty=2, col=cols[1])
abline(fit_5x, lty=2, col=cols[2])
abline(fit_10x, lty=2, col=cols[3])
abline(fit_50x, lty=2, col=cols[4])
abline(fit_100x, lty=2, col=cols[5])
abline(fit_500x, lty=2, col=cols[6])
dev.off()
















###############################################################################################################################################33
#Old code
hs04391_data_1x = read.table("HS04391_Lanes1-23_Coverage_100kblocks_1ReadGenes_1xCoverage.txt", header = FALSE, sep = "", comment.char = "#")
hs04391_data_10x = read.table("HS04391_Lanes1-23_Coverage_100kblocks_10ReadGenes_10xCoverage.txt", header = FALSE, sep = "", comment.char = "#")
colnames(hs04391_data_1x) = c("MappedReads","ExpressedGenes","ExonCoverage")
colnames(hs04391_data_10x) = c("MappedReads","ExpressedGenes","ExonCoverage")

hs04401_data_1x = read.table("HS04401_Lanes1-16_Coverage_100kblocks_1ReadGenes_1xCoverage.txt", header = FALSE, sep = "", comment.char = "#")
hs04401_data_10x = read.table("HS04401_Lanes1-16_Coverage_100kblocks_10ReadGenes_10xCoverage.txt", header = FALSE, sep = "", comment.char = "#")
colnames(hs04401_data_1x) = c("MappedReads","ExpressedGenes","ExonCoverage")
colnames(hs04401_data_10x) = c("MappedReads","ExpressedGenes","ExonCoverage")

#Calculate the number of NEW genes and covered bases after each iteration of data was added
calculateDiff = function(data, old_col, new_col){
  current = 0;
  for(i in 1:length(data[,old_col])){
    diff = data[i, old_col] - current
    data[i,new_col] = diff
    current = data[i,old_col] 
  }
  return(data)
}
hs04391_data_1x = calculateDiff(hs04391_data_1x, "ExpressedGenes", "GenesDiff")
hs04391_data_1x = calculateDiff(hs04391_data_1x, "ExonCoverage", "CoverageDiff")
hs04391_data_10x = calculateDiff(hs04391_data_10x, "ExpressedGenes", "GenesDiff")
hs04391_data_10x = calculateDiff(hs04391_data_10x, "ExonCoverage", "CoverageDiff")

hs04401_data_1x = calculateDiff(hs04401_data_1x, "ExpressedGenes", "GenesDiff")
hs04401_data_1x = calculateDiff(hs04401_data_1x, "ExonCoverage", "CoverageDiff")
hs04401_data_10x = calculateDiff(hs04401_data_10x, "ExpressedGenes", "GenesDiff")
hs04401_data_10x = calculateDiff(hs04401_data_10x, "ExonCoverage", "CoverageDiff")


################################################
#Plot Exonic base coverage versus read depth
total_exonic_bases = 66915856 #Max theoretically possible
mapability_rate = 61503886/65920000 #The overal percentage of mapable bases
total_mapable_exonic_bases = total_exonic_bases*mapability_rate 

#HS04391
x1 = hs04391_data_1x[,"MappedReads"]
y1 = hs04391_data_1x[,"ExonCoverage"]
x2 = hs04391_data_10x[,"MappedReads"]
y2 = hs04391_data_10x[,"ExonCoverage"]

plot (x1, y1, xlab="Mapped Reads", ylab="Exonic Base Coverage", ylim=c(0,70000000),
	type="l", lty=1, lwd=2, col="blue", main="Total exonic bases with increasing read depth (MIP101)") 
lines(x2,y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=total_exonic_bases, col="red", lty=2, lwd=2)
abline(h=total_mapable_exonic_bases, col="orange", lty=2, lwd=2)
lmfit = lm(y1[1900:2099]~x1[1900:2099])
b = lmfit$coefficients
fitted.y = b[1] + (b[2]*x)
lines(x = x, y = fitted.y, lty=2)

legend(locator(n = 1, type = "n"), 
	legend =c ("~210 Million Mapped Reads", "~46.2 Million Exon Bases Sequenced",
		     "Exonic Bases Observed >= 1X", "Exonic Bases Observed >= 10X", "Total EnsEMBL Bases", "Total Mapable Bases"), 
	lty=c(0,0,1,1,2,2), lwd=c(0,0,2,2,2,2), col=c("black","black","blue","dark green","red","orange")) 


#Plot CHANGE IN in exonic base coverage versus read depth
x1 = hs04391_data_1x[,"MappedReads"]
y1 = hs04391_data_1x[,"CoverageDiff"]
x2 = hs04391_data_10x[,"MappedReads"]
y2 = hs04391_data_10x[,"CoverageDiff"]

plot (x2, y2, xlab="Mapped Reads", ylab="New Exonic Bases (>= 10X depth)",
	type="l", lty=1, lwd=2, col="dark green", main="New exonic bases identified with increasing read depth (MIP101)") 
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~210 Million Mapped Reads", "~46.2 Million Exon Bases Sequenced","New Exonic Bases Observed >= 10X"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 


#Plot GENES identified versus read depth
total_genes = 32561 #Max theoretically possible
x1 = hs04391_data_1x[,"MappedReads"]
y1 = hs04391_data_1x[,"ExpressedGenes"]
x2 = hs04391_data_10x[,"MappedReads"]
y2 = hs04391_data_10x[,"ExpressedGenes"]

plot (x1, y1, xlab="Mapped Reads", ylab="Genes Detected", ylim=c(0,35000),
	type="l", lty=1, lwd=2, col="blue", main="Total genes discovered as read depth increases (MIP101)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=total_genes, col="red", lty=2, lwd=2,)

#Linear Fit - excluding the early part of the data
data_length = length(hs04391_data_1x[,1])

lmfit = lm(y1[(data_length-200):data_length]~x1[(data_length-200):data_length])
b = lmfit$coefficients
fitted.y = b[1] + (b[2]*x)
lines(x = x, y = fitted.y, lty=2)

legend(locator(n = 1, type = "n"), 
	legend =c ("~210 Million Mapped Reads", "26,325 Genes Detected","Genes Detected (>= 1X)", "Genes Detected (>= 10X)", "All EnsEMBL Genes"), 
	lty=c(0,0,1,1,2), lwd=c(0,0,2,2,2), col=c("black","black","blue","dark green","red")) 

#Plot change in GENES identified versus read depth
x2 = hs04391_data_10x[,"MappedReads"]
y2 = hs04391_data_10x[,"GenesDiff"]

plot (x2, y2, xlab="Mapped Reads", ylab="New Genes Detected",
	type="l", lty=1, lwd=2, col="blue", main="New genes discovered as read depth increases (MIP101)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~210 Million Mapped Reads", "26,325 Genes Detected","New Genes Detected (>= 10X)"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 

#Zoom in
plot (x2, y2, xlab="Mapped Reads", ylab="New Genes Detected", ylim=c(0,100),
	type="l", lty=1, lwd=2, col="blue", main="New genes discovered as read depth increases (MIP101)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~210 Million Mapped Reads", "26,325 Genes Detected","New Genes Detected (>= 10X)"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 


#HS04401
x1 = hs04401_data_1x[,"MappedReads"]
y1 = hs04401_data_1x[,"ExonCoverage"]
x2 = hs04401_data_10x[,"MappedReads"]
y2 = hs04401_data_10x[,"ExonCoverage"]

data_length = length(hs04401_data_1x[,1])

plot (x1, y1, xlab="Mapped Reads", ylab="Exonic Base Coverage", ylim=c(0,70000000),
	type="l", lty=1, lwd=2, col="blue", main="Total exonic bases with increasing read depth (MIP/5FU)") 
lines(x2,y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=total_exonic_bases, col="red", lty=2, lwd=2)
abline(h=total_mapable_exonic_bases, col="orange", lty=2, lwd=2)
lmfit = lm(y1[(data_length-200):data_length]~x1[(data_length-200):data_length])
b = lmfit$coefficients
fitted.y = b[1] + (b[2]*x)
lines(x = x, y = fitted.y, lty=2)

legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "~44.0 Million Exon Bases Sequenced",
		     "Exonic Bases Observed >= 1X", "Exonic Bases Observed >= 10X", "Total EnsEMBL Bases", "Total Mapable Bases"), 
	lty=c(0,0,1,1,2,2), lwd=c(0,0,2,2,2,2), col=c("black","black","blue","dark green","red","orange")) 


#Plot CHANGE IN in exonic base coverage versus read depth
x1 = hs04401_data_1x[,"MappedReads"]
y1 = hs04401_data_1x[,"CoverageDiff"]
x2 = hs04401_data_10x[,"MappedReads"]
y2 = hs04401_data_10x[,"CoverageDiff"]

plot (x2, y2, xlab="Mapped Reads", ylab="New Exonic Bases (>= 10X depth)", ylim=c(0,110000),
	type="l", lty=1, lwd=2, col="dark green", main="New exonic bases identified with increasing read depth (MIP/5FU)") 
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "~44.0 Million Exon Bases Sequenced","New Exonic Bases Observed >= 10X"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 

#Linear Fit - excluding the early part of the data - assumes growth in coverage has achieved linearity
lmfit = lm(y1[1900:2099]~x1[1900:2099])
b = lmfit$coefficients
fitted.y = b[1] + (b[2]*x)
lines(x = x, y = fitted.y, lty=2)

#What would x have to equal in order to get y to be the maximum number of genes
depth_needed = (total_mapable_exonic_bases - b[1])/b[2]
depth_needed

#Recreate the plot to show this
plot (x1, y1, xlab="Mapped Reads", ylab="Exonic Base Coverage", ylim=c(0,70000000), xlim=c(0,depth_needed+1000000),
	type="l", lty=1, lwd=2, col="blue", main="Total exonic bases with increasing read depth (MIP/5FU)") 
lines(x2,y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=total_exonic_bases, col="red", lty=2, lwd=2)
abline(h=total_mapable_exonic_bases, col="orange", lty=2, lwd=2)
abline(v=depth_needed, lty=2)
test_x = seq(0,depth_needed+1000000,1000000)
test_fitted.y = b[1] + (b[2]*test_x)
lines(x = test_x, y = test_fitted.y, lty=2)


legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "~44.0 Million Exon Bases Sequenced",
		     "Exonic Bases Observed >= 1X", "Exonic Bases Observed >= 10X", "Total EnsEMBL Bases", "Total Mapable Bases"), 
	lty=c(0,0,1,1,2,2), lwd=c(0,0,2,2,2,2), col=c("black","black","blue","dark green","red","orange")) 


#Plot GENES identified versus read depth
total_genes = 32561 #Max theoretically possible
x1 = hs04401_data_1x[,"MappedReads"]
y1 = hs04401_data_1x[,"ExpressedGenes"]
x2 = hs04401_data_10x[,"MappedReads"]
y2 = hs04401_data_10x[,"ExpressedGenes"]

plot (x1, y1, xlab="Mapped Reads", ylab="Genes Detected", ylim=c(0,35000),
	type="l", lty=1, lwd=2, col="blue", main="Total genes discovered as read depth increases (MIP/5FU)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=total_genes, col="red", lty=2, lwd=2,)

#Linear Fit - excluding the early part of the data
data_length = length(hs04401_data_1x[,1])

lmfit = lm(y1[(data_length-200):data_length]~x1[(data_length-200):data_length])
b = lmfit$coefficients
fitted.y = b[1] + (b[2]*x)
lines(x = x, y = fitted.y, lty=2)

legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "25,061 Genes Detected","Genes Detected (>= 1X)", "Genes Detected (>= 10X)", "All EnsEMBL Genes"), 
	lty=c(0,0,1,1,2), lwd=c(0,0,2,2,2), col=c("black","black","blue","dark green","red")) 

#Plot change in GENES identified versus read depth
x2 = hs04401_data_10x[,"MappedReads"]
y2 = hs04401_data_10x[,"GenesDiff"]

plot (x2, y2, xlab="Mapped Reads", ylab="New Genes Detected",
	type="l", lty=1, lwd=2, col="blue", main="New genes discovered as read depth increases (MIP/5FU)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "25,061 Genes Detected","New Genes Detected (>= 10X)"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 

#Zoom in
plot (x2, y2, xlab="Mapped Reads", ylab="New Genes Detected", ylim=c(0,100),
	type="l", lty=1, lwd=2, col="blue", main="New genes discovered as read depth increases (MIP/5FU)")
lines(x2, y2, type="l", lty=1, lwd=2, col="dark green")
abline(h=0, lty=2, lwd=2)
legend(locator(n = 1, type = "n"), 
	legend =c ("~113 Million Mapped Reads", "25,061 Genes Detected","New Genes Detected (>= 10X)"), 
	lty=c(0,0,1), lwd=c(0,0,2), col=c("black","black","dark green")) 







