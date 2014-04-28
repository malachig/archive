#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Import SNP data
datadir = "/projects/malachig/solexa/maq_mapping_files/13288AAXX_Lanes1-2_HS03271/"
dir(datadir)
setwd(datadir)

snp_data = read.table("both_lanes_SNPs.txt", header=F, quote="",sep="\t", comment.char="", as.is=c(1,3,4), na.strings='', nrows=820433)
names(snp_data) = c("chromosome", "position", "ref_base", "consensus_base", "consensus_quality", "read_depth", "avg_number_hits", "highest_map_qual", "allele_quality_diff")

#How many SNPS are there from alignment regions with a depth of at least 6X and a consensus quality >30
length(which(snp_data[,"consensus_quality"] >= 30 & snp_data[,"read_depth"] >= 6))  #2746

quality_snps_i = which(snp_data[,"consensus_quality"] >= 30 & snp_data[,"read_depth"] >= 6)
quality_snps = snp_data[quality_snps_i,]

quality_snps[,c("chromosome", "position", "ref_base", "consensus_base", "read_depth")]
