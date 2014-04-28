#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Using exon regions and junctions, determine the number of genes where multiple expressed transcripts are observed

#Get input files
#HS04391

cd /projects/malachig/solexa/figures_and_stats/HS04391/

GENES="/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_GeneExpression_v53.txt"
EXPRESSED_GENES="ExpressedGenes.txt"
ENSG_COL=$(head -n 1 $GENES | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
EXP_COL=$(head -n 1 $GENES | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Expressed"){print "$count";}}')
cut -f $ENSG_COL,$EXP_COL $GENES | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[1] == 1){print "$line[0]\n"}' > $EXPRESSED_GENES

#For how many genes could we theoretically detect multiple transcripts expressed
TRANS="/projects/malachig/solexa/read_records/HS04391/Transcripts_v53/Summary/HS04391_Lanes1-23_TranscriptGeneSummary_v53.txt"
cut -f 2,8,9,10 $TRANS | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[2] > 1){print "$_\n"}' | wc -l

#How many of these show evidence for expression of at least one transcript?
cut -f 2,8,9,10 $TRANS | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[2] > 1 && $line[3] > 0){print "$_\n"}' | wc -l

#How many of these show evidence for expression of multiple transcripts?
cut -f 2,8,9,10 $TRANS | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[2] > 1 && $line[3] > 1){print "$_\n"}' | wc -l

#How many of these would be dropped if we required that they have evidence for expression at the gene level? -> very few 
grep -f ExpressedGenes.txt GenesWithMutlipleKnownTranscriptsExpressed.txt | wc -l

#Store the list of genes with multiple transcripts
cut -f 2,8,9,10 $TRANS | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[2] > 1 && $line[3] > 1){print "$line[0]\n"}' > "GenesWithMutlipleKnownTranscriptsExpressed.txt"

#Now what about novel splicing events? ....
JUNC="/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_NovelJunctionExpression_v53.txt"
ENSG_COL=$(head -n 1 $JUNC | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
EXP_COL=$(head -n 1 $JUNC | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Expressed"){print "$count";}}')
cut -f $ENSG_COL,$EXP_COL $JUNC | perl -ne 'if ($_ =~ /(\S+)\s+1/){print "$1\n"}' | sort | uniq > "tmp.txt"
join tmp.txt ExpressedGenes.txt > "GenesWithNovelJunctions.txt"

BOUND="/projects/malachig/solexa/read_records/HS04391/Boundaries_v53/Summary/HS04391_Lanes1-23_NovelBoundaryExpression_v53.txt"
ENSG_COL=$(head -n 1 $BOUND | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
EXP_COL=$(head -n 1 $BOUND | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Expressed"){print "$count";}}')
cut -f $ENSG_COL,$EXP_COL $BOUND | perl -ne 'if ($_ =~ /(\S+)\s+1/){print "$1\n"}' | sort | uniq > "tmp.txt"
join tmp.txt ExpressedGenes.txt > "GenesWithNovelBoundaries.txt"

INTRON_A="/projects/malachig/solexa/read_records/HS04391/Introns_v53/Summary/HS04391_Lanes1-23_ActiveIntronRegionExpression_v53.txt"
ENSG_COL=$(head -n 1 $INTRON_A | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
EXP_COL=$(head -n 1 $INTRON_A | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Expressed"){print "$count";}}')
cut -f $ENSG_COL,$EXP_COL $INTRON_A | perl -ne 'if ($_ =~ /(\S+)\s+1/){print "$1\n"}' | sort | uniq > "tmp.txt"
join tmp.txt ExpressedGenes.txt > "GenesWithNovelIntronActiveRegions.txt"

INTRON_S="/projects/malachig/solexa/read_records/HS04391/Introns_v53/Summary/HS04391_Lanes1-23_SilentIntronRegionExpression_v53.txt"
ENSG_COL=$(head -n 1 $INTRON_S | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
EXP_COL=$(head -n 1 $INTRON_S | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Expressed"){print "$count";}}')
cut -f $ENSG_COL,$EXP_COL $INTRON_S | perl -ne 'if ($_ =~ /(\S+)\s+1/){print "$1\n"}' | sort | uniq > "tmp.txt"
join tmp.txt ExpressedGenes.txt > "GenesWithNovelIntronSilentRegions.txt"

rm -f tmp.txt


#Summarize the counts
wc -l GenesWithMutlipleKnownTranscriptsExpressed.txt
cat GenesWithMutlipleKnownTranscriptsExpressed.txt GenesWithNovelJunctions.txt | sort | uniq | wc -l
cat GenesWithMutlipleKnownTranscriptsExpressed.txt GenesWithNovelJunctions.txt GenesWithNovelBoundaries.txt | sort | uniq | wc -l
cat GenesWithMutlipleKnownTranscriptsExpressed.txt GenesWithNovelJunctions.txt GenesWithNovelBoundaries.txt GenesWithNovelIntronActiveRegions.txt | sort | uniq | wc -l
cat GenesWithMutlipleKnownTranscriptsExpressed.txt GenesWithNovelJunctions.txt GenesWithNovelBoundaries.txt GenesWithNovelIntronActiveRegions.txt GenesWithNovelIntronSilentRegions.txt | sort | uniq | wc -l




