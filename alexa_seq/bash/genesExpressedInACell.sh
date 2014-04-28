#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.


#Estimate the number of genes/transcripts in a single cell

#Make some simple assumptions:
#1.) A cell is a self-contained unit.  If it needs a gene, it must express it
#2.) The minimum copy number gene is 1 copy per cell
#3.) Our data has identified all truly expressed genes
#4.) Our data has provides an accurate estimate of the relative dynamic range of expression levels (from the most lowly to the most highly expressed gene)

#Get expression values for all expressed genes - > 12396 genes expressed
#Use RAW values for this estimation
cd /projects/malachig/solexa/figures_and_stats/temp/
cut -f 2,14,21 /projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_GeneExpression_v53.txt | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[2]==1){print "$line[0]\t$line[1]\n";}' > GeneExpressionLevels_RAW.txt

#Get the minimum expression value for an expressed gene: 7.5897035881 (copy number = 1)
sort -n -k 2  GeneExpressionLevels_RAW.txt | head -n 1

#Get the maximum expression value for an expressed gene: 47372.7989623865 (copy number = 6242)
sort -n -k 2  GeneExpressionLevels_RAW.txt | tail -n 1

#Set this minimum to be equivalent to a copy number of '1 per cell'
perl -ne 'chomp($_); @line=split("\t", $_); $copy_number=$line[1]/7.5897035881; print "$_\t$copy_number\n"' GeneExpressionLevels_RAW.txt > GeneCopyNumbers.txt

#Get the grand total copy number -> 460019.371471546
perl -ne 'chomp($_); @line=split("\t", $_); $copy_number_count+=$line[2]; if (eof){print"\n\nGRAND COPY NUMBER: $copy_number_count\n\n"}' GeneCopyNumbers.txt




