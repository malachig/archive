#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#For an alternative expression database.  Determine the number of features that are EnsEMBL supported, mRNA/EST supported and totally novel
#Genes and transcripts are by definition EnSEMBL supported
#Active intron and intergenic regions are by definition sequence supported but not EnsEMBL
#Silent intron and intergenic regions are by definition totally novel

cd /projects/malachig/sequence_databases/hs_53_36o/

#Genes
#Dont count these - only features of them.

#Transcripts
#Dont count these - only features of them

#Junctions
JUNC="exonJunctions/exonJunctions_62mers_annotated.txt.gz"
ENSG_COL=$(zcat $JUNC | head -n 1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EnsEMBL_Count/){print "$c"}}')
MRNA_COL=$(zcat $JUNC | head -n 1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_mRNA_Count/){print "$c"}}')
EST_COL=$(zcat $JUNC | head -n 1  | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EST_Count/){print "$c"}}')

zcat $JUNC | cut -f $ENSG_COL,$MRNA_COL,$EST_COL | perl -ne 'chomp $_; if ($_ =~ /(\d+)\s+(\d+)\s+(\d+)/){$count++; if ($1 >0){$ensg++;}elsif($2 > 0 || $3 > 0){$seq++;}else{$novel++}} if(eof){print "\n\nJUNCTIONS.  count: $count ensg:$ensg seq:$seq novel:$novel"}'

#Boundaries
BOUND="exonBoundaries/exonBoundaries_62mers_annotated.txt.gz"
ENSG_COL=$(zcat $BOUND | head -n 1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EnsEMBL_Count/){print "$c"}}')
MRNA_COL=$(zcat $BOUND | head -n 1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_mRNA_Count/){print "$c"}}')
EST_COL=$(zcat $BOUND | head -n 1  | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EST_Count/){print "$c"}}')

zcat $BOUND | cut -f $ENSG_COL,$MRNA_COL,$EST_COL | perl -ne 'chomp $_; if ($_ =~ /(\d+)\s+(\d+)\s+(\d+)/){$count++; if ($1 >0){$ensg++;}elsif($2 > 0 || $3 > 0){$seq++;}else{$novel++}} if(eof){print "\n\nBOUNDARIES.  count: $count ensg:$ensg seq:$seq novel:$novel"}'





