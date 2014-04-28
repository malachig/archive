#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#HS04401

#Define a simple bash function that gets the column number for the 'Expressed' column from an input file and summarizes the number of expressed and non-expressed events from that file
function getExpStats() {
  echo
  echo $1
  COL_NUM=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')
  cut -f $COL_NUM $1 | perl -ne 'chomp $_; $total++; if ($_ == 1){$exp++;} if ($_ == 0){$not_exp++;} if (eof){$exp_p=sprintf("%.2f", ($exp/$total)*100); $not_exp_p=sprintf("%.2f", ($not_exp/$total)*100); print "\tTotal Events: $total\tExpressed: $exp ($exp_p%)\tNot Expressed: $not_exp ($not_exp_p%)\n\n"}'
}

#Define list of files - Genes, Transcripts, Exon Regions, Junctions, Boundaries, Introns, Active Intron Regions, Silent Intron Regions, Intergenics, Active Intergenic Regions, Silent Intergenic Regions
FILES="
/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_Lanes1-16_GeneExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Transcripts_v53/Summary/HS04401_Lanes1-16_TranscriptExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_Lanes1-16_ExonRegionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_Lanes1-16_JunctionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_Lanes1-16_KnownJunctionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_Lanes1-16_NovelJunctionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Boundaries_v53/Summary/HS04401_Lanes1-16_BoundaryExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Boundaries_v53/Summary/HS04401_Lanes1-16_KnownBoundaryExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Boundaries_v53/Summary/HS04401_Lanes1-16_NovelBoundaryExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Introns_v53/Summary/HS04401_Lanes1-16_IntronExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Introns_v53/Summary/HS04401_Lanes1-16_ActiveIntronRegionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Introns_v53/Summary/HS04401_Lanes1-16_SilentIntronRegionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Intergenics_v53/Summary/HS04401_Lanes1-16_IntergenicExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Intergenics_v53/Summary/HS04401_Lanes1-16_ActiveIntergenicRegionExpression_v53.txt
/projects/malachig/solexa/read_records/HS04401/Intergenics_v53/Summary/HS04401_Lanes1-16_SilentIntergenicRegionExpression_v53.txt
"

for f in $FILES
do
  getExpStats $f
done




