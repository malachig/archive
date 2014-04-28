#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Move to the root directory
cd /projects/malachig/solexa/figures_and_stats

#Define a function to grab a column of feature IDs from a DE feature list
function get_DE_feature_ids(){
  COL_NUM1=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Gene_Name/){print "$c"}}')
  COL_NUM2=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/EnsEMBL_Gene_ID/){print "$c"}}')
  cut -f 1,$COL_NUM1,$COL_NUM2 $1 | perl -ne 'chomp $_; if ($_ =~ /\d+/){print "$_\n";}'
}

#Define a function to grab a column of feature IDs from an SI feature list and filter it so that only features with a percent feature contribution (PFC) are greater 75
function get_SI_feature_ids(){
  COL_NUM1=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/EnsEMBL_Gene_ID/){print "$c"}}')
  COL_NUM2=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/percent_SEQ_Log2_DE/){print "$c"}}')
  cut -f 1,$COL_NUM1,$COL_NUM2 $1 | perl -ne 'chomp $_; if ($_ =~ /(\S+)\s+(\S+)\s+(\S+)/){if ($3 >= 50){print "$1\t$2\n";}}'
}

#Define some temp files
TMP_DE="de_tmp.txt"
TMP_SI="si_tmp.txt"

mkdir "DE_SI_intersect/"

function getIntersection() {
  echo
  echo $1 $2 $3
  COUNT1=$(wc -l $2 | perl -ne 'if ($_ =~ /^(\d+)/){print "$1";}')
  COUNT2=$(wc -l $3 | perl -ne 'if ($_ =~ /^(\d+)/){print "$1";}')
  echo Line counts: $COUNT1 $COUNT2

  #Define an output file for the intersection list
  RESULT_FILE=$(echo "DE_SI_intersect/"$1"_intersection.txt")

  get_DE_feature_ids $2 | sort > $TMP_DE
  get_SI_feature_ids $3 | sort > $TMP_SI

  INT=$(join $TMP_DE $TMP_SI | wc -l)
  echo Intersection between the DE and SI features of type $TYPE = $INT
  join $TMP_DE $TMP_SI > $RESULT_FILE
}

getIntersection "Transcript" DE/Transcripts_v53/Mip101_vs_Mip5FuR_Transcript_DE_Values_Significant_MTC.txt SI/Transcripts_v53/Mip101_vs_Mip5FuR_Transcript_SI_Values_Sorted_Cutoff.txt
getIntersection "ExonRegion" DE/ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_DE_Values_Significant_MTC.txt SI/ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_SI_Values_Sorted_Cutoff.txt
getIntersection "Junction" DE/Junctions_v53/Mip101_vs_Mip5FuR_Junction_DE_Values_Significant_MTC.txt SI/Junctions_v53/Mip101_vs_Mip5FuR_Junction_SI_Values_Sorted_Cutoff.txt
getIntersection "KnownJunction" DE/Junctions_v53/Mip101_vs_Mip5FuR_KnownJunction_DE_Values_Significant_MTC.txt SI/Junctions_v53/Mip101_vs_Mip5FuR_KnownJunction_SI_Values_Sorted_Cutoff.txt
getIntersection "NovelJunction" DE/Junctions_v53/Mip101_vs_Mip5FuR_NovelJunction_DE_Values_Significant_MTC.txt SI/Junctions_v53/Mip101_vs_Mip5FuR_NovelJunction_SI_Values_Sorted_Cutoff.txt
getIntersection "Boundary" DE/Boundaries_v53/Mip101_vs_Mip5FuR_Boundary_DE_Values_Significant_MTC.txt SI/Boundaries_v53/Mip101_vs_Mip5FuR_Boundary_SI_Values_Sorted_Cutoff.txt 
getIntersection "KnownBoundary" DE/Boundaries_v53/Mip101_vs_Mip5FuR_KnownBoundary_DE_Values_Significant_MTC.txt SI/Boundaries_v53/Mip101_vs_Mip5FuR_KnownBoundary_SI_Values_Sorted_Cutoff.txt 
getIntersection "NovelBoundary" DE/Boundaries_v53/Mip101_vs_Mip5FuR_NovelBoundary_DE_Values_Significant_MTC.txt SI/Boundaries_v53/Mip101_vs_Mip5FuR_NovelBoundary_SI_Values_Sorted_Cutoff.txt 
getIntersection "Intron" DE/Introns_v53/Mip101_vs_Mip5FuR_Intron_DE_Values_Significant_MTC.txt SI/Introns_v53/Mip101_vs_Mip5FuR_Intron_SI_Values_Sorted_Cutoff.txt  
getIntersection "SilentIntronRegion" DE/Introns_v53/Mip101_vs_Mip5FuR_SilentIntronRegion_DE_Values_Significant_MTC.txt SI/Introns_v53/Mip101_vs_Mip5FuR_SilentIntronRegion_SI_Values_Sorted_Cutoff.txt
getIntersection "ActiveIntronRegion" DE/Introns_v53/Mip101_vs_Mip5FuR_ActiveIntronRegion_DE_Values_Significant_MTC.txt SI/Introns_v53/Mip101_vs_Mip5FuR_ActiveIntronRegion_SI_Values_Sorted_Cutoff.txt

rm -f $TMP_DE $TMP_SI

