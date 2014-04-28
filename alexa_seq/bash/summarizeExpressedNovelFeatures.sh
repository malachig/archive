#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Summarize the number of expressed novel features of the type: exonJunction and exonBoundary
#Calculate the following, combining both MIP101 and MIP/5FU libraries

#ALL
EJ_FILE_A="/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_JunctionExpression_v53.txt"
EJ_FILE_B="/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_Lanes1-16_JunctionExpression_v53.txt"
EB_FILE_A="/projects/malachig/solexa/read_records/HS04391/Boundaries_v53/Summary/HS04391_Lanes1-23_BoundaryExpression_v53.txt"
EB_FILE_B="/projects/malachig/solexa/read_records/HS04401/Boundaries_v53/Summary/HS04401_Lanes1-16_BoundaryExpression_v53.txt"

#NOVEL
NEJ_FILE_A="/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_NovelJunctionExpression_v53.txt"
NEJ_FILE_B="/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_Lanes1-16_NovelJunctionExpression_v53.txt"
NEB_FILE_A="/projects/malachig/solexa/read_records/HS04391/Boundaries_v53/Summary/HS04391_Lanes1-23_NovelBoundaryExpression_v53.txt"
NEB_FILE_B="/projects/malachig/solexa/read_records/HS04401/Boundaries_v53/Summary/HS04401_Lanes1-16_NovelBoundaryExpression_v53.txt"

#A.) Totals
#Total number of each feature type: FeatureCount
#Total number of each feature type expressed: ExpressedFeatureCount
function get_expressed_features(){
  #echo ARGS: $1 $2 $3
  TYPE=$1
  FEATURE_COUNT=$(grep -v Expressed $2 | wc -l )
  COL_EXP=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')
  EXPRESSED_FEATURE_COUNT=$(cat $2 $3 | cut -f 1,$COL_EXP | perl -ne 'chomp($_); @line=split("\t", $_); if ($line[1]==1){print "$line[0]\n";}' | sort | uniq | wc -l)
  echo
  echo "TYPE = $TYPE   FEATURE_COUNT = $FEATURE_COUNT   EXPRESSED_FEATURE_COUNT = $EXPRESSED_FEATURE_COUNT" 
  echo
}
get_expressed_features "Junction" $EJ_FILE_A $EJ_FILE_B
get_expressed_features "Boundary" $EB_FILE_A $EB_FILE_B
get_expressed_features "NovelJunction" $NEJ_FILE_A $NEJ_FILE_B
get_expressed_features "NovelBoundary" $NEB_FILE_A $NEB_FILE_B


#B.) Percent mRNA/EST supported
#Percent of FeatureCount that are mRNA/EST supported
#Percent of ExpressedFeatureCount that are mRNA/EST supported
#Is there a statistically significant enrichment for those that are expressed versus the total?
function get_expressed_supported_features(){
  #echo ARGS: $1 $2 $3
  TYPE=$1
  COL_MRNA=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_mRNA_Count/){print "$c"}}')
  COL_EST=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EST_Count/){print "$c"}}')
  COL_EXP=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')

  FEATURE_COUNT=$(cat $2 | cut -f 1,$COL_MRNA,$COL_EST | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0 || $line[2]>0)){print "$line[0]\n";}' | wc -l)
  EXPRESSED_FEATURE_COUNT=$(cat $2 $3 | cut -f 1,$COL_MRNA,$COL_EST,$COL_EXP | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0 || $line[2]>0) && $line[3]==1){print "$line[0]\n";}' | sort | uniq | wc -l)
  echo
  echo "TYPE = $TYPE   FEATURE_COUNT_SEQ_SUPPORTED = $FEATURE_COUNT   EXPRESSED_FEATURE_COUNT_SEQ_SUPPORTED = $EXPRESSED_FEATURE_COUNT" 
  echo
}
get_expressed_supported_features "Junction" $EJ_FILE_A $EJ_FILE_B
get_expressed_supported_features "Boundary" $EB_FILE_A $EB_FILE_B
get_expressed_supported_features "NovelJunction" $NEJ_FILE_A $NEJ_FILE_B
get_expressed_supported_features "NovelBoundary" $NEB_FILE_A $NEB_FILE_B


#C.) Percent conserved 
#Percent of FeatureCount that are Conserved
#Percent of ExpressedFeatureCount that are Conserved
#Is there a statistically significant enrichment for those that are expressed versus the total?
function get_expressed_conserved_features(){
  #echo ARGS: $1 $2 $3
  TYPE=$1
  COL_CON=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Conserved_Species_Count/){print "$c"}}')
  COL_EXP=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')

  FEATURE_COUNT=$(cat $2 | cut -f 1,$COL_CON | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0)){print "$line[0]\n";}' | wc -l)
  EXPRESSED_FEATURE_COUNT=$(cat $2 $3 | cut -f 1,$COL_CON,$COL_EXP | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0) && $line[2]==1){print "$line[0]\n";}' | sort | uniq | wc -l)
  echo
  echo "TYPE = $TYPE   FEATURE_COUNT_CONSERVED = $FEATURE_COUNT   EXPRESSED_FEATURE_COUNT_CONSERVED = $EXPRESSED_FEATURE_COUNT" 
  echo
}
get_expressed_conserved_features "Junction" $EJ_FILE_A $EJ_FILE_B
get_expressed_conserved_features "Boundary" $EB_FILE_A $EB_FILE_B
get_expressed_conserved_features "NovelJunction" $NEJ_FILE_A $NEJ_FILE_B
get_expressed_conserved_features "NovelBoundary" $NEB_FILE_A $NEB_FILE_B



#D.) Percent affecting protein coding
#Percent of FeatureCount that affect protein coding content
#Percent of ExpressedFeatureCount that affect protein coding content
#Is there a statistically significant enrichment for those that are expressed versus the total?
function get_expressed_coding_features(){
  #echo ARGS: $1 $2 $3
  TYPE=$1
  COL_COD=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Coding_Base_Count/){print "$c"}}')
  COL_EXP=$(head -n 1 $2 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')

  FEATURE_COUNT=$(cat $2 | cut -f 1,$COL_COD | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0)){print "$line[0]\n";}' | wc -l)
  EXPRESSED_FEATURE_COUNT=$(cat $2 $3 | cut -f 1,$COL_COD,$COL_EXP | perl -ne 'chomp($_); @line=split("\t", $_); if (($line[1]>0) && $line[2]==1){print "$line[0]\n";}' | sort | uniq | wc -l)
  echo
  echo "TYPE = $TYPE   FEATURE_COUNT_CODING = $FEATURE_COUNT   EXPRESSED_FEATURE_COUNT_CODING = $EXPRESSED_FEATURE_COUNT" 
  echo
}
get_expressed_coding_features "Junction" $EJ_FILE_A $EJ_FILE_B
get_expressed_coding_features "Boundary" $EB_FILE_A $EB_FILE_B
get_expressed_coding_features "NovelJunction" $NEJ_FILE_A $NEJ_FILE_B
get_expressed_coding_features "NovelBoundary" $NEB_FILE_A $NEB_FILE_B








