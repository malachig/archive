#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#HS04391
#move to a working directory
#cd /projects/malachig/solexa/figures_and_stats/HS04391/temp
#GENE_EXP="/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_GeneExpression_v53.txt"
#NAME="MIP101"

#HS04401
#move to a working directory
cd /projects/malachig/solexa/figures_and_stats/HS04401/temp
GENE_EXP="/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_Lanes1-16_GeneExpression_v53.txt"
NAME="MIP5FU"


#Define the alexa database to use
DB_VERSION="ALEXA_hs_53_36o"

#Define the gene expression file of interest

#Define the two temp files to be use
TEMP_EXP="expressed_ids.txt"
TEMP_TYPE="type_ids.txt"
SUMMARY="ExpressedGenesByType.txt"

#Clear output file and add header
rm $SUMMARY
echo "Type	TypeCount	ExpressedGeneCount	Intersection	PercentOfTypeIsExpressed	PercentOfExpressedIsType" >> $SUMMARY

#Define a function to get gene ids for all expressed genes from the input file using the 'Expressed' column
function getExpCol() {
  echo
  echo $1
  COL_NUM=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')
  cut -f 1,$COL_NUM $1 | perl -ne 'chomp $_; $total++; if ($_ =~ /(\S+)\s+1/){print "$1\n"}'
}
getExpCol $GENE_EXP | sort > $TEMP_EXP

#Get all possible types from the mysql database
TYPES=$(echo 'SELECT distinct gene_type FROM Gene;'  | mysql -N $DB_VERSION)

#TYPES="protein_coding pseudogene rRNA scRNA miRNA"

#Now for each type, summarize the percentages
for TYPE in $TYPES
do

echo 'SELECT id FROM Gene WHERE gene_type="'$TYPE'";' | mysql -N $DB_VERSION | sort > $TEMP_TYPE

EXPRESSED_AND_TYPE_COUNT=$(join $TEMP_EXP $TEMP_TYPE | wc -l)
EXPRESSED_COUNT=$(cat $TEMP_EXP | wc -l)
TYPE_COUNT=$(cat $TEMP_TYPE | wc -l)

echo
echo $TYPE
echo "Number of genes expressed in "$NAME "=" $EXPRESSED_COUNT
echo "Number of protein coding genes in " $DB_VERSION "=" $TYPE_COUNT
echo "Intersection of these two categories = " $EXPRESSED_AND_TYPE_COUNT
echo


# What percentage of the expressed genes are protein coding?
RESULT1=$(echo "($EXPRESSED_AND_TYPE_COUNT / $EXPRESSED_COUNT)*100" | bc -l)
RESULT1_R=$(printf '%.3f' $RESULT1)
echo "Percentage of genes that were expressed that are this type of gene = " $RESULT1_R

# What percentage of all protein coding genes were found to be expressed?
RESULT2=$(echo "($EXPRESSED_AND_TYPE_COUNT / $TYPE_COUNT)*100" | bc -l)
RESULT2_R=$(printf '%.3f' $RESULT2)
echo "Percentage of all genes of this type that were found to be expressed = " $RESULT2_R
echo

#Add results to the output file
echo $TYPE"	"$TYPE_COUNT"	"$EXPRESSED_COUNT"	"$EXPRESSED_AND_TYPE_COUNT"	"$RESULT2_R"	"$RESULT1_R >> $SUMMARY

done

#Clean up
rm $TEMP_EXP
rm $TEMP_TYPE
cp $SUMMARY ..

