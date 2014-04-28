#!/bin/bash
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

# This script counts the jobs on the cluster for a user every N seconds
E_BADARGS=65
if [ ! $# == 5 ]
then
  echo
  echo "Usage: `basename $0` GENE_NAME LibA_Name LibB_Name Data_Path Ensembl_Version"
  echo "Example: getExampleGeneData.sh UMPS HS04391 HS04401 /projects/malachig/solexa/read_records/ 53"
  echo
  exit $E_BADARGS
fi  

GENE=$1
LIBA=$2
LIBB=$3
BASEDIR=$4
VERSION=$5
echo
echo Searching for exon/junction data for Gene:$GENE in libraries $LIBA and $LIBB in the directory $BASEDIR Ensembl version: $VERSION
echo 

J_FILE_A=$(echo $BASEDIR$LIBA/Junctions\_v$VERSION/Summary/$LIBA\_JunctionExpression\_v$VERSION.txt)
E_FILE_A=$(echo $BASEDIR$LIBA/ENST\_v$VERSION/Summary/$LIBA\_ExonRegionExpression\_v$VERSION.txt)
J_FILE_B=$(echo $BASEDIR$LIBB/Junctions\_v$VERSION/Summary/$LIBB\_JunctionExpression\_v$VERSION.txt)
E_FILE_B=$(echo $BASEDIR$LIBB/ENST\_v$VERSION/Summary/$LIBB\_ExonRegionExpression\_v$VERSION.txt)
RESULT_FILE=$(echo $GENE.txt)

echo processing: $E_FILE_A
cat $E_FILE_A | cut -f 4,9,10,20,25,27,32 | grep $GENE > e_file_a.txt
echo processing: $E_FILE_B
cat $E_FILE_B | cut -f 4,25,27,32 | grep $GENE | cut -f 2,3,4 > e_file_b.txt
paste e_file_a.txt e_file_b.txt > e_file_ab.txt

echo processing: $J_FILE_A
cat $J_FILE_A | cut -f 4,12,13,25,31,33,38 | grep $GENE > j_file_a.txt
echo processing: $J_FILE_B
cat $J_FILE_B | cut -f 4,31,33,38 | grep $GENE | cut -f 2,3,4 > j_file_b.txt
paste j_file_a.txt j_file_b.txt > j_file_ab.txt

#Sort and remove extra data columns
echo "StartCoord	EndCoord	FeatureName	ExpressionLibA	ExpressionLibB" > header.txt
cat e_file_ab.txt j_file_ab.txt | sort -n -k 2 | cut -f 2,3,4,5,8 > sorted.txt
cat header.txt sorted.txt > $RESULT_FILE

rm -f e_file_a.txt e_file_b.txt e_file_ab.txt j_file_a.txt j_file_b.txt j_file_ab.txt sorted.txt header.txt




