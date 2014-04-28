#!/bin/bash
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.


E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /projects/malachig/solexa 53 HS0499  (i.e. base_dir ensembl_version library_id)"
  exit $E_BADARGS
fi  
BASE_DIR=$1
VERSION=$2
LIB=$3

#Print summary of args 
echo
echo "Processing: $1 $2 $3"

#Define the output file
OUT=$(echo $BASE_DIR/figures_and_stats/$LIB/Expression_v$VERSION/$LIB\_ExpressedFeatures.txt)

#Define a simple bash function that gets the column number for the 'Expressed' column from an input file and summarizes the number of expressed and non-expressed events from that file
echo "FeatureType	TotalEvents	Expressed	Expressed_Percent	NotExpressed	NotExpressed_Percent" > $OUT
function getExpStats() {
  FEATURE_TYPE=$2
  export FEATURE_TYPE
  COL_NUM=$(head -n 1 $1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')
  cut -f $COL_NUM $1 | perl -ne '$type=$ENV{"FEATURE_TYPE"}; chomp $_; $total++; if ($_ == 1){$exp++;} if ($_ == 0){$not_exp++;} if (eof){$exp_p=sprintf("%.2f", ($exp/$total)*100); $not_exp_p=sprintf("%.2f", ($not_exp/$total)*100); print "$type\t$total\t$exp\t$exp_p\t$not_exp\t$not_exp_p\n"}'
}

#Define list of files - Genes, Transcripts, Exon Regions, Junctions, Boundaries, Introns, Active Intron Regions, Silent Intron Regions, Intergenics, Active Intergenic Regions, Silent Intergenic Regions
getExpStats $BASE_DIR/read_records/$LIB/ENST_v$VERSION/Summary/$LIB\_GeneExpression_v$VERSION.txt Gene >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Transcripts_v$VERSION/Summary/$LIB\_TranscriptExpression_v$VERSION.txt Transcript >> $OUT 
getExpStats $BASE_DIR/read_records/$LIB/ENST_v$VERSION/Summary/$LIB\_ExonRegionExpression_v$VERSION.txt ExonRegion >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Junctions_v$VERSION/Summary/$LIB\_JunctionExpression_v$VERSION.txt Junction >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Junctions_v$VERSION/Summary/$LIB\_KnownJunctionExpression_v$VERSION.txt KnownJunction >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Junctions_v$VERSION/Summary/$LIB\_NovelJunctionExpression_v$VERSION.txt NovelJunction >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Boundaries_v$VERSION/Summary/$LIB\_BoundaryExpression_v$VERSION.txt Boundary >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Boundaries_v$VERSION/Summary/$LIB\_KnownBoundaryExpression_v$VERSION.txt KnownBoundary >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Boundaries_v$VERSION/Summary/$LIB\_NovelBoundaryExpression_v$VERSION.txt NovelBoundary >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Introns_v$VERSION/Summary/$LIB\_IntronExpression_v$VERSION.txt Intron >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Introns_v$VERSION/Summary/$LIB\_ActiveIntronRegionExpression_v$VERSION.txt ActiveIntronRegion >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Introns_v$VERSION/Summary/$LIB\_SilentIntronRegionExpression_v$VERSION.txt SilentIntronRegion >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Intergenics_v$VERSION/Summary/$LIB\_IntergenicExpression_v$VERSION.txt Intergenic >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Intergenics_v$VERSION/Summary/$LIB\_ActiveIntergenicRegionExpression_v$VERSION.txt ActiveIntergenic >> $OUT
getExpStats $BASE_DIR/read_records/$LIB/Intergenics_v$VERSION/Summary/$LIB\_SilentIntergenicRegionExpression_v$VERSION.txt SilentIntergenic >> $OUT

echo "Done"
echo
