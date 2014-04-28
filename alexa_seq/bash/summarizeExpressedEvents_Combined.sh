#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Get the number of features of each type that are expressed in BOTH HS04391 and HS04401

function getCombinedExpressedFeatures() {
  echo $2
  COL_NUM=$(cat $1 | head -n 1 | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Expressed/){print "$c"}}')
  cat $1 | cut -f 1,$COL_NUM | perl -ne '@line = split("\t", $_); if ($line[1]==1){print "$_";}' | sort | uniq | wc -l
  echo
}

getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/ENST_v53/Summary/*_GeneExpression_v53.txt" "Gene"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Transcripts_v53/Summary/*_TranscriptExpression_v53.txt" "Transcript"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/ENST_v53/Summary/*_ExonRegionExpression_v53.txt" "ExonRegion"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Junctions_v53/Summary/*_JunctionExpression_v53.txt" "Junction"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Junctions_v53/Summary/*_KnownJunctionExpression_v53.txt" "KnownJunction"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Junctions_v53/Summary/*_NovelJunctionExpression_v53.txt" "NovelJunction"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Boundaries_v53/Summary/*_BoundaryExpression_v53.txt" "Boundary"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Boundaries_v53/Summary/*_KnownBoundaryExpression_v53.txt" "KnownBoundary"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Boundaries_v53/Summary/*_NovelBoundaryExpression_v53.txt" "NovelBoundary"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Introns_v53/Summary/*_IntronExpression_v53.txt" "Intron"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Introns_v53/Summary/*_ActiveIntronRegionExpression_v53.txt" "ActiveIntron"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Introns_v53/Summary/*_SilentIntronRegionExpression_v53.txt" "SilentIntron"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Intergenics_v53/Summary/*_IntergenicExpression_v53.txt" "Intergenic"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Intergenics_v53/Summary/*_ActiveIntergenicRegionExpression_v53.txt" "ActiveIntergenic"
getCombinedExpressedFeatures "/projects/malachig/solexa/read_records/*/Intergenics_v53/Summary/*_SilentIntergenicRegionExpression_v53.txt" "SilentIntergenic"



